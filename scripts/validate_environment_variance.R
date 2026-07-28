#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("project-root", "input-dir", "output-dir", "figure-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

as_absolute <- function(path) {
  if (grepl("^/", path)) path else file.path(getwd(), path)
}
project_root <- normalizePath(as_absolute(args[["project-root"]]), mustWork = TRUE)
input_dir <- normalizePath(as_absolute(args[["input-dir"]]), mustWork = TRUE)
output_dir <- as_absolute(args[["output-dir"]])
figure_dir <- as_absolute(args[["figure-dir"]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "digest", "ggplot2", "jsonlite", "permute", "ragg",
  "readr", "scales", "svglite", "systemfonts", "vegan"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

source(file.path(project_root, "R", "theme_pub.R"), local = FALSE)
suppressPackageStartupMessages(library(vegan))

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
  x <- gsub(input_dir, "<INPUT_DIR>", x, fixed = TRUE)
  x <- gsub(output_dir, "<OUTPUT_DIR>", x, fixed = TRUE)
  x <- gsub(figure_dir, "<FIGURE_DIR>", x, fixed = TRUE)
  x <- gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
  x
}

write_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(output_dir, filename), na = "")
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character(),
    col_types = readr::cols(.default = readr::col_character())
  )
  ids <- as.character(x[[1L]])
  if (anyDuplicated(ids)) {
    stop("Duplicated identifiers in ", basename(path), call. = FALSE)
  }
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

check_rows <- list()
add_check <- function(id, category, observed, expected, passed) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check_id = id,
    category = category,
    observed = paste(observed, collapse = ","),
    expected = paste(expected, collapse = ","),
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    check.names = FALSE
  )
}

near <- function(observed, expected, tolerance = 1e-8) {
  observed <- as.numeric(observed)
  expected <- as.numeric(expected)
  length(observed) == length(expected) &&
    all(is.finite(observed)) &&
    all(is.finite(expected)) &&
    max(abs(observed - expected)) <= tolerance
}

command_path <- function(command) {
  path <- Sys.which(command)
  if (!nzchar(path)) "" else normalizePath(path, mustWork = TRUE)
}

run_command <- function(command, args = character()) {
  status <- 0L
  output <- tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      status <<- 1L
      conditionMessage(e)
    }
  )
  command_status <- attr(output, "status")
  if (!is.null(command_status)) status <- as.integer(command_status)
  list(status = status, output = output)
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p <= 0.001, "<=0.001", sprintf("%.3f", p)))
}

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  environment = file.path(input_dir, "environment.tsv"),
  source_summary = file.path(input_dir, "source_summary.json"),
  prepare_script = file.path(project_root, "scripts", "prepare_pilot_data.R"),
  theme_pub = file.path(project_root, "R", "theme_pub.R")
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  environment = "45814f6724e8abd901b8063ba227c42716777cd0749ac130e9cae99f0477f7a8",
  source_summary = "acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b",
  prepare_script = "fb9ce168b70aeb3636ab73c3c3fabf4208934a04f8eb89a547dc93c7d66a1e3f",
  theme_pub = "8d3821a485aeb529b12184bbe977c1ca952131a22d5eb69b802088ce91909beb"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 25 input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

observed_sha256 <- vapply(input_paths, sha256_file, character(1))
input_audit <- data.frame(
  asset = names(input_paths),
  path = vapply(
    input_paths,
    function(path) sanitize_text(normalizePath(path, mustWork = TRUE)),
    character(1)
  ),
  bytes = unname(file.info(input_paths)$size),
  observed_sha256 = unname(observed_sha256),
  expected_sha256 = unname(expected_sha256[names(input_paths)]),
  status = ifelse(
    observed_sha256 == expected_sha256[names(input_paths)],
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)
for (i in seq_len(nrow(input_audit))) {
  add_check(
    paste0("sha256-", input_audit$asset[[i]]),
    "input",
    input_audit$observed_sha256[[i]],
    input_audit$expected_sha256[[i]],
    input_audit$status[[i]] == "PASS"
  )
}

expected_versions <- c(
  digest = "0.6.36",
  ggplot2 = "3.5.2",
  jsonlite = "1.8.8",
  permute = "0.9-7",
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
  vegan = "2.6-6.1"
)
observed_versions <- vapply(
  names(expected_versions),
  function(package) packageDescription(package, fields = "Version"),
  character(1)
)
package_audit <- data.frame(
  package = names(expected_versions),
  expected_version = unname(expected_versions),
  observed_version = unname(observed_versions),
  status = ifelse(observed_versions == expected_versions, "PASS", "FAIL"),
  check.names = FALSE
)
r_version <- paste(R.version$major, R.version$minor, sep = ".")
add_check("r-version", "environment", r_version, "4.4.1", r_version == "4.4.1")
for (i in seq_len(nrow(package_audit))) {
  add_check(
    paste0("package-", package_audit$package[[i]]),
    "environment",
    package_audit$observed_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$status[[i]] == "PASS"
  )
}

command_names <- c("pdfinfo", "pdffonts", "identify")
command_audit <- data.frame(
  command = command_names,
  path = vapply(command_names, command_path, character(1)),
  stringsAsFactors = FALSE
)
command_audit$status <- ifelse(nzchar(command_audit$path), "PASS", "FAIL")
for (i in seq_len(nrow(command_audit))) {
  add_check(
    paste0("command-", command_audit$command[[i]]),
    "environment",
    sanitize_text(command_audit$path[[i]]),
    "available",
    command_audit$status[[i]] == "PASS"
  )
}

font_family <- font_pub
font_info <- systemfonts::font_info(font_family)
font_path <- font_info$path[[1L]]
font_glyphs <- systemfonts::glyph_info(
  paste(
    "Environment Mantel variation partition Hellinger edaphic climate",
    "geographic wetland group correlation adjusted residual 0123456789"
  ),
  family = font_family
)
font_audit <- data.frame(
  requested_family = font_family,
  resolved_family = font_info$family[[1L]],
  style = font_info$style[[1L]],
  path = sanitize_text(font_path),
  file_sha256 = if (file.exists(font_path)) sha256_file(font_path) else "",
  scalable = font_info$scalable[[1L]],
  glyphs_checked = nrow(font_glyphs),
  missing_glyphs = sum(font_glyphs$index == 0L),
  stringsAsFactors = FALSE
)
add_check(
  "font-resolved-family", "font", font_info$family[[1L]], font_family,
  identical(font_info$family[[1L]], font_family)
)
add_check("font-file-exists", "font", file.exists(font_path), TRUE, file.exists(font_path))
add_check("font-scalable", "font", font_info$scalable[[1L]], TRUE, isTRUE(font_info$scalable[[1L]]))
add_check(
  "font-glyph-coverage", "font", font_audit$missing_glyphs, 0,
  font_audit$missing_glyphs == 0L
)

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
environment_frame <- read_keyed_tsv(input_paths[["environment"]])
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])

counts_character <- as.matrix(otu_frame)
counts <- matrix(
  suppressWarnings(as.numeric(counts_character)),
  nrow = nrow(counts_character),
  ncol = ncol(counts_character),
  dimnames = dimnames(counts_character)
)
environment_character <- as.matrix(environment_frame)
environment_all <- matrix(
  suppressWarnings(as.numeric(environment_character)),
  nrow = nrow(environment_character),
  ncol = ncol(environment_character),
  dimnames = dimnames(environment_character)
)

feature_ids <- rownames(counts)
sample_ids <- colnames(counts)
rank_columns <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
environment_variables <- c(
  "Latitude", "Longitude", "Altitude", "Temperature", "Precipitation",
  "TOC", "NH4", "NO3", "pH", "Conductivity", "TN"
)

add_check("feature-count", "data", nrow(counts), 13628, nrow(counts) == 13628L)
add_check("sample-count", "data", ncol(counts), 90, ncol(counts) == 90L)
add_check("read-count", "data", sum(counts), 1619670, sum(counts) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628L)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90L)
add_check(
  "environment-source-sample-count", "data", nrow(environment_all), 200,
  nrow(environment_all) == 200L
)
add_check(
  "environment-variable-count", "data", ncol(environment_all), 11,
  ncol(environment_all) == 11L
)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "taxonomy-feature-ids", "data",
  length(intersect(feature_ids, rownames(taxonomy))), length(feature_ids),
  setequal(feature_ids, rownames(taxonomy))
)
add_check(
  "metadata-sample-ids", "data",
  length(intersect(sample_ids, rownames(metadata))), length(sample_ids),
  setequal(sample_ids, rownames(metadata))
)
add_check(
  "environment-sample-ids", "data",
  length(intersect(sample_ids, rownames(environment_all))), length(sample_ids),
  all(sample_ids %in% rownames(environment_all))
)
add_check("counts-finite", "data", sum(is.finite(counts)), length(counts), all(is.finite(counts)))
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0))
add_check(
  "counts-integer", "data", max(abs(counts - round(counts))), 0,
  all(abs(counts - round(counts)) < .Machine$double.eps^0.5)
)
add_check("nonempty-features", "data", sum(rowSums(counts) > 0), nrow(counts), all(rowSums(counts) > 0))
add_check("nonempty-samples", "data", sum(colSums(counts) > 0), ncol(counts), all(colSums(counts) > 0))
add_check(
  "taxonomy-rank-columns", "data", intersect(rank_columns, colnames(taxonomy)),
  rank_columns, all(rank_columns %in% colnames(taxonomy))
)
add_check(
  "metadata-required-columns", "data",
  intersect(c("Group", "Type", "Saline"), colnames(metadata)),
  c("Group", "Type", "Saline"),
  all(c("Group", "Type", "Saline") %in% colnames(metadata))
)
add_check(
  "environment-required-columns", "data",
  intersect(environment_variables, colnames(environment_all)),
  environment_variables,
  all(environment_variables %in% colnames(environment_all))
)
add_check(
  "source-summary-feature-count", "provenance", source_summary$features, 13628,
  identical(as.integer(source_summary$features), 13628L)
)
add_check(
  "source-summary-common-samples", "provenance",
  source_summary$common_environment_samples, 90,
  identical(as.integer(source_summary$common_environment_samples), 90L)
)

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
environment <- as.data.frame(
  environment_all[sample_ids, environment_variables, drop = FALSE],
  check.names = FALSE
)
group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
metadata$GroupCode <- metadata$Group
metadata$Group <- factor(
  unname(group_labels[metadata$GroupCode]),
  levels = unname(group_labels)
)
metadata$Saline <- factor(metadata$Saline)
metadata$Type <- factor(metadata$Type)

sample_alignment_audit <- data.frame(
  SampleID = sample_ids,
  CountTable = TRUE,
  Metadata = sample_ids %in% rownames(metadata),
  Environment = sample_ids %in% rownames(environment_all),
  LibraryReads = as.numeric(colSums(counts)),
  Group = as.character(metadata$Group),
  Type = as.character(metadata$Type),
  Saline = as.character(metadata$Saline),
  CompleteEnvironment = stats::complete.cases(environment),
  stringsAsFactors = FALSE
)
add_check(
  "aligned-analysis-samples", "alignment", nrow(sample_alignment_audit), 90,
  nrow(sample_alignment_audit) == 90L
)
add_check(
  "all-environment-complete", "alignment",
  sum(sample_alignment_audit$CompleteEnvironment), 90,
  all(sample_alignment_audit$CompleteEnvironment)
)
add_check(
  "group-balance", "design", as.integer(table(metadata$Group)), rep(30L, 3L),
  identical(as.integer(table(metadata$Group)), rep(30L, 3L))
)
group_saline_cells <- nrow(unique(data.frame(
  Group = metadata$Group,
  Saline = metadata$Saline
)))
add_check(
  "group-saline-observed-cells", "design", group_saline_cells, 3,
  group_saline_cells == 3L
)
group_determines_saline <- all(vapply(
  split(as.character(metadata$Saline), metadata$Group),
  function(x) length(unique(x)) == 1L,
  logical(1)
))
add_check(
  "group-determines-salinity", "design", group_determines_saline, TRUE,
  group_determines_saline
)
add_check(
  "no-repeat-design-columns", "design",
  length(intersect(
    c("SubjectID", "PairID", "Block", "SiteID", "Time", "Batch"),
    colnames(metadata)
  )),
  0,
  length(intersect(
    c("SubjectID", "PairID", "Block", "SiteID", "Time", "Batch"),
    colnames(metadata)
  )) == 0L
)

environment_variable_audit <- data.frame(
  Variable = environment_variables,
  Role = c(
    "Geographic", "Geographic", "Geographic context", "Climate", "Climate",
    "Edaphic", "Edaphic", "Edaphic", "Edaphic", "Edaphic", "Edaphic sensitivity"
  ),
  UsedInEnvfit = TRUE,
  UsedInPrimaryMantel = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  UsedInPrimaryVPA = environment_variables %in% c(
    "pH", "TOC", "NH4", "NO3", "Conductivity"
  ),
  Missing = colSums(is.na(environment)),
  UniqueValues = vapply(environment, function(x) length(unique(x)), integer(1)),
  Mean = colMeans(environment),
  SD = vapply(environment, stats::sd, numeric(1)),
  Minimum = vapply(environment, min, numeric(1)),
  Maximum = vapply(environment, max, numeric(1)),
  stringsAsFactors = FALSE
)
add_check(
  "environment-all-numeric", "environment",
  sum(vapply(as.data.frame(environment), is.numeric, logical(1))), 11,
  all(vapply(as.data.frame(environment), is.numeric, logical(1)))
)
add_check(
  "environment-no-zero-variance", "environment",
  sum(environment_variable_audit$SD > 0), 11,
  all(environment_variable_audit$SD > 0)
)

environment_z <- as.data.frame(scale(environment))
add_check(
  "environment-zscore-means", "environment",
  max(abs(colMeans(environment_z))), "<1e-12",
  max(abs(colMeans(environment_z))) < 1e-12
)
add_check(
  "environment-zscore-sd", "environment",
  max(abs(vapply(environment_z, stats::sd, numeric(1)) - 1)), "<1e-12",
  max(abs(vapply(environment_z, stats::sd, numeric(1)) - 1)) < 1e-12
)

correlation_matrix <- stats::cor(environment, method = "pearson")
correlation_index <- expand.grid(
  Variable1 = rownames(correlation_matrix),
  Variable2 = colnames(correlation_matrix),
  stringsAsFactors = FALSE
)
correlation_index$PearsonR <- as.vector(correlation_matrix)
correlation_index$AbsR <- abs(correlation_index$PearsonR)
correlation_index$Diagonal <- correlation_index$Variable1 == correlation_index$Variable2
variable_order <- setNames(seq_along(environment_variables), environment_variables)
correlation_index$UniquePair <- variable_order[correlation_index$Variable1] <
  variable_order[correlation_index$Variable2]
correlation_index$HighCollinearity <- correlation_index$UniquePair &
  correlation_index$AbsR >= 0.80
high_correlation_pairs <- correlation_index[
  correlation_index$HighCollinearity,
  c("Variable1", "Variable2", "PearsonR", "AbsR"),
  drop = FALSE
]
add_check(
  "high-correlation-pair-count", "collinearity",
  nrow(high_correlation_pairs), 4, nrow(high_correlation_pairs) == 4L
)
toc_tn_r <- correlation_matrix["TOC", "TN"]
add_check(
  "toc-tn-correlation", "collinearity", round(toc_tn_r, 9), 0.988019536,
  near(toc_tn_r, 0.988019536, tolerance = 1e-8)
)

primary_seed <- 20260725L
permutation_count <- 999L
permutation_control <- permute::how(
  within = permute::Within(type = "free"),
  nperm = permutation_count
)
set.seed(primary_seed)
permutation_matrix <- permute::shuffleSet(
  n = length(sample_ids),
  nset = permutation_count,
  control = permutation_control
)
set.seed(primary_seed)
permutation_matrix_repeat <- permute::shuffleSet(
  n = length(sample_ids),
  nset = permutation_count,
  control = permutation_control
)
valid_permutation_rows <- apply(
  permutation_matrix,
  1L,
  function(x) identical(sort(as.integer(x)), seq_along(sample_ids))
)
identity_rows <- apply(
  permutation_matrix,
  1L,
  function(x) identical(as.integer(x), seq_along(sample_ids))
)
unique_permutation_rows <- nrow(unique(as.data.frame(permutation_matrix)))
permutation_hash <- digest::digest(permutation_matrix, algo = "sha256")
permutation_audit <- data.frame(
  Seed = primary_seed,
  Samples = length(sample_ids),
  Mode = "Free sample-label permutations",
  Reason = paste(
    "Metadata has no subject, pair, block, site, time or batch identifier;",
    "samples are treated as independent for this public-data example"
  ),
  Requested = permutation_count,
  Generated = nrow(permutation_matrix),
  UniqueRows = unique_permutation_rows,
  IdentityRows = sum(identity_rows),
  LegalRows = sum(valid_permutation_rows),
  MinimumP = 1 / (permutation_count + 1),
  MatrixSHA256 = permutation_hash,
  DeterministicRepeat = identical(permutation_matrix, permutation_matrix_repeat),
  stringsAsFactors = FALSE
)
add_check(
  "permutation-row-count", "permutation", nrow(permutation_matrix), 999,
  nrow(permutation_matrix) == 999L
)
add_check(
  "permutation-unique-rows", "permutation", unique_permutation_rows, 999,
  unique_permutation_rows == 999L
)
add_check(
  "permutation-identity-rows", "permutation", sum(identity_rows), 0,
  sum(identity_rows) == 0L
)
add_check(
  "permutation-legal-rows", "permutation", sum(valid_permutation_rows), 999,
  all(valid_permutation_rows)
)
add_check(
  "permutation-deterministic", "permutation",
  identical(permutation_matrix, permutation_matrix_repeat), TRUE,
  identical(permutation_matrix, permutation_matrix_repeat)
)

min_prevalence <- 0.05
minimum_samples <- ceiling(min_prevalence * ncol(counts))
primary_keep <- rowSums(counts > 0) >= minimum_samples
primary_counts <- counts[primary_keep, , drop = FALSE]
community_counts <- t(primary_counts)
community_relative <- sweep(community_counts, 1L, rowSums(community_counts), "/")
community_hellinger <- vegan::decostand(community_counts, method = "hellinger")
bray_distance <- vegan::vegdist(community_relative, method = "bray")
add_check(
  "primary-feature-count", "standardization", sum(primary_keep), 11113,
  sum(primary_keep) == 11113L
)
add_check(
  "primary-minimum-samples", "standardization", minimum_samples, 5,
  minimum_samples == 5L
)
add_check(
  "relative-row-sums", "standardization",
  max(abs(rowSums(community_relative) - 1)), "<1e-12",
  max(abs(rowSums(community_relative) - 1)) < 1e-12
)
add_check(
  "hellinger-finite", "standardization",
  sum(is.finite(community_hellinger)), length(community_hellinger),
  all(is.finite(community_hellinger))
)
add_check(
  "bray-distance-count", "distance", length(bray_distance), choose(90, 2),
  length(bray_distance) == choose(90, 2)
)
add_check(
  "bray-distance-range", "distance", round(range(bray_distance), 6), "[0,1]",
  min(bray_distance) >= 0 && max(bray_distance) <= 1
)

pcoa <- vegan::wcmdscale(
  bray_distance,
  k = 2,
  eig = TRUE,
  add = "lingoes",
  x.ret = TRUE
)
pcoa_scores <- data.frame(
  SampleID = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1L],
  PCoA2 = pcoa$points[, 2L],
  Group = metadata[rownames(pcoa$points), "Group"],
  stringsAsFactors = FALSE
)
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
pcoa_axis_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)
add_check("pcoa-score-rows", "ordination", nrow(pcoa_scores), 90, nrow(pcoa_scores) == 90L)
add_check(
  "pcoa-finite", "ordination", sum(is.finite(pcoa$points)), length(pcoa$points),
  all(is.finite(pcoa$points))
)
add_check(
  "pcoa-lingoes-nonnegative", "ordination", pcoa$ac, ">=0",
  is.finite(pcoa$ac) && pcoa$ac >= 0
)

envfit_result <- vegan::envfit(
  pcoa$points,
  environment_z,
  permutations = permutation_matrix,
  w = rep(1, nrow(pcoa$points)),
  na.rm = FALSE
)
envfit_scaled_scores <- vegan::scores(envfit_result, display = "vectors")
envfit_results <- data.frame(
  Variable = rownames(envfit_result$vectors$arrows),
  Axis1Direction = envfit_result$vectors$arrows[, 1L],
  Axis2Direction = envfit_result$vectors$arrows[, 2L],
  PlotAxis1 = envfit_scaled_scores[, 1L],
  PlotAxis2 = envfit_scaled_scores[, 2L],
  TwoAxisR2 = as.numeric(envfit_result$vectors$r),
  PValue = as.numeric(envfit_result$vectors$pvals),
  stringsAsFactors = FALSE
)
envfit_results$PAdjustedBH <- stats::p.adjust(envfit_results$PValue, method = "BH")
envfit_results$RejectBH05 <- envfit_results$PAdjustedBH < 0.05
envfit_results$TestFamilySize <- nrow(envfit_results)
envfit_results$Interpretation <- paste(
  "Association with the displayed two-axis PCoA configuration;",
  "not variance explained in the full community table"
)
envfit_results <- envfit_results[
  order(envfit_results$PAdjustedBH, -envfit_results$TwoAxisR2),
  ,
  drop = FALSE
]
rownames(envfit_results) <- NULL
add_check(
  "envfit-variable-count", "envfit", nrow(envfit_results), 11,
  nrow(envfit_results) == 11L
)
add_check(
  "envfit-bh-family", "multiplicity",
  unique(envfit_results$TestFamilySize), 11,
  identical(unique(envfit_results$TestFamilySize), 11L)
)
add_check(
  "envfit-significant-count", "envfit",
  sum(envfit_results$RejectBH05), 10,
  sum(envfit_results$RejectBH05) == 10L
)
temperature_r2 <- envfit_results$TwoAxisR2[
  envfit_results$Variable == "Temperature"
]
add_check(
  "envfit-temperature-r2", "envfit", round(temperature_r2, 9), 0.80081147,
  near(temperature_r2, 0.80081147, tolerance = 1e-8)
)
nh4_rejected <- envfit_results$RejectBH05[envfit_results$Variable == "NH4"]
add_check(
  "envfit-nh4-not-rejected", "envfit", nh4_rejected, FALSE,
  isFALSE(unname(nh4_rejected))
)

edaphic_names <- c("pH", "TOC", "NH4", "NO3", "Conductivity")
climate_names <- c("Temperature", "Precipitation")
edaphic_z <- as.data.frame(scale(environment[, edaphic_names, drop = FALSE]))
climate_z <- as.data.frame(scale(environment[, climate_names, drop = FALSE]))

haversine_dist <- function(latitude, longitude) {
  rad <- pi / 180
  lat <- latitude * rad
  lon <- longitude * rad
  n <- length(lat)
  out <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    dlat <- lat[j] - lat[i]
    dlon <- lon[j] - lon[i]
    a <- sin(dlat / 2)^2 +
      cos(lat[i]) * cos(lat[j]) * sin(dlon / 2)^2
    out[i, j] <- 6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
    out[j, i] <- out[i, j]
  }
  dimnames(out) <- list(sample_ids, sample_ids)
  stats::as.dist(out)
}

environment_distances <- list(
  Edaphic = stats::dist(edaphic_z, method = "euclidean"),
  Climate = stats::dist(climate_z, method = "euclidean"),
  Geographic = haversine_dist(environment$Latitude, environment$Longitude)
)
mantel_results <- do.call(
  rbind,
  lapply(seq_along(environment_distances), function(i) {
    fit <- vegan::mantel(
      bray_distance,
      environment_distances[[i]],
      method = "spearman",
      permutations = permutation_matrix,
      na.rm = FALSE
    )
    data.frame(
      DistanceBlock = names(environment_distances)[[i]],
      CommunityDistance = "Bray-Curtis on relative abundance",
      EnvironmentDistance = c(
        "Euclidean on five Z-scored soil variables",
        "Euclidean on two Z-scored climate variables",
        "Haversine great-circle distance (km)"
      )[[i]],
      Method = "Spearman",
      Statistic = as.numeric(fit$statistic),
      PValue = as.numeric(fit$signif),
      Permutations = permutation_count,
      Alternative = "Positive matrix association",
      stringsAsFactors = FALSE
    )
  })
)
mantel_results$PAdjustedHolm <- stats::p.adjust(mantel_results$PValue, method = "holm")
mantel_results$RejectHolm05 <- mantel_results$PAdjustedHolm < 0.05
add_check(
  "mantel-family-size", "multiplicity", nrow(mantel_results), 3,
  nrow(mantel_results) == 3L
)
add_check(
  "mantel-all-holm-rejected", "mantel",
  sum(mantel_results$RejectHolm05), 3,
  all(mantel_results$RejectHolm05)
)
expected_mantel <- c(
  Edaphic = 0.3318957,
  Climate = 0.3950914,
  Geographic = 0.5869835
)
for (block in names(expected_mantel)) {
  observed <- mantel_results$Statistic[mantel_results$DistanceBlock == block]
  add_check(
    paste0("mantel-", tolower(block), "-statistic"),
    "mantel", round(observed, 7), expected_mantel[[block]],
    near(observed, expected_mantel[[block]], tolerance = 1e-7)
  )
}
add_check(
  "mantel-geographic-largest", "mantel",
  mantel_results$DistanceBlock[which.max(mantel_results$Statistic)],
  "Geographic",
  identical(
    mantel_results$DistanceBlock[which.max(mantel_results$Statistic)],
    "Geographic"
  )
)

partial_mantel <- vegan::mantel.partial(
  bray_distance,
  environment_distances$Edaphic,
  environment_distances$Geographic,
  method = "pearson",
  permutations = permutation_matrix,
  na.rm = FALSE
)
partial_mantel_sensitivity <- data.frame(
  Analysis = "Edaphic distance | geographic distance",
  Method = "Pearson partial Mantel",
  Statistic = as.numeric(partial_mantel$statistic),
  PValue = as.numeric(partial_mantel$signif),
  Permutations = permutation_count,
  PrimaryInference = FALSE,
  Warning = paste(
    "Sensitivity analysis only; partial Mantel can have invalid error rates",
    "under spatial autocorrelation and is not used for causal adjustment"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "partial-mantel-pearson-method", "mantel",
  partial_mantel_sensitivity$Method, "Pearson partial Mantel",
  identical(partial_mantel_sensitivity$Method, "Pearson partial Mantel")
)
add_check(
  "partial-mantel-statistic", "mantel",
  round(partial_mantel_sensitivity$Statistic, 7), 0.1884762,
  near(partial_mantel_sensitivity$Statistic, 0.1884762, tolerance = 1e-7)
)
add_check(
  "partial-mantel-not-primary", "interpretation",
  partial_mantel_sensitivity$PrimaryInference, FALSE,
  identical(partial_mantel_sensitivity$PrimaryInference, FALSE)
)

group_df <- data.frame(Group = metadata$Group, row.names = sample_ids)
variation_partition_object <- vegan::varpart(
  community_hellinger,
  edaphic_z,
  group_df
)
individual_fractions <- variation_partition_object$part$indfract
variation_partition <- data.frame(
  FractionCode = c("[a]", "[b]", "[c]", "[d]"),
  Fraction = c(
    "Pure edaphic | wetland group",
    "Pure wetland group | edaphic",
    "Shared edaphic + wetland group",
    "Unexplained"
  ),
  AdjustedR2 = as.numeric(individual_fractions$Adj.R.squared),
  Testable = as.logical(individual_fractions$Testable),
  Interpretation = c(
    "Unique association of five soil variables after conditioning group",
    "Unique association of three-level wetland group after conditioning soil",
    "Overlap is arithmetic and cannot be assigned uniquely or tested directly",
    "Residual Hellinger community variation"
  ),
  stringsAsFactors = FALSE
)
expected_fractions <- c(0.09673210, 0.11404959, 0.07096358, 0.71825474)
add_check(
  "vpa-fraction-count", "variation-partition", nrow(variation_partition), 4,
  nrow(variation_partition) == 4L
)
add_check(
  "vpa-fraction-values", "variation-partition",
  round(variation_partition$AdjustedR2, 8), expected_fractions,
  near(variation_partition$AdjustedR2, expected_fractions, tolerance = 1e-7)
)
add_check(
  "vpa-fractions-sum-one", "variation-partition",
  sum(variation_partition$AdjustedR2), 1,
  near(sum(variation_partition$AdjustedR2), 1, tolerance = 1e-10)
)
add_check(
  "vpa-testable-pattern", "variation-partition",
  variation_partition$Testable, c(TRUE, TRUE, FALSE, FALSE),
  identical(variation_partition$Testable, c(TRUE, TRUE, FALSE, FALSE))
)

model_data <- cbind(edaphic_z, Group = metadata$Group, Saline = metadata$Saline)
combined_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity + Group,
  data = model_data
)
pure_edaphic_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity +
    Condition(Group),
  data = model_data
)
pure_group_rda <- vegan::rda(
  community_hellinger ~ Group +
    Condition(pH + TOC + NH4 + NO3 + Conductivity),
  data = model_data
)
edaphic_unconditioned_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity,
  data = model_data
)
edaphic_condition_saline_rda <- vegan::rda(
  community_hellinger ~ pH + TOC + NH4 + NO3 + Conductivity +
    Condition(Saline),
  data = model_data
)

combined_r2 <- vegan::RsquareAdj(combined_rda)
pure_edaphic_r2 <- vegan::RsquareAdj(pure_edaphic_rda)
pure_group_r2 <- vegan::RsquareAdj(pure_group_rda)
set.seed(primary_seed)
pure_edaphic_test <- anova(pure_edaphic_rda, permutations = permutation_matrix)
set.seed(primary_seed)
pure_group_test <- anova(pure_group_rda, permutations = permutation_matrix)

variation_fraction_tests <- data.frame(
  Fraction = c(
    "Pure edaphic | wetland group",
    "Pure wetland group | edaphic"
  ),
  Df = c(pure_edaphic_test$Df[[1L]], pure_group_test$Df[[1L]]),
  Variance = c(
    pure_edaphic_test$Variance[[1L]],
    pure_group_test$Variance[[1L]]
  ),
  FValue = c(pure_edaphic_test$F[[1L]], pure_group_test$F[[1L]]),
  AdjustedR2 = c(
    pure_edaphic_r2$adj.r.squared,
    pure_group_r2$adj.r.squared
  ),
  PValue = c(
    pure_edaphic_test$`Pr(>F)`[[1L]],
    pure_group_test$`Pr(>F)`[[1L]]
  ),
  Permutations = permutation_count,
  stringsAsFactors = FALSE
)
variation_fraction_tests$PAdjustedHolm <- stats::p.adjust(
  variation_fraction_tests$PValue,
  method = "holm"
)
variation_fraction_tests$RejectHolm05 <-
  variation_fraction_tests$PAdjustedHolm < 0.05
add_check(
  "combined-adjusted-r2", "variation-partition",
  round(combined_r2$adj.r.squared, 8), 0.2817453,
  near(combined_r2$adj.r.squared, 0.2817453, tolerance = 1e-7)
)
add_check(
  "pure-fraction-test-count", "multiplicity",
  nrow(variation_fraction_tests), 2,
  nrow(variation_fraction_tests) == 2L
)
add_check(
  "pure-fraction-rejections-holm", "variation-partition",
  sum(variation_fraction_tests$RejectHolm05), 2,
  all(variation_fraction_tests$RejectHolm05)
)
add_check(
  "pure-fraction-p-floor", "permutation",
  variation_fraction_tests$PValue, rep(0.001, 2),
  identical(variation_fraction_tests$PValue, rep(0.001, 2))
)

combined_vif <- vegan::vif.cca(combined_rda)
vif_audit <- data.frame(
  Term = names(combined_vif),
  VIF = as.numeric(combined_vif),
  Threshold = 5,
  Status = ifelse(combined_vif < 5, "PASS", "REVIEW"),
  stringsAsFactors = FALSE
)
add_check(
  "combined-vif-max", "collinearity",
  round(max(combined_vif), 6), "<5",
  max(combined_vif) < 5
)

conditioning_sensitivity <- data.frame(
  Conditioning = c("None", "Binary salinity", "Three-level wetland group"),
  RawR2 = c(
    vegan::RsquareAdj(edaphic_unconditioned_rda)$r.squared,
    vegan::RsquareAdj(edaphic_condition_saline_rda)$r.squared,
    pure_edaphic_r2$r.squared
  ),
  AdjustedR2 = c(
    vegan::RsquareAdj(edaphic_unconditioned_rda)$adj.r.squared,
    vegan::RsquareAdj(edaphic_condition_saline_rda)$adj.r.squared,
    pure_edaphic_r2$adj.r.squared
  ),
  Primary = c(FALSE, FALSE, TRUE),
  Boundary = c(
    "Total edaphic association",
    "Salinity is completely nested in Group and is not an independent cause",
    "Primary pure fraction; Group remains a broad composite context"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "conditioning-branch-count", "sensitivity",
  nrow(conditioning_sensitivity), 3,
  nrow(conditioning_sensitivity) == 3L
)
add_check(
  "conditioning-primary-count", "sensitivity",
  sum(conditioning_sensitivity$Primary), 1,
  sum(conditioning_sensitivity$Primary) == 1L
)

prevalence_thresholds <- c(0, 0.05, 0.10, 0.20)
prevalence_sensitivity <- do.call(
  rbind,
  lapply(prevalence_thresholds, function(threshold) {
    keep <- if (threshold == 0) {
      rep(TRUE, nrow(counts))
    } else {
      rowSums(counts > 0) >= ceiling(threshold * ncol(counts))
    }
    branch_counts <- t(counts[keep, , drop = FALSE])
    branch_relative <- sweep(branch_counts, 1L, rowSums(branch_counts), "/")
    branch_hellinger <- vegan::decostand(branch_counts, method = "hellinger")
    branch_vpa <- vegan::varpart(branch_hellinger, edaphic_z, group_df)
    fractions <- as.numeric(branch_vpa$part$indfract$Adj.R.squared)
    branch_mantel <- vegan::mantel(
      vegan::vegdist(branch_relative, method = "bray"),
      environment_distances$Edaphic,
      method = "spearman",
      permutations = permutation_matrix
    )
    data.frame(
      PrevalenceThreshold = threshold,
      MinimumSamples = ifelse(threshold == 0, 1L, ceiling(threshold * 90)),
      Features = sum(keep),
      PureEdaphicAdjustedR2 = fractions[[1L]],
      PureGroupAdjustedR2 = fractions[[2L]],
      SharedAdjustedR2 = fractions[[3L]],
      UnexplainedAdjustedR2 = fractions[[4L]],
      EdaphicMantelRho = as.numeric(branch_mantel$statistic),
      EdaphicMantelP = as.numeric(branch_mantel$signif),
      Primary = threshold == min_prevalence,
      stringsAsFactors = FALSE
    )
  })
)
add_check(
  "prevalence-branch-count", "sensitivity",
  nrow(prevalence_sensitivity), 4,
  nrow(prevalence_sensitivity) == 4L
)
add_check(
  "prevalence-feature-counts", "sensitivity",
  prevalence_sensitivity$Features, c(13628, 11113, 9043, 5700),
  identical(prevalence_sensitivity$Features, c(13628L, 11113L, 9043L, 5700L))
)
add_check(
  "prevalence-pure-edaphic-range", "sensitivity",
  range(prevalence_sensitivity$PureEdaphicAdjustedR2), "0.09--0.11",
  min(prevalence_sensitivity$PureEdaphicAdjustedR2) > 0.09 &&
    max(prevalence_sensitivity$PureEdaphicAdjustedR2) < 0.11
)
add_check(
  "prevalence-mantel-positive", "sensitivity",
  prevalence_sensitivity$EdaphicMantelRho, ">0",
  all(prevalence_sensitivity$EdaphicMantelRho > 0)
)
add_check(
  "prevalence-mantel-p-floor", "sensitivity",
  prevalence_sensitivity$EdaphicMantelP, rep(0.001, 4),
  identical(prevalence_sensitivity$EdaphicMantelP, rep(0.001, 4))
)

correlation_plot_data <- correlation_index[
  variable_order[correlation_index$Variable1] <=
    variable_order[correlation_index$Variable2],
  ,
  drop = FALSE
]
correlation_plot_data$Variable1 <- factor(
  correlation_plot_data$Variable1,
  levels = rev(environment_variables)
)
correlation_plot_data$Variable2 <- factor(
  correlation_plot_data$Variable2,
  levels = environment_variables
)
correlation_plot_data$Label <- ifelse(
  correlation_plot_data$Diagonal,
  "1.00",
  ifelse(
    correlation_plot_data$AbsR >= 0.50,
    sprintf("%.2f", correlation_plot_data$PearsonR),
    ""
  )
)
correlation_plot <- ggplot2::ggplot(
  correlation_plot_data,
  ggplot2::aes(Variable2, Variable1, fill = PearsonR)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_tile(
    data = correlation_plot_data[
      correlation_plot_data$HighCollinearity,
      ,
      drop = FALSE
    ],
    fill = NA,
    colour = "#111111",
    linewidth = 0.85
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    size = 2.5,
    colour = "#111111",
    family = font_family
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Environmental covariates are not independent",
    subtitle = "Black outlines mark |r| >= 0.80; labels are shown for |r| >= 0.50",
    x = NULL,
    y = NULL,
    caption = paste(
      "All 11 variables enter the descriptive envfit family.",
      "TN is excluded from the primary edaphic set because it is nearly redundant with TOC.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 8.4) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, size = 7.4),
    axis.text.y = ggplot2::element_text(size = 7.4),
    legend.position = "right"
  )

group_palette <- c(
  "Inland wetland" = "#0072B2",
  "Coastal wetland" = "#D55E00",
  "Tibetan Plateau" = "#009E73"
)
group_shapes <- c(
  "Inland wetland" = 16,
  "Coastal wetland" = 17,
  "Tibetan Plateau" = 15
)
envfit_plot_arrows <- envfit_results[
  envfit_results$RejectBH05,
  ,
  drop = FALSE
]
envfit_plot_arrows <- head(
  envfit_plot_arrows[order(-envfit_plot_arrows$TwoAxisR2), , drop = FALSE],
  7L
)
arrow_multiplier <- min(
  diff(range(pcoa_scores$PCoA1)) * 0.40 /
    max(abs(envfit_plot_arrows$PlotAxis1)),
  diff(range(pcoa_scores$PCoA2)) * 0.40 /
    max(abs(envfit_plot_arrows$PlotAxis2))
)
envfit_plot_arrows$ArrowX <- envfit_plot_arrows$PlotAxis1 * arrow_multiplier
envfit_plot_arrows$ArrowY <- envfit_plot_arrows$PlotAxis2 * arrow_multiplier
envfit_label_offsets <- data.frame(
  Variable = c(
    "Temperature", "Latitude", "Altitude", "Longitude",
    "Conductivity", "Precipitation", "pH"
  ),
  OffsetX = c(0.025, -0.015, -0.025, 0.005, -0.060, -0.035, 0.020),
  OffsetY = c(0.030, -0.045, 0.035, -0.030, 0.035, -0.055, 0.040),
  Hjust = c(0, 1, 0.5, 0.5, 1, 0.5, 0),
  stringsAsFactors = FALSE
)
offset_index <- match(envfit_plot_arrows$Variable, envfit_label_offsets$Variable)
envfit_plot_arrows$LabelX <-
  envfit_plot_arrows$ArrowX + envfit_label_offsets$OffsetX[offset_index]
envfit_plot_arrows$LabelY <-
  envfit_plot_arrows$ArrowY + envfit_label_offsets$OffsetY[offset_index]
envfit_plot_arrows$LabelHjust <- envfit_label_offsets$Hjust[offset_index]
envfit_plot <- ggplot2::ggplot(
  pcoa_scores,
  ggplot2::aes(PCoA1, PCoA2, colour = Group, shape = Group)
) +
  ggplot2::stat_ellipse(
    type = "t",
    level = 0.95,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 2.4, alpha = 0.83) +
  ggplot2::geom_segment(
    data = envfit_plot_arrows,
    ggplot2::aes(x = 0, y = 0, xend = ArrowX, yend = ArrowY),
    inherit.aes = FALSE,
    colour = "#222222",
    linewidth = 0.55,
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  ) +
  ggplot2::geom_segment(
    data = envfit_plot_arrows,
    ggplot2::aes(x = ArrowX, y = ArrowY, xend = LabelX, yend = LabelY),
    inherit.aes = FALSE,
    colour = "#777777",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = envfit_plot_arrows,
    ggplot2::aes(
      x = LabelX, y = LabelY, label = Variable, hjust = LabelHjust
    ),
    inherit.aes = FALSE,
    colour = "#111111",
    family = font_family,
    fontface = "bold",
    size = 2.8
  ) +
  ggplot2::scale_colour_manual(values = group_palette, name = "Wetland group") +
  ggplot2::scale_shape_manual(values = group_shapes, name = "Wetland group") +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(
    title = "envfit overlays external gradients on an unconstrained PCoA",
    subtitle = paste(
      "Top 7 of",
      sum(envfit_results$RejectBH05),
      "BH-significant variables; all 11 tests remain in the audit table"
    ),
    x = sprintf("PCoA1 (%.1f%% positive inertia)", pcoa_axis_percent[[1L]]),
    y = sprintf("PCoA2 (%.1f%% positive inertia)", pcoa_axis_percent[[2L]]),
    caption = paste(
      "Arrow direction and length describe association with this two-axis display.",
      "They are not causal effects or full-community explained fractions.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

mantel_plot_data <- data.frame(
  Label = c(
    "Geographic distance",
    "Climate distance",
    "Edaphic distance",
    "Edaphic | geographic"
  ),
  Statistic = c(
    mantel_results$Statistic[mantel_results$DistanceBlock == "Geographic"],
    mantel_results$Statistic[mantel_results$DistanceBlock == "Climate"],
    mantel_results$Statistic[mantel_results$DistanceBlock == "Edaphic"],
    partial_mantel_sensitivity$Statistic
  ),
  PValue = c(
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Geographic"],
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Climate"],
    mantel_results$PAdjustedHolm[mantel_results$DistanceBlock == "Edaphic"],
    partial_mantel_sensitivity$PValue
  ),
  Family = c(rep("Primary Mantel family", 3L), "Partial Mantel sensitivity"),
  stringsAsFactors = FALSE
)
mantel_plot_data$Label <- factor(
  mantel_plot_data$Label,
  levels = rev(mantel_plot_data$Label)
)
mantel_plot_data$Annotation <- sprintf(
  "r = %.3f; P %s",
  mantel_plot_data$Statistic,
  format_p(mantel_plot_data$PValue)
)
mantel_plot <- ggplot2::ggplot(
  mantel_plot_data,
  ggplot2::aes(Statistic, Label, colour = Family, shape = Family)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Statistic, yend = Label),
    linewidth = 0.8,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::geom_text(
    ggplot2::aes(label = Annotation),
    hjust = -0.08,
    colour = "#111111",
    family = font_family,
    size = 3.1,
    show.legend = FALSE
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.4) +
  ggplot2::scale_colour_manual(
    values = c(
      "Primary Mantel family" = "#0072B2",
      "Partial Mantel sensitivity" = "#D55E00"
    ),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Primary Mantel family" = 16,
      "Partial Mantel sensitivity" = 17
    ),
    name = NULL
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 0.72), clip = "off") +
  ggplot2::labs(
    title = "Mantel statistics compare two distance structures",
    subtitle = "Primary Spearman tests use Holm correction across three prespecified blocks",
    x = "Matrix correlation",
    y = NULL,
    caption = paste(
      "The partial Mantel row uses Pearson correlation and is a sensitivity analysis only.",
      "Distance correlation does not identify an independent environmental cause.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "top",
    plot.margin = ggplot2::margin(5.5, 38, 5.5, 5.5)
  )

variation_plot_data <- variation_partition
variation_plot_data$Fraction <- factor(
  variation_plot_data$Fraction,
  levels = rev(variation_plot_data$Fraction)
)
variation_palette <- c(
  "Pure edaphic | wetland group" = "#0072B2",
  "Pure wetland group | edaphic" = "#D55E00",
  "Shared edaphic + wetland group" = "#CC79A7",
  "Unexplained" = "#9CA3AF"
)
variation_plot <- ggplot2::ggplot(
  variation_plot_data,
  ggplot2::aes(AdjustedR2, Fraction, fill = Fraction)
) +
  ggplot2::geom_col(width = 0.64, colour = "#333333", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(AdjustedR2, accuracy = 0.1),
      hjust = ifelse(AdjustedR2 >= 0, -0.08, 1.08)
    ),
    family = font_family,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = variation_palette, guide = "none") +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.02, 0.13))
  ) +
  ggplot2::labs(
    title = "Adjusted variation partition separates unique and shared fractions",
    subtitle = sprintf(
      "Combined adjusted R2 = %.3f; both pure fractions Holm P <= 0.002",
      combined_r2$adj.r.squared
    ),
    x = "Adjusted fraction of Hellinger community variation",
    y = NULL,
    caption = paste(
      "Shared variation is not directly testable or attributable to one predictor set.",
      "Wetland group is a broad context variable and is completely confounded with salinity here.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 18, 5.5, 5.5)
  )

figure_specification <- data.frame(
  stem = c(
    "25-environment-correlation",
    "25-envfit-pcoa",
    "25-mantel-distance-blocks",
    "25-variation-partition"
  ),
  width_mm = rep(183, 4L),
  height_mm = c(126, 122, 96, 100),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  correlation_plot,
  envfit_plot,
  mantel_plot,
  variation_plot
)
for (i in seq_len(nrow(figure_specification))) {
  save_pub(
    plot = plot_objects[[i]],
    file_base = file.path(figure_dir, figure_specification$stem[[i]]),
    width = figure_specification$width_mm[[i]],
    height = figure_specification$height_mm[[i]],
    units = "mm",
    dpi = figure_specification$raster_ppi[[i]],
    write_svg = TRUE,
    write_tiff = TRUE,
    base_family = font_family
  )
}

parse_pdf_info <- function(path) {
  result <- run_command("pdfinfo", shQuote(path))
  if (result$status != 0L) stop("pdfinfo failed for ", path, call. = FALSE)
  page_line <- grep("^Pages:", result$output, value = TRUE)
  size_line <- grep("^Page size:", result$output, value = TRUE)
  size_match <- regexec(
    "Page size:[[:space:]]+([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts",
    size_line
  )
  size_parts <- regmatches(size_line, size_match)[[1L]]
  if (length(page_line) != 1L || length(size_parts) != 3L) {
    stop("Unexpected pdfinfo output for ", path, call. = FALSE)
  }
  data.frame(
    pages = as.integer(sub("^Pages:[[:space:]]+", "", page_line)),
    width_pt = as.numeric(size_parts[[2L]]),
    height_pt = as.numeric(size_parts[[3L]]),
    stringsAsFactors = FALSE
  )
}

parse_svg_info <- function(path) {
  svg_text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  root_match <- regexec(
    "<svg[^>]*width='([0-9.]+)pt'[^>]*height='([0-9.]+)pt'",
    svg_text
  )
  root_parts <- regmatches(svg_text, root_match)[[1L]]
  if (length(root_parts) != 3L) {
    stop("Unexpected SVG root dimensions for ", path, call. = FALSE)
  }
  data.frame(
    width_pt = as.numeric(root_parts[[2L]]),
    height_pt = as.numeric(root_parts[[3L]]),
    text_nodes = lengths(regmatches(svg_text, gregexpr("<text\\b", svg_text, perl = TRUE))),
    font_family_present = grepl('font-family: "DejaVu Sans"', svg_text, fixed = TRUE),
    chinese_characters = lengths(regmatches(svg_text, gregexpr("[一-鿿]", svg_text, perl = TRUE))),
    stringsAsFactors = FALSE
  )
}

identify_raster <- function(path) {
  result <- run_command(
    "identify",
    c("-format", shQuote("%m|%w|%h|%x|%y|%U|%C"), shQuote(path))
  )
  if (result$status != 0L || length(result$output) == 0L) {
    stop("identify failed for ", path, call. = FALSE)
  }
  fields <- strsplit(result$output[[1L]], "\\|", fixed = FALSE)[[1L]]
  if (length(fields) != 7L) stop("Unexpected identify output for ", path, call. = FALSE)
  data.frame(
    detected_format = fields[[1L]],
    pixel_width = as.integer(fields[[2L]]),
    pixel_height = as.integer(fields[[3L]]),
    density_x = as.numeric(fields[[4L]]),
    density_y = as.numeric(fields[[5L]]),
    density_unit = fields[[6L]],
    compression = fields[[7L]],
    stringsAsFactors = FALSE
  )
}

file_signature_ok <- function(path, extension) {
  raw <- readBin(path, what = "raw", n = 16L)
  if (extension == "pdf") {
    return(identical(rawToChar(raw[seq_len(4L)]), "%PDF"))
  }
  if (extension == "png") {
    return(identical(
      as.integer(raw[seq_len(8L)]),
      c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)
    ))
  }
  if (extension == "tiff") {
    signature <- as.integer(raw[seq_len(4L)])
    return(
      identical(signature, c(73L, 73L, 42L, 0L)) ||
        identical(signature, c(77L, 77L, 0L, 42L))
    )
  }
  if (extension == "svg") {
    prefix <- rawToChar(raw)
    return(startsWith(prefix, "<?xml") || startsWith(prefix, "<svg"))
  }
  FALSE
}

format_rows <- list()
for (spec_i in seq_len(nrow(figure_specification))) {
  spec <- figure_specification[spec_i, , drop = FALSE]
  for (extension in c("pdf", "svg", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(spec$stem, ".", extension))
    exists <- file.exists(path)
    bytes <- if (exists) file.info(path)$size else 0
    add_check(paste0(spec$stem, "-", extension, "-exists"), "export", exists, TRUE, exists)
    add_check(
      paste0(spec$stem, "-", extension, "-nonempty"),
      "export", bytes, ">1000", bytes > 1000
    )
    signature_ok <- exists && file_signature_ok(path, extension)
    add_check(
      paste0(spec$stem, "-", extension, "-signature"),
      "export", signature_ok, TRUE, signature_ok
    )

    actual_width_mm <- NA_real_
    actual_height_mm <- NA_real_
    pixel_width <- NA_integer_
    pixel_height <- NA_integer_
    effective_ppi_x <- NA_real_
    effective_ppi_y <- NA_real_
    detected_format <- toupper(extension)
    compression <- NA_character_
    pages <- NA_integer_
    text_nodes <- NA_integer_
    chinese_characters <- NA_integer_
    font_status <- NA_character_

    if (extension == "pdf") {
      pdf_info <- parse_pdf_info(path)
      actual_width_mm <- pdf_info$width_pt * 25.4 / 72
      actual_height_mm <- pdf_info$height_pt * 25.4 / 72
      pages <- pdf_info$pages
      font_result <- run_command("pdffonts", shQuote(path))
      font_output <- paste(font_result$output, collapse = "\n")
      font_embedded <- font_result$status == 0L &&
        grepl("DejaVuSans", font_output, fixed = TRUE) &&
        grepl("yes[[:space:]]+yes[[:space:]]+yes", font_output)
      font_status <- ifelse(font_embedded, "embedded", "missing")
      add_check(
        paste0(spec$stem, "-pdf-single-page"), "vector", pages, 1,
        pages == 1L
      )
      add_check(
        paste0(spec$stem, "-pdf-font-embedded"), "vector",
        font_status, "embedded", font_embedded
      )
    } else if (extension == "svg") {
      svg_info <- parse_svg_info(path)
      actual_width_mm <- svg_info$width_pt * 25.4 / 72
      actual_height_mm <- svg_info$height_pt * 25.4 / 72
      text_nodes <- svg_info$text_nodes
      chinese_characters <- svg_info$chinese_characters
      font_status <- ifelse(svg_info$font_family_present, "declared", "missing")
      add_check(
        paste0(spec$stem, "-svg-text-nodes"), "vector",
        text_nodes, ">0", text_nodes > 0L
      )
      add_check(
        paste0(spec$stem, "-svg-font-declared"), "vector",
        font_status, "declared", svg_info$font_family_present
      )
      add_check(
        paste0(spec$stem, "-svg-english-only"), "vector",
        chinese_characters, 0, chinese_characters == 0L
      )
    } else {
      raster_info <- identify_raster(path)
      detected_format <- raster_info$detected_format
      pixel_width <- raster_info$pixel_width
      pixel_height <- raster_info$pixel_height
      compression <- raster_info$compression
      actual_width_mm <- spec$width_mm
      actual_height_mm <- spec$height_mm
      effective_ppi_x <- pixel_width / (spec$width_mm / 25.4)
      effective_ppi_y <- pixel_height / (spec$height_mm / 25.4)
      expected_width_px <- floor(spec$width_mm / 25.4 * spec$raster_ppi)
      expected_height_px <- floor(spec$height_mm / 25.4 * spec$raster_ppi)
      add_check(
        paste0(spec$stem, "-", extension, "-pixel-dimensions"),
        "raster", c(pixel_width, pixel_height),
        c(expected_width_px, expected_height_px),
        all(c(pixel_width, pixel_height) == c(expected_width_px, expected_height_px))
      )
      add_check(
        paste0(spec$stem, "-", extension, "-effective-ppi"),
        "raster", round(c(effective_ppi_x, effective_ppi_y), 3),
        paste0(spec$raster_ppi, " +/- 0.5"),
        all(abs(c(effective_ppi_x, effective_ppi_y) - spec$raster_ppi) <= 0.5)
      )
      add_check(
        paste0(spec$stem, "-", extension, "-format"),
        "raster", detected_format, toupper(extension),
        identical(detected_format, toupper(extension))
      )
      if (extension == "tiff") {
        add_check(
          paste0(spec$stem, "-tiff-lzw"), "raster",
          compression, "LZW", identical(compression, "LZW")
        )
      }
    }

    dimension_pass <-
      abs(actual_width_mm - spec$width_mm) <= 0.4 &&
      abs(actual_height_mm - spec$height_mm) <= 0.4
    add_check(
      paste0(spec$stem, "-", extension, "-physical-dimensions"),
      "export", round(c(actual_width_mm, actual_height_mm), 3),
      c(spec$width_mm, spec$height_mm), dimension_pass
    )

    format_rows[[length(format_rows) + 1L]] <- data.frame(
      figure = spec$stem,
      extension = extension,
      path = file.path("figures", paste0(spec$stem, ".", extension)),
      detected_format = detected_format,
      bytes = bytes,
      target_width_mm = spec$width_mm,
      target_height_mm = spec$height_mm,
      actual_width_mm = actual_width_mm,
      actual_height_mm = actual_height_mm,
      pixel_width = pixel_width,
      pixel_height = pixel_height,
      effective_ppi_x = effective_ppi_x,
      effective_ppi_y = effective_ppi_y,
      compression = compression,
      pages = pages,
      text_nodes = text_nodes,
      chinese_characters = chinese_characters,
      font_status = font_status,
      dimension_status = ifelse(dimension_pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
}
format_audit <- do.call(rbind, format_rows)
add_check(
  "primary-figure-count", "export",
  length(unique(format_audit$figure)), 4,
  length(unique(format_audit$figure)) == 4L
)
add_check(
  "primary-format-count", "export", nrow(format_audit), 16,
  nrow(format_audit) == 16L
)
add_check(
  "vector-file-count", "export",
  sum(format_audit$extension %in% c("pdf", "svg")), 8,
  sum(format_audit$extension %in% c("pdf", "svg")) == 8L
)
add_check(
  "raster-file-count", "export",
  sum(format_audit$extension %in% c("png", "tiff")), 8,
  sum(format_audit$extension %in% c("png", "tiff")) == 8L
)
add_check(
  "all-physical-dimensions-pass", "export",
  sum(format_audit$dimension_status == "PASS"), 16,
  all(format_audit$dimension_status == "PASS")
)
add_check(
  "all-svg-english-only", "export",
  sum(format_audit$chinese_characters[
    format_audit$extension == "svg"
  ], na.rm = TRUE),
  0,
  all(format_audit$chinese_characters[
    format_audit$extension == "svg"
  ] == 0L)
)

validation_checks <- do.call(rbind, check_rows)
checks_total <- nrow(validation_checks)
checks_passed <- sum(validation_checks$status == "PASS")
checks_failed <- checks_total - checks_passed

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(command_audit, "command-audit.tsv")
write_tsv(font_audit, "font-audit.tsv")
write_tsv(sample_alignment_audit, "sample-alignment-audit.tsv")
write_tsv(environment_variable_audit, "environmental-variable-audit.tsv")
write_tsv(correlation_index, "environmental-correlation.tsv")
write_tsv(high_correlation_pairs, "high-correlation-pairs.tsv")
write_tsv(permutation_audit, "permutation-audit.tsv")
write_tsv(envfit_results, "envfit-results.tsv")
write_tsv(mantel_results, "mantel-results.tsv")
write_tsv(partial_mantel_sensitivity, "partial-mantel-sensitivity.tsv")
write_tsv(variation_partition, "variation-partition.tsv")
write_tsv(variation_fraction_tests, "variation-fraction-tests.tsv")
write_tsv(vif_audit, "vif-audit.tsv")
write_tsv(conditioning_sensitivity, "conditioning-sensitivity.tsv")
write_tsv(prevalence_sensitivity, "prevalence-sensitivity.tsv")
write_tsv(pcoa_scores, "pcoa-scores.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 25L,
  title = "Environmental fitting, Mantel tests, and variation partitioning",
  dataset = list(
    source = "microeco 2.0.0 packaged wetland 16S tables",
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    environment_source_samples = nrow(environment_all),
    matched_samples = nrow(environment),
    environment_variables = ncol(environment),
    groups = nlevels(metadata$Group)
  ),
  standardization = list(
    primary_prevalence = min_prevalence,
    minimum_samples = minimum_samples,
    retained_features = sum(primary_keep),
    envfit_mantel = "Relative-abundance Bray-Curtis",
    variation_partition = "Hellinger-transformed community matrix"
  ),
  permutations = list(
    seed = primary_seed,
    mode = "Free sample-label permutations",
    count = permutation_count,
    minimum_p = 1 / (permutation_count + 1),
    unique_rows = unique_permutation_rows,
    identity_rows = sum(identity_rows),
    matrix_sha256 = permutation_hash
  ),
  environment = list(
    envfit_family = nrow(envfit_results),
    envfit_rejections_bh = sum(envfit_results$RejectBH05),
    high_correlation_pairs_abs_r_ge_080 = nrow(high_correlation_pairs),
    edaphic_variables = edaphic_names,
    climate_variables = climate_names,
    group_saline_completely_confounded = group_determines_saline
  ),
  mantel = list(
    family_size = nrow(mantel_results),
    rejections_holm = sum(mantel_results$RejectHolm05),
    edaphic_rho = mantel_results$Statistic[
      mantel_results$DistanceBlock == "Edaphic"
    ],
    climate_rho = mantel_results$Statistic[
      mantel_results$DistanceBlock == "Climate"
    ],
    geographic_rho = mantel_results$Statistic[
      mantel_results$DistanceBlock == "Geographic"
    ],
    partial_edaphic_given_geographic_pearson =
      partial_mantel_sensitivity$Statistic,
    partial_mantel_primary = FALSE
  ),
  variation_partition = list(
    combined_raw_r2 = combined_r2$r.squared,
    combined_adjusted_r2 = combined_r2$adj.r.squared,
    pure_edaphic_adjusted_r2 = variation_partition$AdjustedR2[[1L]],
    pure_group_adjusted_r2 = variation_partition$AdjustedR2[[2L]],
    shared_adjusted_r2 = variation_partition$AdjustedR2[[3L]],
    unexplained_adjusted_r2 = variation_partition$AdjustedR2[[4L]],
    testable_fractions = sum(variation_partition$Testable),
    pure_fraction_rejections_holm = sum(
      variation_fraction_tests$RejectHolm05
    ),
    maximum_vif = max(combined_vif)
  ),
  sensitivity = list(
    prevalence_branches = nrow(prevalence_sensitivity),
    conditioning_branches = nrow(conditioning_sensitivity),
    pure_edaphic_adjusted_r2_range = range(
      prevalence_sensitivity$PureEdaphicAdjustedR2
    ),
    edaphic_mantel_rho_range = range(
      prevalence_sensitivity$EdaphicMantelRho
    )
  ),
  graphics = list(
    font_family = font_family,
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    format_files = nrow(format_audit),
    raster_ppi = unique(figure_specification$raster_ppi),
    svg_files_with_text = sum(
      format_audit$extension == "svg" & format_audit$text_nodes > 0
    ),
    svg_chinese_characters = sum(
      format_audit$chinese_characters[
        format_audit$extension == "svg"
      ],
      na.rm = TRUE
    ),
    pdf_files_with_embedded_font = sum(
      format_audit$extension == "pdf" &
        format_audit$font_status == "embedded"
    ),
    lzw_tiff_files = sum(
      format_audit$extension == "tiff" &
        format_audit$compression == "LZW"
    ),
    dimension_checks_passed = sum(
      format_audit$dimension_status == "PASS"
    )
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "environment-variance-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

validation_lines <- c(
  "Article 25 environment-variance validation",
  paste0("Status: ", ifelse(checks_failed == 0L, "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total),
  paste0("Samples/features: 90/", sum(primary_keep)),
  paste0("Permutation matrix: ", permutation_count, " unique free rows"),
  paste0(
    "envfit BH rejections: ",
    sum(envfit_results$RejectBH05),
    "/",
    nrow(envfit_results)
  ),
  paste0(
    "Mantel rho (edaphic/climate/geographic): ",
    paste(sprintf("%.4f", mantel_results$Statistic), collapse = "/")
  ),
  paste0(
    "VPA adjusted fractions (edaphic/group/shared/residual): ",
    paste(sprintf("%.4f", variation_partition$AdjustedR2), collapse = "/")
  ),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("R version: ", r_version),
  paste0("vegan version: ", observed_versions[["vegan"]])
)
writeLines(validation_lines, file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed_ids <- validation_checks$check_id[
    validation_checks$status == "FAIL"
  ]
  stop(
    "Article 25 validation failed: ",
    paste(failed_ids, collapse = ", "),
    call. = FALSE
  )
}

cat(
  sprintf(
    "Article 25 validation passed: %d/%d checks.\n",
    checks_passed,
    checks_total
  )
)
