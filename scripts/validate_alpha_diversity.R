#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260721)

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

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "digest", "dplyr", "ggplot2", "iNEXT", "jsonlite", "phyloseq",
  "ragg", "readr", "scales", "svglite", "systemfonts", "tibble",
  "tidyr", "vegan"
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
    na = character()
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
  if (!is.null(command_status)) {
    status <- as.integer(command_status)
  }
  list(status = status, output = output)
}

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  source_summary = file.path(input_dir, "source_summary.json")
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  source_summary = "acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 19 input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

observed_sha256 <- vapply(input_paths, sha256_file, character(1))
input_audit <- data.frame(
  asset = names(input_paths),
  file = basename(input_paths),
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
  dplyr = "1.1.4",
  ggplot2 = "3.5.2",
  iNEXT = "3.0.2",
  jsonlite = "1.8.8",
  phyloseq = "1.48.0",
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
  tibble = "3.2.1",
  tidyr = "1.3.1",
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
add_check(
  "r-version",
  "environment",
  paste(R.version$major, R.version$minor, sep = "."),
  "4.4.1",
  identical(paste(R.version$major, R.version$minor, sep = "."), "4.4.1")
)
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
    command_audit$path[[i]],
    "available",
    command_audit$status[[i]] == "PASS"
  )
}

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])

otu_matrix <- suppressWarnings(data.matrix(otu_frame))
storage.mode(otu_matrix) <- "numeric"
feature_ids <- rownames(otu_matrix)
sample_ids <- colnames(otu_matrix)

add_check("feature-count", "data", nrow(otu_matrix), 13628, nrow(otu_matrix) == 13628L)
add_check("sample-count", "data", ncol(otu_matrix), 90, ncol(otu_matrix) == 90L)
add_check("read-count", "data", sum(otu_matrix), 1619670, sum(otu_matrix) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628L)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90L)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "taxonomy-feature-ids",
  "data",
  length(intersect(feature_ids, rownames(taxonomy))),
  length(feature_ids),
  setequal(feature_ids, rownames(taxonomy))
)
add_check(
  "metadata-sample-ids",
  "data",
  length(intersect(sample_ids, rownames(metadata))),
  length(sample_ids),
  setequal(sample_ids, rownames(metadata))
)
add_check("counts-finite", "data", all(is.finite(otu_matrix)), TRUE, all(is.finite(otu_matrix)))
add_check("counts-nonnegative", "data", min(otu_matrix), ">=0", min(otu_matrix) >= 0)
add_check(
  "counts-integer",
  "data",
  max(abs(otu_matrix - round(otu_matrix))),
  0,
  all(abs(otu_matrix - round(otu_matrix)) < .Machine$double.eps^0.5)
)
add_check(
  "source-summary-shape",
  "data",
  c(source_summary$features, source_summary$samples),
  c(13628, 90),
  identical(
    as.integer(c(source_summary$features, source_summary$samples)),
    c(13628L, 90L)
  )
)
add_check(
  "metadata-columns",
  "data",
  intersect(c("Group", "Type", "Saline"), colnames(metadata)),
  c("Group", "Type", "Saline"),
  all(c("Group", "Type", "Saline") %in% colnames(metadata))
)

metadata <- metadata[sample_ids, , drop = FALSE]
community <- t(otu_matrix)
group_counts <- table(metadata$Group)
add_check(
  "group-levels",
  "data",
  sort(names(group_counts)),
  c("CW", "IW", "TW"),
  identical(sort(names(group_counts)), c("CW", "IW", "TW"))
)
add_check(
  "group-balance",
  "data",
  as.integer(group_counts[c("CW", "IW", "TW")]),
  c(30, 30, 30),
  identical(as.integer(group_counts[c("CW", "IW", "TW")]), c(30L, 30L, 30L))
)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
group_palette <- c(
  "Inland wetland" = pal_pub[["blue"]],
  "Coastal wetland" = pal_pub[["vermillion"]],
  "Tibetan Plateau" = pal_pub[["green"]]
)
group_shapes <- c(
  "Inland wetland" = 21,
  "Coastal wetland" = 22,
  "Tibetan Plateau" = 24
)

font_family <- font_pub
font_info <- systemfonts::font_info(font_family)
font_path <- font_info$path[[1L]]
font_glyphs <- systemfonts::glyph_info(
  "Alpha diversity Hill numbers coverage 0123456789",
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
  "font-resolved-family",
  "font",
  font_info$family[[1L]],
  font_family,
  identical(font_info$family[[1L]], font_family)
)
add_check("font-file-exists", "font", file.exists(font_path), TRUE, file.exists(font_path))
add_check("font-scalable", "font", font_info$scalable[[1L]], TRUE, isTRUE(font_info$scalable[[1L]]))
add_check(
  "font-glyph-coverage",
  "font",
  font_audit$missing_glyphs,
  0,
  font_audit$missing_glyphs == 0L
)

library_size <- rowSums(community)
observed <- vegan::specnumber(community)
frequency_one <- rowSums(community == 1)
frequency_two <- rowSums(community == 2)
goods_coverage <- 1 - frequency_one / library_size
estimate_r <- vegan::estimateR(community)
shannon <- vegan::diversity(community, index = "shannon")
simpson <- vegan::diversity(community, index = "simpson")
inverse_simpson <- vegan::diversity(community, index = "invsimpson")
pielou <- shannon / log(observed)

ps <- phyloseq::phyloseq(phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE))
phyloseq_alpha <- suppressWarnings(
  phyloseq::estimate_richness(
    ps,
    measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
  )
)
phyloseq_alpha <- phyloseq_alpha[sample_ids, , drop = FALSE]

alpha_raw <- data.frame(
  SampleID = sample_ids,
  Group = metadata$Group,
  GroupLabel = unname(group_labels[metadata$Group]),
  LibrarySize = as.numeric(library_size),
  Singletons = as.integer(frequency_one),
  Doubletons = as.integer(frequency_two),
  GoodsCoverage = as.numeric(goods_coverage),
  Observed = as.numeric(observed),
  Chao1 = as.numeric(estimate_r["S.chao1", ]),
  Chao1SE = as.numeric(estimate_r["se.chao1", ]),
  Shannon = as.numeric(shannon),
  Simpson = as.numeric(simpson),
  Pielou = as.numeric(pielou),
  Hill0 = as.numeric(observed),
  Hill1 = as.numeric(exp(shannon)),
  Hill2 = as.numeric(inverse_simpson),
  stringsAsFactors = FALSE
)
alpha_raw$GroupLabel <- factor(alpha_raw$GroupLabel, levels = names(group_palette))

depth_observed_spearman <- suppressWarnings(
  cor(alpha_raw$LibrarySize, alpha_raw$Observed, method = "spearman")
)
add_check("library-minimum", "alpha", min(library_size), 10364, min(library_size) == 10364)
add_check("library-maximum", "alpha", max(library_size), 37374, max(library_size) == 37374)
add_check(
  "depth-observed-spearman",
  "alpha",
  round(depth_observed_spearman, 6),
  ">=0.75",
  depth_observed_spearman >= 0.75
)
add_check("raw-row-count", "alpha", nrow(alpha_raw), 90, nrow(alpha_raw) == 90L)
add_check("raw-no-missing", "alpha", sum(!is.finite(as.matrix(alpha_raw[c(
  "LibrarySize", "GoodsCoverage", "Observed", "Chao1", "Shannon",
  "Simpson", "Pielou", "Hill0", "Hill1", "Hill2"
)]))), 0, all(is.finite(as.matrix(alpha_raw[c(
  "LibrarySize", "GoodsCoverage", "Observed", "Chao1", "Shannon",
  "Simpson", "Pielou", "Hill0", "Hill1", "Hill2"
)]))))
add_check(
  "singleton-preserved",
  "alpha",
  sum(alpha_raw$Singletons),
  ">0",
  sum(alpha_raw$Singletons) > 0
)
add_check(
  "observed-parity-phyloseq",
  "parity",
  max(abs(alpha_raw$Observed - phyloseq_alpha$Observed)),
  0,
  max(abs(alpha_raw$Observed - phyloseq_alpha$Observed)) < 1e-10
)
add_check(
  "chao1-parity-phyloseq",
  "parity",
  max(abs(alpha_raw$Chao1 - phyloseq_alpha$Chao1)),
  "<1e-8",
  max(abs(alpha_raw$Chao1 - phyloseq_alpha$Chao1)) < 1e-8
)
add_check(
  "shannon-parity-phyloseq",
  "parity",
  max(abs(alpha_raw$Shannon - phyloseq_alpha$Shannon)),
  "<1e-10",
  max(abs(alpha_raw$Shannon - phyloseq_alpha$Shannon)) < 1e-10
)
add_check(
  "simpson-parity-phyloseq",
  "parity",
  max(abs(alpha_raw$Simpson - phyloseq_alpha$Simpson)),
  "<1e-10",
  max(abs(alpha_raw$Simpson - phyloseq_alpha$Simpson)) < 1e-10
)
add_check(
  "inverse-simpson-parity-phyloseq",
  "parity",
  max(abs(alpha_raw$Hill2 - phyloseq_alpha$InvSimpson)),
  "<1e-10",
  max(abs(alpha_raw$Hill2 - phyloseq_alpha$InvSimpson)) < 1e-10
)
add_check(
  "hill-one-formula",
  "parity",
  max(abs(alpha_raw$Hill1 - exp(alpha_raw$Shannon))),
  0,
  max(abs(alpha_raw$Hill1 - exp(alpha_raw$Shannon))) < 1e-10
)
add_check(
  "hill-two-formula",
  "parity",
  max(abs(alpha_raw$Hill2 - 1 / (1 - alpha_raw$Simpson))),
  "<1e-8",
  max(abs(alpha_raw$Hill2 - 1 / (1 - alpha_raw$Simpson))) < 1e-8
)

abundance_list <- setNames(
  lapply(seq_len(nrow(community)), function(i) {
    as.numeric(community[i, community[i, ] > 0])
  }),
  sample_ids
)
data_info <- iNEXT::DataInfo(abundance_list, datatype = "abundance")
data_info <- data_info[match(sample_ids, data_info$Assemblage), , drop = FALSE]
alpha_raw$iNEXTCoverage <- data_info$SC
add_check(
  "coverage-estimator-parity",
  "coverage",
  max(abs(alpha_raw$GoodsCoverage - alpha_raw$iNEXTCoverage)),
  "<=0.0001",
  max(abs(alpha_raw$GoodsCoverage - alpha_raw$iNEXTCoverage)) <= 0.0001
)

size_target <- 10000L
set.seed(20260721)
size_estimate <- iNEXT::estimateD(
  abundance_list,
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "size",
  level = size_target,
  nboot = 0
)
alpha_size <- data.frame(
  SampleID = size_estimate$Assemblage,
  Group = metadata[size_estimate$Assemblage, "Group"],
  GroupLabel = unname(group_labels[metadata[size_estimate$Assemblage, "Group"]]),
  TargetReads = as.numeric(size_estimate$m),
  Method = size_estimate$Method,
  Order = as.integer(size_estimate$Order.q),
  SampleCoverage = size_estimate$SC,
  HillD = size_estimate$qD,
  stringsAsFactors = FALSE
)
alpha_size$GroupLabel <- factor(alpha_size$GroupLabel, levels = names(group_palette))
add_check("size-row-count", "standardization", nrow(alpha_size), 270, nrow(alpha_size) == 270L)
add_check(
  "size-target",
  "standardization",
  sort(unique(alpha_size$TargetReads)),
  size_target,
  identical(sort(unique(alpha_size$TargetReads)), as.numeric(size_target))
)
add_check(
  "size-samples-retained",
  "standardization",
  length(unique(alpha_size$SampleID)),
  90,
  length(unique(alpha_size$SampleID)) == 90L
)
add_check(
  "size-method",
  "standardization",
  unique(alpha_size$Method),
  "Rarefaction",
  identical(unique(alpha_size$Method), "Rarefaction")
)
add_check(
  "size-orders",
  "standardization",
  sort(unique(alpha_size$Order)),
  c(0, 1, 2),
  identical(sort(unique(alpha_size$Order)), 0:2)
)

coverage_target <- 0.90
set.seed(20260721)
coverage_estimate <- iNEXT::estimateD(
  abundance_list,
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = coverage_target,
  nboot = 0
)
alpha_coverage <- data.frame(
  SampleID = coverage_estimate$Assemblage,
  Group = metadata[coverage_estimate$Assemblage, "Group"],
  GroupLabel = unname(group_labels[metadata[coverage_estimate$Assemblage, "Group"]]),
  ImpliedReads = coverage_estimate$m,
  ObservedReads = as.numeric(library_size[coverage_estimate$Assemblage]),
  EffortRatio = coverage_estimate$m / as.numeric(library_size[coverage_estimate$Assemblage]),
  Method = coverage_estimate$Method,
  Order = as.integer(coverage_estimate$Order.q),
  TargetCoverage = coverage_estimate$SC,
  HillD = coverage_estimate$qD,
  stringsAsFactors = FALSE
)
alpha_coverage$GroupLabel <- factor(alpha_coverage$GroupLabel, levels = names(group_palette))
coverage_q0 <- alpha_coverage[alpha_coverage$Order == 0L, , drop = FALSE]
coverage_rarefaction_n <- sum(coverage_q0$Method == "Rarefaction")
coverage_extrapolation_n <- sum(coverage_q0$Method == "Extrapolation")
coverage_observed_n <- sum(coverage_q0$Method == "Observed")
coverage_max_ratio <- max(coverage_q0$EffortRatio)
add_check("coverage-row-count", "standardization", nrow(alpha_coverage), 270, nrow(alpha_coverage) == 270L)
add_check(
  "coverage-target",
  "standardization",
  unique(alpha_coverage$TargetCoverage),
  coverage_target,
  max(abs(alpha_coverage$TargetCoverage - coverage_target)) < 1e-8
)
add_check(
  "coverage-rarefaction-samples",
  "standardization",
  coverage_rarefaction_n,
  56,
  coverage_rarefaction_n == 56L
)
add_check(
  "coverage-extrapolation-samples",
  "standardization",
  coverage_extrapolation_n,
  34,
  coverage_extrapolation_n == 34L
)
add_check(
  "coverage-observed-samples",
  "standardization",
  coverage_observed_n,
  0,
  coverage_observed_n == 0L
)
add_check(
  "coverage-max-effort-ratio",
  "standardization",
  round(coverage_max_ratio, 6),
  "<=1.5",
  coverage_max_ratio <= 1.5
)
add_check(
  "coverage-orders",
  "standardization",
  sort(unique(alpha_coverage$Order)),
  c(0, 1, 2),
  identical(sort(unique(alpha_coverage$Order)), 0:2)
)

coverage_levels <- c(0.88, 0.90, 0.92)
coverage_sensitivity_rows <- lapply(coverage_levels, function(level) {
  set.seed(20260721)
  estimate <- iNEXT::estimateD(
    abundance_list,
    q = 0,
    datatype = "abundance",
    base = "coverage",
    level = level,
    nboot = 0
  )
  observed_reads <- as.numeric(library_size[estimate$Assemblage])
  data.frame(
    TargetCoverage = level,
    Samples = nrow(estimate),
    RarefactionSamples = sum(estimate$Method == "Rarefaction"),
    ObservedSamples = sum(estimate$Method == "Observed"),
    ExtrapolationSamples = sum(estimate$Method == "Extrapolation"),
    MedianImpliedReads = median(estimate$m),
    MaximumEffortRatio = max(estimate$m / observed_reads),
    stringsAsFactors = FALSE
  )
})
coverage_sensitivity <- do.call(rbind, coverage_sensitivity_rows)
add_check(
  "coverage-sensitivity-levels",
  "sensitivity",
  coverage_sensitivity$TargetCoverage,
  coverage_levels,
  identical(coverage_sensitivity$TargetCoverage, coverage_levels)
)
add_check(
  "coverage-sensitivity-monotone-effort",
  "sensitivity",
  coverage_sensitivity$MedianImpliedReads,
  "strictly increasing",
  all(diff(coverage_sensitivity$MedianImpliedReads) > 0)
)
add_check(
  "coverage-090-bounded-extrapolation",
  "sensitivity",
  coverage_sensitivity$MaximumEffortRatio[coverage_sensitivity$TargetCoverage == 0.90],
  "<=1.5",
  coverage_sensitivity$MaximumEffortRatio[coverage_sensitivity$TargetCoverage == 0.90] <= 1.5
)

depth_grid <- seq(1000L, size_target, by = 1000L)
rarefaction_matrix <- vapply(
  depth_grid,
  function(depth) as.numeric(vegan::rarefy(community, sample = depth)),
  numeric(nrow(community))
)
rownames(rarefaction_matrix) <- sample_ids
colnames(rarefaction_matrix) <- as.character(depth_grid)
rarefaction_curve_data <- data.frame(
  SampleID = rep(sample_ids, times = length(depth_grid)),
  TargetReads = rep(depth_grid, each = length(sample_ids)),
  ExpectedObserved = as.numeric(rarefaction_matrix),
  stringsAsFactors = FALSE
)
rarefaction_curve_data$Group <- metadata[rarefaction_curve_data$SampleID, "Group"]
rarefaction_curve_data$GroupLabel <- unname(group_labels[rarefaction_curve_data$Group])
rarefaction_curve_data$GroupLabel <- factor(
  rarefaction_curve_data$GroupLabel,
  levels = names(group_palette)
)
monotone_by_sample <- vapply(
  split(rarefaction_curve_data, rarefaction_curve_data$SampleID),
  function(x) all(diff(x$ExpectedObserved[order(x$TargetReads)]) >= -1e-8),
  logical(1)
)
add_check(
  "rarefaction-curve-rows",
  "rarefaction",
  nrow(rarefaction_curve_data),
  900,
  nrow(rarefaction_curve_data) == 900L
)
add_check(
  "rarefaction-curves-monotone",
  "rarefaction",
  sum(monotone_by_sample),
  90,
  all(monotone_by_sample)
)
add_check(
  "rarefaction-curves-finite",
  "rarefaction",
  all(is.finite(rarefaction_curve_data$ExpectedObserved)),
  TRUE,
  all(is.finite(rarefaction_curve_data$ExpectedObserved))
)

raw_hill <- rbind(
  data.frame(SampleID = alpha_raw$SampleID, Order = 0L, HillD = alpha_raw$Hill0),
  data.frame(SampleID = alpha_raw$SampleID, Order = 1L, HillD = alpha_raw$Hill1),
  data.frame(SampleID = alpha_raw$SampleID, Order = 2L, HillD = alpha_raw$Hill2)
)
raw_hill$Group <- metadata[raw_hill$SampleID, "Group"]
raw_hill$GroupLabel <- unname(group_labels[raw_hill$Group])

safe_spearman <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman"))
}
rank_rows <- list()
for (order_q in 0:2) {
  raw_q <- raw_hill[raw_hill$Order == order_q, c("SampleID", "HillD")]
  names(raw_q)[[2L]] <- "RawHillD"
  size_q <- alpha_size[alpha_size$Order == order_q, c("SampleID", "HillD")]
  names(size_q)[[2L]] <- "SizeHillD"
  coverage_q <- alpha_coverage[alpha_coverage$Order == order_q, c("SampleID", "HillD")]
  names(coverage_q)[[2L]] <- "CoverageHillD"
  joined <- Reduce(function(x, y) merge(x, y, by = "SampleID", sort = FALSE), list(raw_q, size_q, coverage_q))
  rank_rows[[length(rank_rows) + 1L]] <- data.frame(
    Order = order_q,
    Comparison = c("Raw vs 10,000 reads", "Raw vs 90% coverage", "10,000 reads vs 90% coverage"),
    SpearmanRho = c(
      safe_spearman(joined$RawHillD, joined$SizeHillD),
      safe_spearman(joined$RawHillD, joined$CoverageHillD),
      safe_spearman(joined$SizeHillD, joined$CoverageHillD)
    ),
    stringsAsFactors = FALSE
  )
}
rank_stability <- do.call(rbind, rank_rows)
add_check("rank-stability-rows", "sensitivity", nrow(rank_stability), 9, nrow(rank_stability) == 9L)
add_check(
  "rank-stability-finite",
  "sensitivity",
  all(is.finite(rank_stability$SpearmanRho)),
  TRUE,
  all(is.finite(rank_stability$SpearmanRho))
)
add_check(
  "rank-stability-range",
  "sensitivity",
  range(rank_stability$SpearmanRho),
  "[-1,1]",
  all(abs(rank_stability$SpearmanRho) <= 1)
)

standardized_hill <- rbind(
  transform(
    raw_hill[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "Raw counts"
  ),
  transform(
    alpha_size[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "10,000 reads"
  ),
  transform(
    alpha_coverage[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "90% coverage"
  )
)
standardized_hill$Standardization <- factor(
  standardized_hill$Standardization,
  levels = c("Raw counts", "10,000 reads", "90% coverage")
)
group_descriptive_summary <- standardized_hill |>
  dplyr::group_by(Standardization, Group, GroupLabel, Order) |>
  dplyr::summarise(
    Samples = dplyr::n(),
    Median = stats::median(HillD),
    Q25 = stats::quantile(HillD, 0.25),
    Q75 = stats::quantile(HillD, 0.75),
    .groups = "drop"
  ) |>
  as.data.frame()
add_check(
  "group-description-rows",
  "summary",
  nrow(group_descriptive_summary),
  27,
  nrow(group_descriptive_summary) == 27L
)
add_check(
  "group-description-n",
  "summary",
  sort(unique(group_descriptive_summary$Samples)),
  30,
  identical(sort(unique(group_descriptive_summary$Samples)), 30L)
)

metric_definition <- data.frame(
  Metric = c(
    "Observed", "Chao1", "Shannon", "Simpson", "Pielou",
    "Hill q=0", "Hill q=1", "Hill q=2", "Good coverage"
  ),
  Definition = c(
    "Number of detected features",
    "Singleton/doubleton-corrected richness lower-bound estimator",
    "Shannon entropy",
    "One minus Simpson concentration",
    "Shannon divided by log observed richness",
    "Observed feature richness",
    "Exponential of Shannon entropy",
    "Inverse Simpson concentration",
    "One minus singleton reads divided by library size"
  ),
  Unit = c(
    "features", "features", "entropy", "probability", "ratio",
    "effective features", "effective features", "effective features", "proportion"
  ),
  RareFeatureWeight = c(
    "high", "very high", "moderate", "low", "moderate",
    "high", "moderate", "low", "singleton-sensitive"
  ),
  stringsAsFactors = FALSE
)
add_check("metric-definition-rows", "definition", nrow(metric_definition), 9, nrow(metric_definition) == 9L)

depth_plot_data <- rbind(
  data.frame(
    SampleID = alpha_raw$SampleID,
    GroupLabel = alpha_raw$GroupLabel,
    LibrarySize = alpha_raw$LibrarySize,
    Metric = "Observed features (count)",
    Value = alpha_raw$Observed
  ),
  data.frame(
    SampleID = alpha_raw$SampleID,
    GroupLabel = alpha_raw$GroupLabel,
    LibrarySize = alpha_raw$LibrarySize,
    Metric = "Estimated sample coverage (%)",
    Value = 100 * alpha_raw$iNEXTCoverage
  )
)
depth_plot_data$Metric <- factor(
  depth_plot_data$Metric,
  levels = c("Observed features (count)", "Estimated sample coverage (%)")
)

depth_plot <- ggplot2::ggplot(
  depth_plot_data,
  ggplot2::aes(x = LibrarySize, y = Value)
) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "#4D4D4D",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = GroupLabel, shape = GroupLabel),
    colour = "#1A1A1A",
    stroke = 0.3,
    size = 2.0,
    alpha = 0.86
  ) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", nrow = 1L) +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma(accuracy = 1)) +
  ggplot2::labs(
    title = "Sequencing depth and completeness are separate diagnostics",
    subtitle = paste0(
      "Raw observed richness remains depth-associated (Spearman rho = ",
      sprintf("%.3f", depth_observed_spearman), ")"
    ),
    x = "Library size (reads)",
    y = NULL,
    fill = "Wetland group",
    shape = "Wetland group",
    caption = "Coverage describes the processed feature table; it is not absolute cell or species coverage."
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

curve_summary <- rarefaction_curve_data |>
  dplyr::group_by(GroupLabel, TargetReads) |>
  dplyr::summarise(
    Median = stats::median(ExpectedObserved),
    Q25 = stats::quantile(ExpectedObserved, 0.25),
    Q75 = stats::quantile(ExpectedObserved, 0.75),
    .groups = "drop"
  ) |>
  as.data.frame()

rarefaction_plot <- ggplot2::ggplot(
  rarefaction_curve_data,
  ggplot2::aes(x = TargetReads, y = ExpectedObserved)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = SampleID, colour = GroupLabel),
    linewidth = 0.25,
    alpha = 0.18
  ) +
  ggplot2::geom_ribbon(
    data = curve_summary,
    ggplot2::aes(
      x = TargetReads,
      ymin = Q25,
      ymax = Q75,
      fill = GroupLabel
    ),
    inherit.aes = FALSE,
    alpha = 0.16,
    colour = NA
  ) +
  ggplot2::geom_line(
    data = curve_summary,
    ggplot2::aes(y = Median, colour = GroupLabel),
    linewidth = 0.9
  ) +
  ggplot2::geom_vline(xintercept = size_target, linetype = 2, colour = "#4D4D4D", linewidth = 0.45) +
  ggplot2::facet_wrap(~GroupLabel, nrow = 1L) +
  ggplot2::scale_colour_manual(values = group_palette, guide = "none") +
  ggplot2::scale_fill_manual(values = group_palette, guide = "none") +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Sample-size rarefaction exposes unsaturated richness",
    subtitle = "Thin lines are samples; thick lines and ribbons are descriptive medians and IQRs",
    x = "Standardized reads",
    y = "Expected observed features",
    caption = "The dashed line marks 10,000 reads; no pooled curve is used as a replicate."
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "none")

raw_for_plot <- raw_hill[c("SampleID", "Order", "HillD")]
names(raw_for_plot)[[3L]] <- "RawHillD"
size_for_plot <- alpha_size[c("SampleID", "GroupLabel", "Order", "HillD")]
size_for_plot$Standardization <- "10,000 reads"
coverage_for_plot <- alpha_coverage[c("SampleID", "GroupLabel", "Order", "HillD")]
coverage_for_plot$Standardization <- "90% coverage"
comparison_plot_data <- rbind(size_for_plot, coverage_for_plot)
comparison_plot_data <- merge(
  comparison_plot_data,
  raw_for_plot,
  by = c("SampleID", "Order"),
  sort = FALSE
)
comparison_plot_data$OrderLabel <- factor(
  paste0("q = ", comparison_plot_data$Order),
  levels = paste0("q = ", 0:2)
)
comparison_plot_data$Standardization <- factor(
  comparison_plot_data$Standardization,
  levels = c("10,000 reads", "90% coverage")
)

standardization_plot <- ggplot2::ggplot(
  comparison_plot_data,
  ggplot2::aes(x = RawHillD, y = HillD)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#7A7A7A", linewidth = 0.45) +
  ggplot2::geom_point(
    ggplot2::aes(fill = GroupLabel, shape = Standardization),
    colour = "#1A1A1A",
    stroke = 0.28,
    size = 1.65,
    alpha = 0.78
  ) +
  ggplot2::facet_wrap(~OrderLabel, scales = "free", nrow = 1L) +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_shape_manual(values = c("10,000 reads" = 21, "90% coverage" = 24)) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(shape = 21, size = 2.4, alpha = 1)
    ),
    shape = ggplot2::guide_legend(
      override.aes = list(fill = "white", size = 2.4, alpha = 1)
    )
  ) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Standardization changes rare-feature-sensitive diversity most",
    subtitle = "Each point compares one standardized estimate with its raw plug-in estimate",
    x = "Raw effective features",
    y = "Standardized effective features",
    fill = "Wetland group",
    shape = "Standardization"
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

coverage_profile <- alpha_coverage
coverage_profile$OrderLabel <- coverage_profile$Order
coverage_profile_summary <- coverage_profile |>
  dplyr::group_by(GroupLabel, OrderLabel) |>
  dplyr::summarise(
    Median = stats::median(HillD),
    Q25 = stats::quantile(HillD, 0.25),
    Q75 = stats::quantile(HillD, 0.75),
    .groups = "drop"
  ) |>
  as.data.frame()

hill_plot <- ggplot2::ggplot(
  coverage_profile,
  ggplot2::aes(x = OrderLabel, y = HillD)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = SampleID, colour = GroupLabel),
    linewidth = 0.25,
    alpha = 0.16
  ) +
  ggplot2::geom_ribbon(
    data = coverage_profile_summary,
    ggplot2::aes(
      x = OrderLabel,
      ymin = Q25,
      ymax = Q75,
      fill = GroupLabel,
      group = GroupLabel
    ),
    inherit.aes = FALSE,
    alpha = 0.14,
    colour = NA
  ) +
  ggplot2::geom_line(
    data = coverage_profile_summary,
    ggplot2::aes(y = Median, group = GroupLabel, colour = GroupLabel),
    linewidth = 0.95
  ) +
  ggplot2::geom_point(
    data = coverage_profile_summary,
    ggplot2::aes(y = Median, fill = GroupLabel),
    shape = 21,
    colour = "#1A1A1A",
    stroke = 0.3,
    size = 2.2
  ) +
  ggplot2::scale_colour_manual(values = group_palette, guide = "none") +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_x_continuous(
    breaks = 0:2,
    labels = c("q = 0\nrare taxa", "q = 1\nbalanced", "q = 2\ndominant taxa")
  ) +
  ggplot2::scale_y_log10(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Hill profiles put all three orders in effective-feature units",
    subtitle = "Diversity standardized to 90% sample coverage",
    x = "Hill order",
    y = "Effective features (log scale)",
    fill = "Wetland group",
    caption = "Lines and IQRs are descriptive; group hypothesis tests are not performed here."
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

figure_specification <- data.frame(
  stem = c(
    "19-depth-completeness-audit",
    "19-rarefaction-curves",
    "19-standardization-comparison",
    "19-hill-profile"
  ),
  width_mm = c(183, 183, 183, 183),
  height_mm = c(92, 100, 100, 92),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(depth_plot, rarefaction_plot, standardization_plot, hill_plot)
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
  if (result$status != 0L) {
    stop("pdfinfo failed for ", path, call. = FALSE)
  }
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
  if (length(fields) != 7L) {
    stop("Unexpected identify output for ", path, call. = FALSE)
  }
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
      "export",
      bytes,
      ">1000",
      bytes > 1000
    )
    signature_ok <- exists && file_signature_ok(path, extension)
    add_check(
      paste0(spec$stem, "-", extension, "-signature"),
      "export",
      signature_ok,
      TRUE,
      signature_ok
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
      add_check(paste0(spec$stem, "-pdf-single-page"), "vector", pages, 1, pages == 1L)
      add_check(
        paste0(spec$stem, "-pdf-font-embedded"),
        "vector",
        font_status,
        "embedded",
        font_embedded
      )
    } else if (extension == "svg") {
      svg_info <- parse_svg_info(path)
      actual_width_mm <- svg_info$width_pt * 25.4 / 72
      actual_height_mm <- svg_info$height_pt * 25.4 / 72
      text_nodes <- svg_info$text_nodes
      font_status <- ifelse(svg_info$font_family_present, "declared", "missing")
      add_check(
        paste0(spec$stem, "-svg-text-nodes"),
        "vector",
        text_nodes,
        ">0",
        text_nodes > 0L
      )
      add_check(
        paste0(spec$stem, "-svg-font-declared"),
        "vector",
        font_status,
        "declared",
        svg_info$font_family_present
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
        "raster",
        c(pixel_width, pixel_height),
        c(expected_width_px, expected_height_px),
        all(c(pixel_width, pixel_height) == c(expected_width_px, expected_height_px))
      )
      add_check(
        paste0(spec$stem, "-", extension, "-effective-ppi"),
        "raster",
        round(c(effective_ppi_x, effective_ppi_y), 3),
        paste0(spec$raster_ppi, " +/- 0.5"),
        all(abs(c(effective_ppi_x, effective_ppi_y) - spec$raster_ppi) <= 0.5)
      )
      add_check(
        paste0(spec$stem, "-", extension, "-format"),
        "raster",
        detected_format,
        toupper(extension),
        identical(detected_format, toupper(extension))
      )
      if (extension == "tiff") {
        add_check(
          paste0(spec$stem, "-tiff-lzw"),
          "raster",
          compression,
          "LZW",
          identical(compression, "LZW")
        )
      }
    }

    dimension_pass <-
      abs(actual_width_mm - spec$width_mm) <= 0.4 &&
      abs(actual_height_mm - spec$height_mm) <= 0.4
    add_check(
      paste0(spec$stem, "-", extension, "-physical-dimensions"),
      "export",
      round(c(actual_width_mm, actual_height_mm), 3),
      c(spec$width_mm, spec$height_mm),
      dimension_pass
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
      font_status = font_status,
      dimension_status = ifelse(dimension_pass, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
}
format_audit <- do.call(rbind, format_rows)
add_check(
  "primary-figure-count",
  "export",
  length(unique(format_audit$figure)),
  4,
  length(unique(format_audit$figure)) == 4L
)
add_check("primary-format-count", "export", nrow(format_audit), 16, nrow(format_audit) == 16L)
add_check(
  "vector-file-count",
  "export",
  sum(format_audit$extension %in% c("pdf", "svg")),
  8,
  sum(format_audit$extension %in% c("pdf", "svg")) == 8L
)
add_check(
  "raster-file-count",
  "export",
  sum(format_audit$extension %in% c("png", "tiff")),
  8,
  sum(format_audit$extension %in% c("png", "tiff")) == 8L
)
add_check(
  "all-physical-dimensions-pass",
  "export",
  sum(format_audit$dimension_status == "PASS"),
  16,
  all(format_audit$dimension_status == "PASS")
)

validation_checks <- do.call(rbind, check_rows)
checks_total <- nrow(validation_checks)
checks_passed <- sum(validation_checks$status == "PASS")
checks_failed <- checks_total - checks_passed

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(command_audit, "command-audit.tsv")
write_tsv(font_audit, "font-audit.tsv")
write_tsv(alpha_raw, "alpha-metrics-raw.tsv")
write_tsv(
  alpha_raw[c(
    "SampleID", "Group", "GroupLabel", "LibrarySize", "Singletons",
    "Doubletons", "GoodsCoverage", "iNEXTCoverage"
  )],
  "library-depth-audit.tsv"
)
write_tsv(alpha_size, "alpha-metrics-size-10000.tsv")
write_tsv(alpha_coverage, "alpha-metrics-coverage-090.tsv")
write_tsv(rarefaction_curve_data, "rarefaction-curve-data.tsv")
write_tsv(coverage_sensitivity, "coverage-sensitivity.tsv")
write_tsv(rank_stability, "standardization-rank-stability.tsv")
write_tsv(group_descriptive_summary, "group-descriptive-summary.tsv")
write_tsv(metric_definition, "metric-definition.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 19L,
  title = "Alpha diversity, rarefaction, coverage, and Hill numbers",
  dataset = list(
    source = "microeco 2.0.0 wetland 16S example",
    features = ncol(community),
    samples = nrow(community),
    reads = sum(community),
    groups = length(group_counts),
    samples_per_group = as.integer(group_counts[c("CW", "IW", "TW")]),
    min_depth = min(library_size),
    median_depth = unname(stats::median(library_size)),
    max_depth = max(library_size)
  ),
  raw = list(
    depth_observed_spearman = unname(depth_observed_spearman),
    coverage_min = min(alpha_raw$iNEXTCoverage),
    coverage_median = unname(stats::median(alpha_raw$iNEXTCoverage)),
    coverage_max = max(alpha_raw$iNEXTCoverage),
    observed_min = min(alpha_raw$Observed),
    observed_median = unname(stats::median(alpha_raw$Observed)),
    observed_max = max(alpha_raw$Observed),
    singleton_reads = sum(alpha_raw$Singletons),
    doubleton_features = sum(alpha_raw$Doubletons)
  ),
  size_standardization = list(
    target_reads = size_target,
    samples_retained = length(unique(alpha_size$SampleID)),
    samples_dropped = nrow(community) - length(unique(alpha_size$SampleID)),
    total_reads_retained = size_target * length(unique(alpha_size$SampleID)),
    total_reads_discarded = sum(community) - size_target * length(unique(alpha_size$SampleID)),
    fraction_reads_discarded = 1 - size_target * length(unique(alpha_size$SampleID)) / sum(community)
  ),
  coverage_standardization = list(
    target_coverage = coverage_target,
    rarefaction_samples = coverage_rarefaction_n,
    observed_samples = coverage_observed_n,
    extrapolation_samples = coverage_extrapolation_n,
    median_implied_reads = unname(stats::median(coverage_q0$ImpliedReads)),
    max_effort_ratio = unname(coverage_max_ratio)
  ),
  hill = list(
    orders = 3L,
    q0 = "observed richness",
    q1 = "exp(Shannon entropy)",
    q2 = "inverse Simpson concentration",
    minimum_rank_spearman = min(rank_stability$SpearmanRho)
  ),
  graphics = list(
    font_family = font_family,
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    format_files = nrow(format_audit),
    raster_ppi = 600L,
    svg_files_with_text = sum(format_audit$extension == "svg" & format_audit$text_nodes > 0),
    pdf_files_with_embedded_font = sum(format_audit$extension == "pdf" & format_audit$font_status == "embedded"),
    lzw_tiff_files = sum(format_audit$extension == "tiff" & format_audit$compression == "LZW"),
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "alpha-diversity-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 10
)

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))
log_lines <- c(
  "Article 19 alpha-diversity validation",
  paste0("R version: ", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("Dataset: ", ncol(community), " features x ", nrow(community), " samples"),
  paste0("Reads: ", sum(community)),
  paste0("Library depth: ", min(library_size), "-", max(library_size)),
  paste0("Raw depth-observed Spearman rho: ", round(depth_observed_spearman, 6)),
  paste0("Size standardization: ", size_target, " reads; ", length(unique(alpha_size$SampleID)), " samples"),
  paste0(
    "Coverage standardization: ", coverage_target,
    "; rarefaction=", coverage_rarefaction_n,
    "; extrapolation=", coverage_extrapolation_n,
    "; max effort ratio=", round(coverage_max_ratio, 6)
  ),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("Checks passed: ", checks_passed, "/", checks_total),
  "",
  paste(validation_checks$check_id, validation_checks$status, validation_checks$observed, sep = "\t")
)
writeLines(sanitize_text(log_lines), file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop("Article 19 validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
}

message("Article 19 validation passed: ", checks_passed, "/", checks_total, " checks.")
