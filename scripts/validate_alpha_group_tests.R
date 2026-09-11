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
  "digest", "dplyr", "ggplot2", "iNEXT", "jsonlite", "ragg",
  "readr", "scales", "svglite", "systemfonts", "tidyr", "vegan"
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

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

sample_skewness <- function(x) {
  n <- length(x)
  s <- stats::sd(x)
  if (n < 3L || !is.finite(s) || s == 0) {
    return(NA_real_)
  }
  n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / s)^3)
}

cliffs_delta <- function(x, y) {
  differences <- outer(x, y, "-")
  (sum(differences > 0) - sum(differences < 0)) / length(differences)
}

bootstrap_delta <- function(x, y, replicates, seed) {
  set.seed(seed)
  values <- vapply(
    seq_len(replicates),
    function(i) {
      xb <- sample(x, length(x), replace = TRUE)
      yb <- sample(y, length(y), replace = TRUE)
      cliffs_delta(xb, yb)
    },
    numeric(1)
  )
  c(
    lower = unname(stats::quantile(values, 0.025, names = FALSE)),
    upper = unname(stats::quantile(values, 0.975, names = FALSE)),
    valid = sum(is.finite(values))
  )
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
    "Missing Article 20 input(s): ",
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
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0",
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
  "design",
  sort(names(group_counts)),
  c("CW", "IW", "TW"),
  identical(sort(names(group_counts)), c("CW", "IW", "TW"))
)
add_check(
  "group-balance",
  "design",
  as.integer(group_counts[c("CW", "IW", "TW")]),
  c(30, 30, 30),
  identical(as.integer(group_counts[c("CW", "IW", "TW")]), c(30L, 30L, 30L))
)

subject_candidates <- c("SubjectID", "Subject", "ParticipantID")
pair_candidates <- c("PairID", "Pair", "MatchedSet")
time_candidates <- c("Time", "Timepoint", "Visit")
block_candidates <- c("Block", "Batch", "SiteID", "Center")
present_subject <- intersect(subject_candidates, colnames(metadata))
present_pair <- intersect(pair_candidates, colnames(metadata))
present_time <- intersect(time_candidates, colnames(metadata))
present_block <- intersect(block_candidates, colnames(metadata))
selected_branch <- "Independent three-group"

design_audit <- data.frame(
  branch = c(
    "Independent groups",
    "Paired two-condition",
    "Repeated measures",
    "Blocked or clustered"
  ),
  required_metadata = c(
    "Group with independent experimental units",
    "SubjectID or PairID plus two conditions",
    "SubjectID plus Time and Condition",
    "Block, SiteID, Center, Family, Cage, or Batch"
  ),
  columns_detected = c(
    "Group",
    paste(c(present_subject, present_pair), collapse = ", "),
    paste(c(present_subject, present_time), collapse = ", "),
    paste(present_block, collapse = ", ")
  ),
  available = c(TRUE, FALSE, FALSE, FALSE),
  selected = c(TRUE, FALSE, FALSE, FALSE),
  action = c(
    "Kruskal-Wallis plus planned pairwise rank comparisons",
    "Do not run without a pairing identifier",
    "Do not run without subject and time identifiers",
    "Do not infer exchangeability without the clustering column"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "no-subject-identifier",
  "design",
  length(present_subject),
  0,
  length(present_subject) == 0L
)
add_check(
  "no-pair-identifier",
  "design",
  length(present_pair),
  0,
  length(present_pair) == 0L
)
add_check(
  "no-time-identifier",
  "design",
  length(present_time),
  0,
  length(present_time) == 0L
)
add_check(
  "selected-design-branch",
  "design",
  selected_branch,
  "Independent three-group",
  identical(selected_branch, "Independent three-group")
)
add_check(
  "single-selected-branch",
  "design",
  sum(design_audit$selected),
  1,
  sum(design_audit$selected) == 1L
)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
group_levels <- unname(group_labels[c("IW", "CW", "TW")])
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
  paste(
    "Independent paired repeated alpha diversity Hill numbers",
    "Kruskal Wallis Holm Cliff delta bootstrap 0123456789"
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
shannon <- vegan::diversity(community, index = "shannon")
inverse_simpson <- vegan::diversity(community, index = "invsimpson")

alpha_raw <- rbind(
  data.frame(SampleID = sample_ids, Order = 0L, HillD = as.numeric(observed)),
  data.frame(SampleID = sample_ids, Order = 1L, HillD = as.numeric(exp(shannon))),
  data.frame(SampleID = sample_ids, Order = 2L, HillD = as.numeric(inverse_simpson))
)
alpha_raw$Group <- metadata[alpha_raw$SampleID, "Group"]
alpha_raw$GroupLabel <- factor(
  unname(group_labels[alpha_raw$Group]),
  levels = group_levels
)
alpha_raw$ObservedReads <- as.numeric(library_size[alpha_raw$SampleID])
alpha_raw$Standardization <- "Raw counts"

abundance_list <- setNames(
  lapply(seq_len(nrow(community)), function(i) {
    as.numeric(community[i, community[i, ] > 0])
  }),
  sample_ids
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
  ObservedReads = as.numeric(library_size[size_estimate$Assemblage]),
  ImpliedReads = size_estimate$m,
  Method = size_estimate$Method,
  Order = as.integer(size_estimate$Order.q),
  SampleCoverage = size_estimate$SC,
  HillD = size_estimate$qD,
  Standardization = "10,000 reads",
  stringsAsFactors = FALSE
)
alpha_size$GroupLabel <- factor(alpha_size$GroupLabel, levels = group_levels)

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
  ObservedReads = as.numeric(library_size[coverage_estimate$Assemblage]),
  ImpliedReads = coverage_estimate$m,
  EffortRatio = coverage_estimate$m /
    as.numeric(library_size[coverage_estimate$Assemblage]),
  Method = coverage_estimate$Method,
  Order = as.integer(coverage_estimate$Order.q),
  TargetCoverage = coverage_estimate$SC,
  HillD = coverage_estimate$qD,
  Standardization = "90% coverage",
  stringsAsFactors = FALSE
)
alpha_coverage$GroupLabel <- factor(alpha_coverage$GroupLabel, levels = group_levels)

coverage_q0 <- alpha_coverage[alpha_coverage$Order == 0L, , drop = FALSE]
coverage_rarefaction_n <- sum(coverage_q0$Method == "Rarefaction")
coverage_extrapolation_n <- sum(coverage_q0$Method == "Extrapolation")
coverage_max_ratio <- max(coverage_q0$EffortRatio)

add_check("raw-row-count", "alpha", nrow(alpha_raw), 270, nrow(alpha_raw) == 270L)
add_check("size-row-count", "alpha", nrow(alpha_size), 270, nrow(alpha_size) == 270L)
add_check("coverage-row-count", "alpha", nrow(alpha_coverage), 270, nrow(alpha_coverage) == 270L)
add_check(
  "alpha-all-samples",
  "alpha",
  length(unique(alpha_coverage$SampleID)),
  90,
  length(unique(alpha_coverage$SampleID)) == 90L
)
add_check(
  "alpha-orders",
  "alpha",
  sort(unique(alpha_coverage$Order)),
  0:2,
  identical(sort(unique(alpha_coverage$Order)), 0:2)
)
add_check(
  "coverage-target",
  "alpha",
  unique(alpha_coverage$TargetCoverage),
  coverage_target,
  max(abs(alpha_coverage$TargetCoverage - coverage_target)) < 1e-8
)
add_check(
  "coverage-rarefaction-samples",
  "alpha",
  coverage_rarefaction_n,
  56,
  coverage_rarefaction_n == 56L
)
add_check(
  "coverage-extrapolation-samples",
  "alpha",
  coverage_extrapolation_n,
  34,
  coverage_extrapolation_n == 34L
)
add_check(
  "coverage-bounded-effort",
  "alpha",
  round(coverage_max_ratio, 6),
  "<=1.5",
  coverage_max_ratio <= 1.5
)
add_check(
  "alpha-finite-positive",
  "alpha",
  all(is.finite(alpha_coverage$HillD) & alpha_coverage$HillD > 0),
  TRUE,
  all(is.finite(alpha_coverage$HillD) & alpha_coverage$HillD > 0)
)

group_descriptive_summary <- alpha_coverage |>
  dplyr::group_by(Order, Group, GroupLabel) |>
  dplyr::summarise(
    Samples = dplyr::n(),
    Mean = mean(HillD),
    SD = stats::sd(HillD),
    Median = stats::median(HillD),
    Q25 = stats::quantile(HillD, 0.25),
    Q75 = stats::quantile(HillD, 0.75),
    Minimum = min(HillD),
    Maximum = max(HillD),
    .groups = "drop"
  ) |>
  as.data.frame()
add_check(
  "descriptive-summary-rows",
  "summary",
  nrow(group_descriptive_summary),
  9,
  nrow(group_descriptive_summary) == 9L
)
add_check(
  "descriptive-summary-balanced",
  "summary",
  unique(group_descriptive_summary$Samples),
  30,
  identical(unique(group_descriptive_summary$Samples), 30L)
)

distribution_rows <- lapply(
  split(alpha_coverage, interaction(alpha_coverage$Order, alpha_coverage$GroupLabel)),
  function(x) {
    q25 <- unname(stats::quantile(x$HillD, 0.25))
    q75 <- unname(stats::quantile(x$HillD, 0.75))
    iqr <- q75 - q25
    data.frame(
      Order = unique(x$Order),
      Group = unique(x$Group),
      GroupLabel = as.character(unique(x$GroupLabel)),
      Samples = nrow(x),
      UniqueValues = length(unique(x$HillD)),
      Skewness = sample_skewness(x$HillD),
      IQR = iqr,
      LowerOutliers = sum(x$HillD < q25 - 1.5 * iqr),
      UpperOutliers = sum(x$HillD > q75 + 1.5 * iqr),
      stringsAsFactors = FALSE
    )
  }
)
distribution_diagnostics <- do.call(rbind, distribution_rows)
rownames(distribution_diagnostics) <- NULL
add_check(
  "distribution-diagnostic-rows",
  "diagnostics",
  nrow(distribution_diagnostics),
  9,
  nrow(distribution_diagnostics) == 9L
)
add_check(
  "distribution-no-normality-gate",
  "diagnostics",
  "Shapiro-Wilk not run",
  "Shapiro-Wilk not run",
  !exists("shapiro_results", inherits = FALSE)
)

omnibus_rows <- lapply(0:2, function(order_q) {
  x <- alpha_coverage[alpha_coverage$Order == order_q, , drop = FALSE]
  fit <- stats::kruskal.test(HillD ~ GroupLabel, data = x)
  h <- unname(fit$statistic)
  k <- nlevels(droplevels(x$GroupLabel))
  n <- nrow(x)
  data.frame(
    Order = order_q,
    Metric = paste0("Hill q=", order_q),
    Test = "Kruskal-Wallis rank-sum",
    H = h,
    DF = unname(fit$parameter),
    PValue = fit$p.value,
    EpsilonSquared = max(0, (h - k + 1) / (n - k)),
    Interpretation = "Rank-distribution difference; location only if shapes are comparable",
    stringsAsFactors = FALSE
  )
})
omnibus_tests <- do.call(rbind, omnibus_rows)
omnibus_tests$PAdjustedHolm <- stats::p.adjust(omnibus_tests$PValue, method = "holm")
omnibus_tests$RejectHolm05 <- omnibus_tests$PAdjustedHolm < 0.05
add_check("omnibus-row-count", "inference", nrow(omnibus_tests), 3, nrow(omnibus_tests) == 3L)
add_check(
  "omnibus-holm-family",
  "inference",
  length(omnibus_tests$PAdjustedHolm),
  3,
  length(omnibus_tests$PAdjustedHolm) == 3L
)
add_check(
  "omnibus-values-finite",
  "inference",
  all(is.finite(as.matrix(omnibus_tests[c("H", "DF", "PValue", "EpsilonSquared", "PAdjustedHolm")]))),
  TRUE,
  all(is.finite(as.matrix(omnibus_tests[c("H", "DF", "PValue", "EpsilonSquared", "PAdjustedHolm")])) )
)

welch_rows <- lapply(0:2, function(order_q) {
  x <- alpha_coverage[alpha_coverage$Order == order_q, , drop = FALSE]
  fit <- stats::oneway.test(HillD ~ GroupLabel, data = x, var.equal = FALSE)
  data.frame(
    Order = order_q,
    Metric = paste0("Hill q=", order_q),
    Test = "Welch one-way test",
    Statistic = unname(fit$statistic),
    DF1 = unname(fit$parameter[[1L]]),
    DF2 = unname(fit$parameter[[2L]]),
    PValue = fit$p.value,
    stringsAsFactors = FALSE
  )
})
welch_sensitivity <- do.call(rbind, welch_rows)
welch_sensitivity$PAdjustedHolm <- stats::p.adjust(
  welch_sensitivity$PValue,
  method = "holm"
)
welch_sensitivity$RejectHolm05 <- welch_sensitivity$PAdjustedHolm < 0.05
add_check("welch-row-count", "sensitivity", nrow(welch_sensitivity), 3, nrow(welch_sensitivity) == 3L)

bootstrap_replicates <- 5000L
pair_definitions <- utils::combn(group_levels, 2L, simplify = FALSE)
pairwise_rows <- list()
row_id <- 0L
for (order_q in 0:2) {
  order_data <- alpha_coverage[alpha_coverage$Order == order_q, , drop = FALSE]
  for (pair in pair_definitions) {
    row_id <- row_id + 1L
    x <- order_data$HillD[order_data$GroupLabel == pair[[1L]]]
    y <- order_data$HillD[order_data$GroupLabel == pair[[2L]]]
    wilcox_fit <- suppressWarnings(
      stats::wilcox.test(x, y, exact = FALSE, correct = TRUE)
    )
    delta <- cliffs_delta(x, y)
    delta_ci <- bootstrap_delta(
      x,
      y,
      replicates = bootstrap_replicates,
      seed = 20260721L + 100L * order_q + row_id
    )
    pairwise_rows[[row_id]] <- data.frame(
      Order = order_q,
      Metric = paste0("Hill q=", order_q),
      GroupA = pair[[1L]],
      GroupB = pair[[2L]],
      Contrast = paste(pair[[1L]], "minus", pair[[2L]]),
      NGroupA = length(x),
      NGroupB = length(y),
      MedianGroupA = stats::median(x),
      MedianGroupB = stats::median(y),
      MedianDifference = stats::median(x) - stats::median(y),
      WilcoxonW = unname(wilcox_fit$statistic),
      PValue = wilcox_fit$p.value,
      CliffsDelta = delta,
      ProbabilitySuperiority = (delta + 1) / 2,
      DeltaLower95 = delta_ci[["lower"]],
      DeltaUpper95 = delta_ci[["upper"]],
      BootstrapReplicates = bootstrap_replicates,
      BootstrapValid = as.integer(delta_ci[["valid"]]),
      stringsAsFactors = FALSE
    )
  }
}
pairwise_tests <- do.call(rbind, pairwise_rows)
pairwise_tests$PAdjustedHolm <- stats::p.adjust(
  pairwise_tests$PValue,
  method = "holm"
)
pairwise_tests$RejectHolm05 <- pairwise_tests$PAdjustedHolm < 0.05
add_check("pairwise-row-count", "inference", nrow(pairwise_tests), 9, nrow(pairwise_tests) == 9L)
add_check(
  "pairwise-comparisons-per-order",
  "inference",
  as.integer(table(pairwise_tests$Order)),
  c(3, 3, 3),
  identical(as.integer(table(pairwise_tests$Order)), c(3L, 3L, 3L))
)
add_check(
  "pairwise-global-holm-family",
  "inference",
  length(pairwise_tests$PAdjustedHolm),
  9,
  length(pairwise_tests$PAdjustedHolm) == 9L
)
add_check(
  "pairwise-bootstrap-count",
  "inference",
  unique(pairwise_tests$BootstrapReplicates),
  bootstrap_replicates,
  identical(unique(pairwise_tests$BootstrapReplicates), bootstrap_replicates)
)
add_check(
  "pairwise-bootstrap-valid",
  "inference",
  unique(pairwise_tests$BootstrapValid),
  bootstrap_replicates,
  all(pairwise_tests$BootstrapValid == bootstrap_replicates)
)
add_check(
  "pairwise-delta-range",
  "inference",
  range(pairwise_tests$CliffsDelta),
  "[-1,1]",
  all(abs(pairwise_tests$CliffsDelta) <= 1)
)
add_check(
  "pairwise-ci-ordered",
  "inference",
  all(pairwise_tests$DeltaLower95 <= pairwise_tests$CliffsDelta &
    pairwise_tests$CliffsDelta <= pairwise_tests$DeltaUpper95),
  TRUE,
  all(pairwise_tests$DeltaLower95 <= pairwise_tests$CliffsDelta &
    pairwise_tests$CliffsDelta <= pairwise_tests$DeltaUpper95)
)

standardized_hill <- rbind(
  alpha_raw[c("SampleID", "Group", "GroupLabel", "Order", "HillD", "Standardization")],
  alpha_size[c("SampleID", "Group", "GroupLabel", "Order", "HillD", "Standardization")],
  alpha_coverage[c("SampleID", "Group", "GroupLabel", "Order", "HillD", "Standardization")]
)
standardized_hill$Standardization <- factor(
  standardized_hill$Standardization,
  levels = c("Raw counts", "10,000 reads", "90% coverage")
)
sensitivity_rows <- lapply(
  split(
    standardized_hill,
    interaction(standardized_hill$Standardization, standardized_hill$Order)
  ),
  function(x) {
    fit <- stats::kruskal.test(HillD ~ GroupLabel, data = x)
    h <- unname(fit$statistic)
    k <- nlevels(droplevels(x$GroupLabel))
    n <- nrow(x)
    data.frame(
      Standardization = as.character(unique(x$Standardization)),
      Order = unique(x$Order),
      H = h,
      DF = unname(fit$parameter),
      PValue = fit$p.value,
      EpsilonSquared = max(0, (h - k + 1) / (n - k)),
      stringsAsFactors = FALSE
    )
  }
)
standardization_sensitivity <- do.call(rbind, sensitivity_rows)
rownames(standardization_sensitivity) <- NULL
standardization_sensitivity$Standardization <- factor(
  standardization_sensitivity$Standardization,
  levels = c("Raw counts", "10,000 reads", "90% coverage")
)
standardization_sensitivity <- standardization_sensitivity[
  order(
    standardization_sensitivity$Standardization,
    standardization_sensitivity$Order
  ),
  ,
  drop = FALSE
]
standardization_sensitivity$PAdjustedHolm9 <- stats::p.adjust(
  standardization_sensitivity$PValue,
  method = "holm"
)
standardization_sensitivity$RejectHolm9At05 <-
  standardization_sensitivity$PAdjustedHolm9 < 0.05
add_check(
  "standardization-sensitivity-rows",
  "sensitivity",
  nrow(standardization_sensitivity),
  9,
  nrow(standardization_sensitivity) == 9L
)
raw_q0_p <- standardization_sensitivity$PValue[
  standardization_sensitivity$Standardization == "Raw counts" &
    standardization_sensitivity$Order == 0L
]
size_q0_p <- standardization_sensitivity$PValue[
  standardization_sensitivity$Standardization == "10,000 reads" &
    standardization_sensitivity$Order == 0L
]
coverage_q0_p <- standardization_sensitivity$PValue[
  standardization_sensitivity$Standardization == "90% coverage" &
    standardization_sensitivity$Order == 0L
]
add_check(
  "raw-q0-nominal-signal",
  "sensitivity",
  round(raw_q0_p, 6),
  "<0.05",
  raw_q0_p < 0.05
)
add_check(
  "size-q0-no-nominal-signal",
  "sensitivity",
  round(size_q0_p, 6),
  ">=0.05",
  size_q0_p >= 0.05
)
add_check(
  "coverage-q0-no-nominal-signal",
  "sensitivity",
  round(coverage_q0_p, 6),
  ">=0.05",
  coverage_q0_p >= 0.05
)
add_check(
  "sensitivity-global-holm-none",
  "sensitivity",
  sum(standardization_sensitivity$RejectHolm9At05),
  0,
  sum(standardization_sensitivity$RejectHolm9At05) == 0L
)

design_plot_data <- design_audit
design_plot_data$branch <- factor(
  design_plot_data$branch,
  levels = rev(design_plot_data$branch)
)
design_plot_data$status <- ifelse(
  design_plot_data$selected,
  "Selected",
  "Metadata absent"
)
design_plot_data$dot_fill <- ifelse(
  design_plot_data$selected,
  pal_pub[["green"]],
  "#B8B8B8"
)
design_plot_data$detail <- c(
  "Group: IW / CW / TW; n = 30 each",
  "Needs SubjectID or PairID",
  "Needs SubjectID and Time",
  "Needs a declared clustering column"
)

design_plot <- ggplot2::ggplot(
  design_plot_data,
  ggplot2::aes(y = branch)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0.08, xend = 0.92, yend = branch),
    colour = "#D0D0D0",
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = 0.08, fill = status),
    shape = 21,
    size = 4.2,
    stroke = 0.45,
    colour = "#1A1A1A"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0.15, label = detail),
    hjust = 0,
    size = 3.2,
    colour = "#1A1A1A"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0.92, label = status),
    hjust = 1,
    size = 3.0,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = c("Selected" = pal_pub[["green"]], "Metadata absent" = "#B8B8B8"),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::labs(
    title = "Choose the test from the sampling design",
    subtitle = "Selected branch: independent three-group comparison",
    x = NULL,
    y = NULL,
    caption = "No subject, pair, time, or declared block identifier is present in the demonstration metadata."
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold", size = 8.5),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

order_descriptions <- c(
  "0" = "q = 0 · rare features",
  "1" = "q = 1 · balanced",
  "2" = "q = 2 · dominant features"
)
panel_lookup <- setNames(
  paste0(
    order_descriptions[as.character(omnibus_tests$Order)],
    "\nKruskal-Wallis Holm P ",
    format_p(omnibus_tests$PAdjustedHolm)
  ),
  omnibus_tests$Order
)
distribution_plot_data <- alpha_coverage
distribution_plot_data$Panel <- factor(
  unname(panel_lookup[as.character(distribution_plot_data$Order)]),
  levels = unname(panel_lookup[as.character(0:2)])
)

distribution_plot <- ggplot2::ggplot(
  distribution_plot_data,
  ggplot2::aes(x = GroupLabel, y = HillD, fill = GroupLabel)
) +
  ggplot2::geom_violin(
    width = 0.88,
    alpha = 0.16,
    colour = NA,
    trim = FALSE
  ) +
  ggplot2::geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    alpha = 0.60,
    colour = "#1A1A1A",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = GroupLabel),
    position = ggplot2::position_jitter(width = 0.11, height = 0, seed = 20260721),
    colour = "#1A1A1A",
    stroke = 0.28,
    size = 1.45,
    alpha = 0.68
  ) +
  ggplot2::facet_wrap(~Panel, scales = "free_y", nrow = 1L) +
  ggplot2::scale_fill_manual(values = group_palette, guide = "none") +
  ggplot2::scale_shape_manual(values = group_shapes, guide = "none") +
  ggplot2::scale_y_log10(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Alpha diversity by wetland group",
    subtitle = "Hill diversity standardized to 90% sample coverage; points are independent samples",
    x = NULL,
    y = "Effective features (log scale)",
    caption = "Boxes show medians and IQRs. Kruskal-Wallis tests rank distributions, not medians alone."
  ) +
  theme_pub(base_size = 8.6) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1, size = 7.5),
    strip.text = ggplot2::element_text(size = 7.5),
    legend.position = "none"
  )

sensitivity_plot_data <- standardization_sensitivity
sensitivity_plot_data$OrderLabel <- factor(
  paste0("q = ", sensitivity_plot_data$Order),
  levels = paste0("q = ", 0:2)
)
sensitivity_plot_data$MinusLog10P <- -log10(sensitivity_plot_data$PValue)
sensitivity_plot_data$PLabel <- paste0("P=", format_p(sensitivity_plot_data$PValue))

sensitivity_plot <- ggplot2::ggplot(
  sensitivity_plot_data,
  ggplot2::aes(
    x = Standardization,
    y = MinusLog10P,
    group = OrderLabel,
    colour = OrderLabel
  )
) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    colour = "#666666",
    linewidth = 0.45
  ) +
  ggplot2::geom_line(linewidth = 0.7, alpha = 0.75) +
  ggplot2::geom_point(
    ggplot2::aes(size = EpsilonSquared),
    shape = 21,
    fill = "white",
    stroke = 0.75
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = PLabel),
    nudge_y = 0.10,
    size = 2.55,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "q = 0" = pal_pub[["blue"]],
      "q = 1" = pal_pub[["orange"]],
      "q = 2" = pal_pub[["purple"]]
    )
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.2, 6.2),
    limits = c(0, max(sensitivity_plot_data$EpsilonSquared)),
    name = expression(epsilon^2)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.22)),
    breaks = c(0, -log10(0.05), 2),
    labels = c("1", "0.05", "0.01")
  ) +
  ggplot2::labs(
    title = "The nominal raw-richness signal is not robust to standardization",
    subtitle = "Kruskal-Wallis diagnostics across prespecified measurement scales",
    x = NULL,
    y = "Unadjusted P value (reverse log scale)",
    colour = "Hill order",
    caption = "Point size is epsilon-squared. All nine sensitivity P values are non-significant after global Holm adjustment."
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

forest_plot_data <- pairwise_tests
forest_plot_data$OrderPanel <- factor(
  paste0("q = ", forest_plot_data$Order),
  levels = paste0("q = ", 0:2)
)
forest_plot_data$ContrastLabel <- factor(
  forest_plot_data$Contrast,
  levels = rev(unique(forest_plot_data$Contrast))
)
forest_plot_data$PLabel <- paste0("Holm P=", format_p(forest_plot_data$PAdjustedHolm))

forest_plot <- ggplot2::ggplot(
  forest_plot_data,
  ggplot2::aes(x = CliffsDelta, y = ContrastLabel)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "#666666",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = DeltaLower95, xmax = DeltaUpper95),
    orientation = "y",
    width = 0.16,
    linewidth = 0.6,
    colour = "#3A3A3A"
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 2.6,
    stroke = 0.45,
    fill = pal_pub[["blue"]],
    colour = "#1A1A1A"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = PLabel),
    nudge_y = 0.23,
    size = 2.45,
    colour = "#333333"
  ) +
  ggplot2::facet_wrap(~OrderPanel, nrow = 1L) +
  ggplot2::scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  ggplot2::labs(
    title = "Pairwise dominance effects remain uncertain",
    subtitle = "Cliff's delta with stratified bootstrap 95% percentile intervals",
    x = "Cliff's delta (Group A minus Group B)",
    y = NULL,
    caption = paste(
      "Positive values mean that a random Group A sample tends to have higher diversity than Group B.",
      "Holm adjustment covers all nine planned contrasts.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 7.3)
  )

figure_specification <- data.frame(
  stem = c(
    "20-design-selection-map",
    "20-alpha-group-distributions",
    "20-standardization-sensitivity",
    "20-pairwise-effect-forest"
  ),
  width_mm = c(183, 183, 183, 183),
  height_mm = c(84, 108, 92, 106),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  design_plot,
  distribution_plot,
  sensitivity_plot,
  forest_plot
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
  svg_text <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
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
    add_check(
      paste0(spec$stem, "-", extension, "-exists"),
      "export",
      exists,
      TRUE,
      exists
    )
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
      add_check(
        paste0(spec$stem, "-pdf-single-page"),
        "vector",
        pages,
        1,
        pages == 1L
      )
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
add_check(
  "primary-format-count",
  "export",
  nrow(format_audit),
  16,
  nrow(format_audit) == 16L
)
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
write_tsv(design_audit, "design-audit.tsv")
write_tsv(alpha_raw, "alpha-metrics-raw.tsv")
write_tsv(alpha_size, "alpha-metrics-size-10000.tsv")
write_tsv(alpha_coverage, "alpha-metrics-coverage-090.tsv")
write_tsv(group_descriptive_summary, "group-descriptive-summary.tsv")
write_tsv(distribution_diagnostics, "distribution-diagnostics.tsv")
write_tsv(omnibus_tests, "omnibus-tests.tsv")
write_tsv(pairwise_tests, "pairwise-tests.tsv")
write_tsv(welch_sensitivity, "welch-sensitivity.tsv")
write_tsv(standardization_sensitivity, "standardization-sensitivity.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 20L,
  title = "Alpha diversity group tests and publication graphics",
  dataset = list(
    source = "microeco 2.0.0 wetland 16S example",
    features = ncol(community),
    samples = nrow(community),
    reads = sum(community),
    min_depth = min(library_size),
    median_depth = unname(stats::median(library_size)),
    max_depth = max(library_size)
  ),
  design = list(
    selected_branch = selected_branch,
    groups = length(group_counts),
    samples_per_group = as.integer(group_counts[c("IW", "CW", "TW")]),
    subject_identifier_present = length(present_subject) > 0L,
    pair_identifier_present = length(present_pair) > 0L,
    time_identifier_present = length(present_time) > 0L,
    declared_block_present = length(present_block) > 0L
  ),
  standardization = list(
    primary = "90% sample coverage",
    target_coverage = coverage_target,
    rarefaction_samples = coverage_rarefaction_n,
    extrapolation_samples = coverage_extrapolation_n,
    max_effort_ratio = coverage_max_ratio,
    size_sensitivity_reads = size_target
  ),
  analysis = list(
    primary_metrics = 3L,
    omnibus_tests = nrow(omnibus_tests),
    omnibus_rejections_holm = sum(omnibus_tests$RejectHolm05),
    pairwise_tests = nrow(pairwise_tests),
    pairwise_rejections_holm = sum(pairwise_tests$RejectHolm05),
    bootstrap_replicates = bootstrap_replicates,
    raw_q0_p = unname(raw_q0_p),
    size_q0_p = unname(size_q0_p),
    coverage_q0_p = unname(coverage_q0_p),
    sensitivity_rejections_holm9 = sum(standardization_sensitivity$RejectHolm9At05)
  ),
  graphics = list(
    font_family = font_family,
    primary_figures = length(unique(format_audit$figure)),
    vector_files = sum(format_audit$extension %in% c("pdf", "svg")),
    raster_files = sum(format_audit$extension %in% c("png", "tiff")),
    format_files = nrow(format_audit),
    raster_ppi = 600L,
    svg_files_with_text = sum(
      format_audit$extension == "svg" & format_audit$text_nodes > 0
    ),
    pdf_files_with_embedded_font = sum(
      format_audit$extension == "pdf" & format_audit$font_status == "embedded"
    ),
    lzw_tiff_files = sum(
      format_audit$extension == "tiff" & format_audit$compression == "LZW"
    ),
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "alpha-group-tests-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 10
)

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))
log_lines <- c(
  "Article 20 alpha-group-tests validation",
  paste0("R version: ", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("Dataset: ", ncol(community), " features x ", nrow(community), " samples"),
  paste0("Reads: ", sum(community)),
  paste0("Selected branch: ", selected_branch),
  paste0(
    "Coverage standardization: ", coverage_target,
    "; rarefaction=", coverage_rarefaction_n,
    "; extrapolation=", coverage_extrapolation_n,
    "; max effort ratio=", round(coverage_max_ratio, 6)
  ),
  paste0(
    "Primary omnibus Holm rejections: ",
    sum(omnibus_tests$RejectHolm05),
    "/", nrow(omnibus_tests)
  ),
  paste0(
    "Pairwise Holm rejections: ",
    sum(pairwise_tests$RejectHolm05),
    "/", nrow(pairwise_tests)
  ),
  paste0(
    "Raw / 10k / coverage q0 P: ",
    paste(round(c(raw_q0_p, size_q0_p, coverage_q0_p), 6), collapse = " / ")
  ),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("Checks passed: ", checks_passed, "/", checks_total),
  "",
  paste(
    validation_checks$check_id,
    validation_checks$status,
    validation_checks$observed,
    sep = "\t"
  )
)
writeLines(sanitize_text(log_lines), file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 20 validation failed: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Article 20 validation passed: ",
  checks_passed,
  "/",
  checks_total,
  " checks."
)
