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
  "digest", "ggplot2", "ggrepel", "jsonlite", "ragg", "readr",
  "scales", "svglite", "systemfonts", "vegan"
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
  gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
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

clean_taxon <- function(x, rank_name) {
  x <- trimws(as.character(x))
  x <- sub("^[A-Za-z]__", "", x)
  missing <- is.na(x) | x == "" | tolower(x) %in% c(
    "unassigned", "unclassified", "unknown", "na", "nan"
  )
  x[missing] <- paste("Unknown", tolower(rank_name))
  x
}

jaccard <- function(a, b) {
  union_n <- sum(a | b)
  if (union_n == 0L) 1 else sum(a & b) / union_n
}

group_codes <- c("IW", "CW", "TW")
group_display <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan wetland"
)
rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_prefix <- c(
  Kingdom = "K", Phylum = "P", Class = "C", Order = "O",
  Family = "F", Genus = "G", Species = "S"
)
primary_seed <- 20260728L
primary_detection <- 0.0001
primary_prevalence <- 0.80
primary_rare_threshold <- 0.001
rarefaction_depth <- 10000L
set.seed(primary_seed)

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  source_summary = file.path(input_dir, "source_summary.json"),
  prepare_script = file.path(project_root, "scripts", "prepare_pilot_data.R"),
  theme_pub = file.path(project_root, "R", "theme_pub.R")
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64",
  source_summary = "acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b",
  prepare_script = "fb9ce168b70aeb3636ab73c3c3fabf4208934a04f8eb89a547dc93c7d66a1e3f",
  theme_pub = "8d3821a485aeb529b12184bbe977c1ca952131a22d5eb69b802088ce91909beb"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 28 input(s): ",
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
  ggrepel = "0.9.5",
  jsonlite = "1.8.8",
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
    "core rare occupancy prevalence abundance threshold depth wetland",
    "inland coastal Tibetan detected shared sampling limited 0123456789"
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
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])

counts <- as.matrix(otu_frame)
storage.mode(counts) <- "integer"
feature_ids <- rownames(counts)
sample_ids <- colnames(counts)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
metadata$SampleID <- rownames(metadata)
metadata$GroupDisplay <- unname(group_display[metadata$Group])
metadata$GroupDisplay <- factor(
  metadata$GroupDisplay,
  levels = unname(group_display[group_codes])
)

sample_depths <- colSums(counts)
relative_feature <- sweep(counts, 2L, sample_depths, "/")
feature_equal_mean <- rowMeans(relative_feature)
feature_pooled <- rowSums(counts) / sum(counts)
feature_maximum <- apply(relative_feature, 1L, max)
positive_samples <- rowSums(counts > 0L)
detected_primary <- relative_feature >= primary_detection
detected_samples <- rowSums(detected_primary)
global_prevalence <- detected_samples / ncol(counts)
global_core <- global_prevalence >= primary_prevalence

add_check("feature-count", "data", nrow(counts), 13628, nrow(counts) == 13628L)
add_check("sample-count", "data", ncol(counts), 90, ncol(counts) == 90L)
add_check("read-count", "data", sum(counts), 1619670, sum(counts) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628L)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90L)
add_check(
  "taxonomy-seven-ranks", "data", paste(colnames(taxonomy), collapse = ","),
  paste(rank_names, collapse = ","), identical(colnames(taxonomy), rank_names)
)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "data", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check(
  "feature-id-alignment", "data", sum(feature_ids == rownames(taxonomy)),
  length(feature_ids), identical(feature_ids, rownames(taxonomy))
)
add_check(
  "sample-id-alignment", "data", sum(sample_ids == rownames(metadata)),
  length(sample_ids), identical(sample_ids, rownames(metadata))
)
add_check("counts-finite", "data", sum(is.finite(counts)), length(counts), all(is.finite(counts)))
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0L))
add_check("counts-integer", "data", max(abs(counts - round(counts))), 0, all(counts == round(counts)))
add_check("nonempty-features", "data", sum(rowSums(counts) > 0L), nrow(counts), all(rowSums(counts) > 0L))
add_check("nonempty-samples", "data", sum(sample_depths > 0L), ncol(counts), all(sample_depths > 0L))
add_check(
  "group-levels", "design", paste(sort(unique(metadata$Group)), collapse = ","),
  "CW,IW,TW", identical(sort(unique(metadata$Group)), c("CW", "IW", "TW"))
)
add_check(
  "group-balanced", "design", paste(as.integer(table(metadata$Group)), collapse = ","),
  "30,30,30", all(table(metadata$Group) == 30L)
)
add_check(
  "metadata-columns", "design", paste(colnames(metadata)[1:3], collapse = ","),
  "Group,Type,Saline", identical(colnames(metadata)[1:3], c("Group", "Type", "Saline"))
)
add_check(
  "salinity-visible", "design", length(unique(metadata$Saline)), 2,
  length(unique(metadata$Saline)) == 2L
)
add_check(
  "type-visible", "design", length(unique(metadata$Type)), 6,
  length(unique(metadata$Type)) == 6L
)
add_check("source-summary-samples", "provenance", source_summary$samples, 90, as.numeric(source_summary$samples) == 90)
add_check("source-summary-features", "provenance", source_summary$features, 13628, as.numeric(source_summary$features) == 13628)
add_check("library-depth-min", "data", min(sample_depths), 10364, min(sample_depths) == 10364)
add_check("library-depth-max", "data", max(sample_depths), 37374, max(sample_depths) == 37374)
add_check(
  "feature-relative-closure", "composition",
  max(abs(colSums(relative_feature) - 1)), "<=1e-12",
  max(abs(colSums(relative_feature) - 1)) <= 1e-12
)

sample_alignment_audit <- data.frame(
  SampleID = metadata$SampleID,
  CountTable = metadata$SampleID %in% sample_ids,
  Metadata = metadata$SampleID %in% rownames(metadata),
  GroupCode = metadata$Group,
  Group = as.character(metadata$GroupDisplay),
  Type = metadata$Type,
  Saline = metadata$Saline,
  LibrarySize = as.numeric(sample_depths[metadata$SampleID]),
  stringsAsFactors = FALSE
)

taxonomy_clean <- as.data.frame(
  setNames(
    lapply(rank_names, function(rank_name) clean_taxon(taxonomy[[rank_name]], rank_name)),
    rank_names
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(taxonomy_clean) <- feature_ids
known_matrix <- vapply(
  rank_names,
  function(rank_name) taxonomy_clean[[rank_name]] != paste("Unknown", tolower(rank_name)),
  logical(nrow(taxonomy_clean))
)
best_index <- apply(known_matrix, 1L, function(x) {
  known <- which(x)
  if (length(known) == 0L) 1L else max(known)
})
best_rank <- rank_names[best_index]
best_taxon <- vapply(
  seq_along(feature_ids),
  function(i) taxonomy_clean[[best_rank[[i]]]][[i]],
  character(1)
)
display_label <- paste0(unname(rank_prefix[best_rank]), " · ", best_taxon)
lineage <- apply(taxonomy_clean, 1L, paste, collapse = " > ")

group_core <- list()
group_prevalence <- list()
group_detected_samples <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  group_detected_samples[[group_code]] <- rowSums(detected_primary[, ids, drop = FALSE])
  group_prevalence[[group_code]] <- group_detected_samples[[group_code]] / length(ids)
  group_core[[group_code]] <- group_prevalence[[group_code]] >= primary_prevalence
}

feature_occupancy_abundance <- data.frame(
  FeatureID = feature_ids,
  PositiveSamples = positive_samples,
  RawPresence = positive_samples / ncol(counts),
  DetectedSamples = detected_samples,
  Prevalence = global_prevalence,
  EqualSampleMean = feature_equal_mean,
  PooledReadFraction = feature_pooled,
  MaximumSampleAbundance = feature_maximum,
  TotalReads = rowSums(counts),
  PrimaryCore = global_core,
  IWPrevalence = group_prevalence[["IW"]],
  CWPrevalence = group_prevalence[["CW"]],
  TWPrevalence = group_prevalence[["TW"]],
  IWCore = group_core[["IW"]],
  CWCore = group_core[["CW"]],
  TWCore = group_core[["TW"]],
  BestRank = best_rank,
  BestTaxon = best_taxon,
  DisplayLabel = display_label,
  Lineage = lineage,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_member_rows <- list()
scope_vectors <- c(list(Global = global_core), group_core)
scope_samples <- c(Global = ncol(counts), IW = 30L, CW = 30L, TW = 30L)
scope_detected <- c(list(Global = detected_samples), group_detected_samples)
scope_prevalence <- c(list(Global = global_prevalence), group_prevalence)
for (scope_code in names(scope_vectors)) {
  mask <- scope_vectors[[scope_code]]
  core_member_rows[[scope_code]] <- data.frame(
    ScopeCode = scope_code,
    Scope = ifelse(scope_code == "Global", "All samples", unname(group_display[scope_code])),
    Samples = unname(scope_samples[[scope_code]]),
    FeatureID = feature_ids[mask],
    DetectedSamples = scope_detected[[scope_code]][mask],
    Prevalence = scope_prevalence[[scope_code]][mask],
    EqualSampleMean = feature_equal_mean[mask],
    BestRank = best_rank[mask],
    BestTaxon = best_taxon[mask],
    Lineage = lineage[mask],
    stringsAsFactors = FALSE
  )
}
core_members <- do.call(rbind, core_member_rows)
rownames(core_members) <- NULL

global_core_summary <- data.frame(
  DetectionThreshold = primary_detection,
  PrevalenceThreshold = primary_prevalence,
  CoreFeatures = sum(global_core),
  EqualSampleAbundanceMass = sum(feature_equal_mean[global_core]),
  PooledReadAbundanceMass = sum(feature_pooled[global_core]),
  CoreFeaturesBelowPointOnePercentMean = sum(
    global_core & feature_equal_mean < primary_rare_threshold
  ),
  stringsAsFactors = FALSE
)

group_core_summary <- do.call(rbind, lapply(group_codes, function(group_code) {
  mask <- group_core[[group_code]]
  ids <- metadata$SampleID[metadata$Group == group_code]
  data.frame(
    GroupCode = group_code,
    Group = unname(group_display[group_code]),
    Samples = length(ids),
    CoreFeatures = sum(mask),
    GroupCoreAbundanceMass = mean(colSums(relative_feature[mask, ids, drop = FALSE])),
    GlobalCoreFeaturesInGroupCore = sum(mask & global_core),
    stringsAsFactors = FALSE
  )
}))

membership_flags <- data.frame(
  IW = group_core[["IW"]],
  CW = group_core[["CW"]],
  TW = group_core[["TW"]],
  stringsAsFactors = FALSE
)
membership_key <- apply(membership_flags, 1L, function(x) {
  active <- group_codes[as.logical(x[group_codes])]
  if (length(active) == 0L) "None" else paste(active, collapse = " + ")
})
membership_levels <- c(
  "IW", "CW", "TW", "IW + CW", "IW + TW", "CW + TW", "IW + CW + TW", "None"
)
membership_counts <- table(factor(membership_key, levels = membership_levels))
core_membership_patterns <- data.frame(
  Pattern = membership_levels,
  IW = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
  CW = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE),
  TW = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
  CoreFeatures = as.integer(membership_counts),
  stringsAsFactors = FALSE
)

detection_thresholds <- c(0, 0.0001, 0.0005, 0.001)
detection_labels <- c("Presence", "0.01%", "0.05%", "0.1%")
prevalence_thresholds <- c(0.50, 0.80, 0.90, 1.00)
threshold_rows <- list()
for (detection_i in seq_along(detection_thresholds)) {
  threshold <- detection_thresholds[[detection_i]]
  detected <- if (threshold == 0) counts > 0L else relative_feature >= threshold
  prevalence <- rowMeans(detected)
  for (prevalence_threshold in prevalence_thresholds) {
    mask <- prevalence >= prevalence_threshold
    threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
      DetectionThreshold = threshold,
      DetectionLabel = detection_labels[[detection_i]],
      PrevalenceThreshold = prevalence_threshold,
      PrevalenceLabel = scales::percent(prevalence_threshold, accuracy = 1),
      CoreFeatures = sum(mask),
      EqualSampleAbundanceMass = sum(feature_equal_mean[mask]),
      stringsAsFactors = FALSE
    )
  }
}
core_threshold_sensitivity <- do.call(rbind, threshold_rows)

detection_count_audit <- do.call(rbind, lapply(seq_along(sample_ids), function(i) {
  data.frame(
    SampleID = sample_ids[[i]],
    GroupCode = metadata[sample_ids[[i]], "Group"],
    LibrarySize = sample_depths[[i]],
    DetectionThreshold = detection_thresholds,
    DetectionLabel = detection_labels,
    MinimumReads = ifelse(
      detection_thresholds == 0,
      1L,
      ceiling(sample_depths[[i]] * detection_thresholds)
    ),
    stringsAsFactors = FALSE
  )
}))

set.seed(primary_seed)
rarefied_sample_feature <- vegan::rrarefy(t(counts), sample = rarefaction_depth)
rarefied_counts <- t(rarefied_sample_feature)
set.seed(primary_seed)
rarefied_repeat <- t(vegan::rrarefy(t(counts), sample = rarefaction_depth))
rarefaction_reproducible <- identical(rarefied_counts, rarefied_repeat)
rarefied_one_read_prevalence <- rowMeans(rarefied_counts > 0L)
rarefied_two_read_prevalence <- rowMeans(rarefied_counts >= 2L)
raw_presence_core <- rowMeans(counts > 0L) >= primary_prevalence
rarefied_one_read_core <- rarefied_one_read_prevalence >= primary_prevalence
rarefied_two_read_core <- rarefied_two_read_prevalence >= primary_prevalence
rarefaction_ledger <- data.frame(
  SampleID = sample_ids,
  GroupCode = metadata[sample_ids, "Group"],
  OriginalLibrarySize = as.numeric(sample_depths),
  RarefactionDepth = rarefaction_depth,
  RarefiedLibrarySize = as.numeric(colSums(rarefied_counts)),
  ReadsDiscarded = as.numeric(sample_depths - colSums(rarefied_counts)),
  Seed = primary_seed,
  stringsAsFactors = FALSE
)

depth_definitions <- list(
  `Raw presence` = raw_presence_core,
  `Raw 0.01%` = global_core,
  `Rarefied one read` = rarefied_one_read_core,
  `Rarefied two reads` = rarefied_two_read_core
)
depth_sensitivity <- do.call(rbind, lapply(names(depth_definitions), function(definition) {
  mask <- depth_definitions[[definition]]
  data.frame(
    Definition = definition,
    CoreFeatures = sum(mask),
    SharedWithPrimary = sum(mask & global_core),
    UnionWithPrimary = sum(mask | global_core),
    JaccardWithPrimary = jaccard(mask, global_core),
    EqualSampleAbundanceMass = sum(feature_equal_mean[mask]),
    stringsAsFactors = FALSE
  )
}))

feature_state <- ifelse(
  global_core,
  "Primary core",
  ifelse(
    rowSums(counts) <= 10L & positive_samples <= 2L,
    "Sampling-limited",
    ifelse(
      feature_equal_mean < primary_rare_threshold & feature_maximum >= 0.01,
      "Conditionally rare",
      ifelse(feature_maximum < primary_rare_threshold, "Persistently rare", "Intermediate")
    )
  )
)
state_levels <- c(
  "Primary core", "Conditionally rare", "Persistently rare",
  "Sampling-limited", "Intermediate"
)
feature_state <- factor(feature_state, levels = state_levels)
feature_occupancy_abundance$FeatureState <- as.character(feature_state)
feature_state_classification <- data.frame(
  feature_occupancy_abundance[, c(
    "FeatureID", "FeatureState", "TotalReads", "PositiveSamples", "RawPresence",
    "DetectedSamples", "Prevalence", "EqualSampleMean", "PooledReadFraction",
    "MaximumSampleAbundance", "BestRank", "BestTaxon", "Lineage"
  )],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
feature_state_summary <- do.call(rbind, lapply(state_levels, function(state) {
  mask <- feature_state == state
  data.frame(
    FeatureState = state,
    Features = sum(mask),
    FeatureFraction = mean(mask),
    EqualSampleAbundanceMass = sum(feature_equal_mean[mask]),
    PooledReadAbundanceMass = sum(feature_pooled[mask]),
    stringsAsFactors = FALSE
  )
}))

state_relative <- rowsum(relative_feature, group = feature_state, reorder = FALSE)
state_relative <- state_relative[state_levels, , drop = FALSE]
group_state_mass <- do.call(rbind, lapply(group_codes, function(group_code) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  data.frame(
    GroupCode = group_code,
    Group = unname(group_display[group_code]),
    FeatureState = state_levels,
    MeanAbundanceMass = rowMeans(state_relative[, ids, drop = FALSE]),
    MinimumSampleMass = apply(state_relative[, ids, drop = FALSE], 1L, min),
    MaximumSampleMass = apply(state_relative[, ids, drop = FALSE], 1L, max),
    stringsAsFactors = FALSE
  )
}))

rare_thresholds <- c(0.0001, 0.0005, 0.001, 0.005)
rare_threshold_labels <- c("0.01%", "0.05%", "0.1%", "0.5%")
rare_sample_rows <- list()
for (threshold_i in seq_along(rare_thresholds)) {
  threshold <- rare_thresholds[[threshold_i]]
  rare_mass <- colSums(relative_feature * (relative_feature > 0 & relative_feature < threshold))
  rare_sample_rows[[threshold_i]] <- data.frame(
    SampleID = sample_ids,
    GroupCode = metadata[sample_ids, "Group"],
    Group = unname(group_display[metadata[sample_ids, "Group"]]),
    Type = metadata[sample_ids, "Type"],
    Saline = metadata[sample_ids, "Saline"],
    RareThreshold = threshold,
    ThresholdLabel = rare_threshold_labels[[threshold_i]],
    RareMass = rare_mass,
    stringsAsFactors = FALSE
  )
}
rare_mass_all_thresholds <- do.call(rbind, rare_sample_rows)
rare_mass_by_sample <- rare_mass_all_thresholds[
  rare_mass_all_thresholds$RareThreshold == primary_rare_threshold,
  , drop = FALSE
]
rare_threshold_sensitivity <- do.call(rbind, lapply(seq_along(rare_thresholds), function(i) {
  threshold_data <- rare_mass_all_thresholds[
    rare_mass_all_thresholds$RareThreshold == rare_thresholds[[i]],
    , drop = FALSE
  ]
  rows <- list(data.frame(
    GroupCode = "ALL",
    Group = "All samples",
    RareThreshold = rare_thresholds[[i]],
    ThresholdLabel = rare_threshold_labels[[i]],
    Samples = nrow(threshold_data),
    MeanRareMass = mean(threshold_data$RareMass),
    MinimumRareMass = min(threshold_data$RareMass),
    MaximumRareMass = max(threshold_data$RareMass),
    stringsAsFactors = FALSE
  ))
  for (group_code in group_codes) {
    values <- threshold_data$RareMass[threshold_data$GroupCode == group_code]
    rows[[length(rows) + 1L]] <- data.frame(
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      RareThreshold = rare_thresholds[[i]],
      ThresholdLabel = rare_threshold_labels[[i]],
      Samples = length(values),
      MeanRareMass = mean(values),
      MinimumRareMass = min(values),
      MaximumRareMass = max(values),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}))

core_feature_taxonomy <- data.frame(
  FeatureID = feature_ids[global_core],
  taxonomy_clean[global_core, , drop = FALSE],
  DetectedSamples = detected_samples[global_core],
  Prevalence = global_prevalence[global_core],
  EqualSampleMean = feature_equal_mean[global_core],
  PooledReadFraction = feature_pooled[global_core],
  MaximumSampleAbundance = feature_maximum[global_core],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
core_feature_taxonomy <- core_feature_taxonomy[
  order(core_feature_taxonomy$EqualSampleMean, decreasing = TRUE),
  , drop = FALSE
]

closure_audit <- data.frame(
  SampleID = sample_ids,
  GroupCode = metadata[sample_ids, "Group"],
  RelativeAbundanceClosure = colSums(relative_feature),
  FeatureStateClosure = colSums(state_relative),
  PrimaryRareMass = rare_mass_by_sample$RareMass[match(sample_ids, rare_mass_by_sample$SampleID)],
  stringsAsFactors = FALSE
)
closure_audit$MaximumAbsoluteClosureError <- pmax(
  abs(closure_audit$RelativeAbundanceClosure - 1),
  abs(closure_audit$FeatureStateClosure - 1)
)
for (i in seq_len(nrow(closure_audit))) {
  add_check(
    paste0("sample-closure-", closure_audit$SampleID[[i]]),
    "composition",
    closure_audit$MaximumAbsoluteClosureError[[i]],
    "<=1e-12",
    closure_audit$MaximumAbsoluteClosureError[[i]] <= 1e-12
  )
}

expected_threshold_counts <- matrix(
  c(
    1499L, 334L, 166L, 27L,
    404L, 97L, 53L, 10L,
    48L, 9L, 5L, 0L,
    16L, 0L, 0L, 0L
  ),
  nrow = 4L,
  byrow = TRUE,
  dimnames = list(detection_labels, scales::percent(prevalence_thresholds, accuracy = 1))
)
observed_threshold_counts <- matrix(
  core_threshold_sensitivity$CoreFeatures,
  nrow = 4L,
  byrow = TRUE,
  dimnames = dimnames(expected_threshold_counts)
)
for (detection_label in detection_labels) {
  for (prevalence_label in colnames(expected_threshold_counts)) {
    observed <- observed_threshold_counts[detection_label, prevalence_label]
    expected <- expected_threshold_counts[detection_label, prevalence_label]
    add_check(
      paste0(
        "threshold-", gsub("[^a-z0-9]+", "-", tolower(detection_label)), "-",
        gsub("%", "", prevalence_label, fixed = TRUE)
      ),
      "sensitivity", observed, expected, observed == expected
    )
  }
}

expected_group_core <- c(IW = 166L, CW = 160L, TW = 278L)
for (group_code in group_codes) {
  observed <- group_core_summary$CoreFeatures[group_core_summary$GroupCode == group_code]
  add_check(
    paste0("group-core-", tolower(group_code)), "core",
    observed, expected_group_core[[group_code]], observed == expected_group_core[[group_code]]
  )
}
expected_pattern_counts <- c(
  IW = 79L,
  CW = 60L,
  TW = 166L,
  `IW + CW` = 8L,
  `IW + TW` = 20L,
  `CW + TW` = 33L,
  `IW + CW + TW` = 59L,
  None = 13203L
)
for (pattern in names(expected_pattern_counts)) {
  observed <- core_membership_patterns$CoreFeatures[
    core_membership_patterns$Pattern == pattern
  ]
  add_check(
    paste0("membership-", gsub("[^a-z0-9]+", "-", tolower(pattern))),
    "core", observed, expected_pattern_counts[[pattern]],
    observed == expected_pattern_counts[[pattern]]
  )
}

shared_three <- sum(group_core[["IW"]] & group_core[["CW"]] & group_core[["TW"]])
group_union <- sum(group_core[["IW"]] | group_core[["CW"]] | group_core[["TW"]])
add_check("global-core-count", "core", sum(global_core), 97, sum(global_core) == 97L)
add_check(
  "global-core-abundance-mass", "core", sum(feature_equal_mean[global_core]),
  0.170361211031733, near(sum(feature_equal_mean[global_core]), 0.170361211031733, 1e-12)
)
add_check(
  "core-below-point-one-percent", "core",
  global_core_summary$CoreFeaturesBelowPointOnePercentMean, 34,
  global_core_summary$CoreFeaturesBelowPointOnePercentMean == 34L
)
add_check("shared-three-group-core", "core", shared_three, 59, shared_three == 59L)
add_check("group-core-union", "core", group_union, 425, group_union == 425L)
add_check(
  "global-versus-shared-distinct", "core", sum(global_core) - shared_three, 38,
  sum(global_core) - shared_three == 38L
)
add_check(
  "primary-detection-read-min", "depth",
  min(detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]),
  2, min(detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]) == 2L
)
add_check(
  "primary-detection-read-max", "depth",
  max(detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]),
  4, max(detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]) == 4L
)
add_check(
  "rarefaction-reproducible", "depth", rarefaction_reproducible, TRUE,
  rarefaction_reproducible
)
add_check(
  "rarefaction-all-depths", "depth", paste(unique(colSums(rarefied_counts)), collapse = ","),
  rarefaction_depth, all(colSums(rarefied_counts) == rarefaction_depth)
)
expected_depth_counts <- c(
  `Raw presence` = 334L,
  `Raw 0.01%` = 97L,
  `Rarefied one read` = 167L,
  `Rarefied two reads` = 52L
)
for (definition in names(expected_depth_counts)) {
  observed <- depth_sensitivity$CoreFeatures[depth_sensitivity$Definition == definition]
  add_check(
    paste0("depth-core-", gsub("[^a-z0-9]+", "-", tolower(definition))),
    "depth", observed, expected_depth_counts[[definition]], observed == expected_depth_counts[[definition]]
  )
}
primary_rarefied_jaccard <- jaccard(global_core, rarefied_one_read_core)
add_check(
  "primary-rarefied-one-read-jaccard", "depth", primary_rarefied_jaccard,
  0.580838323353293, near(primary_rarefied_jaccard, 0.580838323353293, 1e-12)
)
add_check(
  "primary-features-retained-rarefied-one-read", "depth",
  sum(global_core & rarefied_one_read_core), 97,
  sum(global_core & rarefied_one_read_core) == 97L
)

expected_state_counts <- c(
  `Primary core` = 97L,
  `Conditionally rare` = 200L,
  `Persistently rare` = 8407L,
  `Sampling-limited` = 1155L,
  Intermediate = 3769L
)
expected_state_mass <- c(
  `Primary core` = 0.170361211031733,
  `Conditionally rare` = 0.0938860476416336,
  `Persistently rare` = 0.172682530092607,
  `Sampling-limited` = 0.00158394392112484,
  Intermediate = 0.561486267312901
)
for (state in state_levels) {
  observed_count <- feature_state_summary$Features[feature_state_summary$FeatureState == state]
  observed_mass <- feature_state_summary$EqualSampleAbundanceMass[
    feature_state_summary$FeatureState == state
  ]
  add_check(
    paste0("state-count-", gsub("[^a-z0-9]+", "-", tolower(state))),
    "rare", observed_count, expected_state_counts[[state]],
    observed_count == expected_state_counts[[state]]
  )
  add_check(
    paste0("state-mass-", gsub("[^a-z0-9]+", "-", tolower(state))),
    "rare", observed_mass, expected_state_mass[[state]],
    near(observed_mass, expected_state_mass[[state]], 1e-12)
  )
}
add_check(
  "state-features-exhaustive", "rare", sum(feature_state_summary$Features), nrow(counts),
  sum(feature_state_summary$Features) == nrow(counts)
)
add_check(
  "state-mass-closure", "rare", sum(feature_state_summary$EqualSampleAbundanceMass), 1,
  near(sum(feature_state_summary$EqualSampleAbundanceMass), 1, 1e-12)
)

expected_rare_means <- c(
  `0.01%` = 0.104676641055023,
  `0.05%` = 0.318640788866651,
  `0.1%` = 0.441501476980241,
  `0.5%` = 0.736507158030592
)
for (threshold_label in names(expected_rare_means)) {
  observed <- rare_threshold_sensitivity$MeanRareMass[
    rare_threshold_sensitivity$GroupCode == "ALL" &
      rare_threshold_sensitivity$ThresholdLabel == threshold_label
  ]
  add_check(
    paste0("rare-mass-", gsub("%", "", threshold_label, fixed = TRUE)),
    "rare", observed, expected_rare_means[[threshold_label]],
    near(observed, expected_rare_means[[threshold_label]], 1e-12)
  )
}
add_check(
  "primary-rare-mass-min", "rare", min(rare_mass_by_sample$RareMass),
  0.233295332193299, near(min(rare_mass_by_sample$RareMass), 0.233295332193299, 1e-12)
)
add_check(
  "primary-rare-mass-max", "rare", max(rare_mass_by_sample$RareMass),
  0.594348466586313, near(max(rare_mass_by_sample$RareMass), 0.594348466586313, 1e-12)
)
add_check(
  "rare-threshold-monotone", "rare",
  paste(expected_rare_means, collapse = ","), "strictly increasing",
  all(diff(expected_rare_means) > 0)
)
add_check(
  "group-state-closure", "composition",
  max(abs(tapply(group_state_mass$MeanAbundanceMass, group_state_mass$GroupCode, sum) - 1)),
  "<=1e-12",
  max(abs(tapply(group_state_mass$MeanAbundanceMass, group_state_mass$GroupCode, sum) - 1)) <= 1e-12
)

state_colours <- c(
  `Primary core` = "#0072B2",
  `Conditionally rare` = "#D55E00",
  `Persistently rare` = "#56B4E9",
  `Sampling-limited` = "#9CA3AF",
  Intermediate = "#CC79A7"
)
group_colours <- c(
  `Inland wetland` = "#E69F00",
  `Coastal wetland` = "#0072B2",
  `Tibetan wetland` = "#009E73"
)

occupancy_plot_data <- feature_occupancy_abundance
occupancy_plot_data$FeatureState <- factor(
  occupancy_plot_data$FeatureState,
  levels = state_levels
)
core_label_ids <- head(
  feature_ids[global_core][order(feature_equal_mean[global_core], decreasing = TRUE)],
  5L
)
conditional_mask <- feature_state == "Conditionally rare"
conditional_label_ids <- head(
  feature_ids[conditional_mask][order(feature_maximum[conditional_mask], decreasing = TRUE)],
  4L
)
occupancy_labels <- occupancy_plot_data[
  occupancy_plot_data$FeatureID %in% c(core_label_ids, conditional_label_ids),
  , drop = FALSE
]
occupancy_labels$PointLabel <- paste0(
  occupancy_labels$DisplayLabel,
  " (",
  occupancy_labels$FeatureID,
  ")"
)
occupancy_plot <- ggplot2::ggplot(
  occupancy_plot_data,
  ggplot2::aes(
    x = Prevalence,
    y = EqualSampleMean,
    colour = FeatureState
  )
) +
  ggplot2::geom_hline(
    yintercept = primary_rare_threshold,
    colour = "#6B7280",
    linewidth = 0.45,
    linetype = "22"
  ) +
  ggplot2::geom_vline(
    xintercept = primary_prevalence,
    colour = "#1F2937",
    linewidth = 0.55,
    linetype = "22"
  ) +
  ggplot2::geom_point(size = 1.15, alpha = 0.58) +
  ggrepel::geom_text_repel(
    data = occupancy_labels,
    ggplot2::aes(label = PointLabel),
    seed = primary_seed,
    family = font_family,
    size = 2.25,
    colour = "#17202A",
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = state_colours, drop = FALSE) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.005, 0.01))
  ) +
  ggplot2::scale_y_log10(
    labels = scales::label_percent(accuracy = 0.001),
    expand = ggplot2::expansion(mult = c(0.04, 0.16))
  ) +
  ggplot2::annotation_logticks(sides = "l", linewidth = 0.25) +
  ggplot2::labs(
    title = "Core membership depends on an explicit abundance–occupancy rule",
    subtitle = "Primary core: at least 0.01% abundance in at least 80% of 90 samples",
    x = "Prevalence at the 0.01% detection threshold",
    y = "Equal-sample mean relative abundance",
    colour = "Feature state",
    caption = paste(
      "Each point is one feature/OTU; all-read denominators and equal sample weights are retained.",
      "The horizontal guide marks 0.1% mean abundance, not a core boundary."
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

membership_plot_data <- core_membership_patterns[
  core_membership_patterns$Pattern != "None",
  , drop = FALSE
]
membership_plot_data <- membership_plot_data[
  order(membership_plot_data$CoreFeatures, decreasing = TRUE),
  , drop = FALSE
]
membership_plot_data$X <- seq_len(nrow(membership_plot_data))
membership_plot_data$PatternLabel <- gsub(" \\+ ", " · ", membership_plot_data$Pattern)
membership_y <- c(IW = -18, CW = -39, TW = -60)
membership_grid <- do.call(rbind, lapply(seq_len(nrow(membership_plot_data)), function(i) {
  data.frame(
    X = membership_plot_data$X[[i]],
    GroupCode = group_codes,
    Y = unname(membership_y[group_codes]),
    Active = as.logical(membership_plot_data[i, group_codes]),
    stringsAsFactors = FALSE
  )
}))
membership_segments <- do.call(rbind, lapply(seq_len(nrow(membership_plot_data)), function(i) {
  active_y <- membership_y[group_codes[as.logical(membership_plot_data[i, group_codes])]]
  data.frame(
    X = membership_plot_data$X[[i]],
    YMin = min(active_y),
    YMax = max(active_y),
    stringsAsFactors = FALSE
  )
}))
membership_plot <- ggplot2::ggplot(membership_plot_data, ggplot2::aes(x = X)) +
  ggplot2::geom_col(
    ggplot2::aes(y = CoreFeatures),
    width = 0.68,
    fill = "#0072B2"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = CoreFeatures, label = CoreFeatures),
    vjust = -0.35,
    family = font_family,
    size = 2.8
  ) +
  ggplot2::geom_segment(
    data = membership_segments,
    ggplot2::aes(x = X, xend = X, y = YMin, yend = YMax),
    inherit.aes = FALSE,
    linewidth = 0.75,
    colour = "#374151"
  ) +
  ggplot2::geom_point(
    data = membership_grid,
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    size = 3.1,
    shape = 21,
    stroke = 0.45,
    colour = "#B8C0CC",
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = membership_grid[membership_grid$Active, , drop = FALSE],
    ggplot2::aes(x = X, y = Y),
    inherit.aes = FALSE,
    size = 3.1,
    colour = "#17202A"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.38,
    y = unname(membership_y),
    label = group_codes,
    hjust = 0,
    family = font_family,
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::scale_x_continuous(
    breaks = membership_plot_data$X,
    labels = membership_plot_data$PatternLabel,
    limits = c(0.3, nrow(membership_plot_data) + 0.55),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(0, 50, 100, 150),
    labels = c("0", "50", "100", "150"),
    limits = c(-72, 190),
    expand = ggplot2::expansion(mult = c(0, 0.015))
  ) +
  ggplot2::labs(
    title = "Group cores overlap, but they are not the global core",
    subtitle = "The same 0.01% detection and 80% prevalence rules are applied within each 30-sample group",
    x = "Group-core membership pattern",
    y = "Features",
    caption = paste(
      "The seven non-empty intersections contain 425 features; 59 are core in all three groups.",
      "IW, inland; CW, coastal; TW, Tibetan wetland."
    )
  ) +
  theme_pub(base_size = 8.5) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1),
    axis.ticks.y = ggplot2::element_blank()
  )

fixed_depth_rows <- list()
for (prevalence_threshold in prevalence_thresholds) {
  fixed_depth_rows[[length(fixed_depth_rows) + 1L]] <- data.frame(
    DetectionLabel = c("One read", "Two reads"),
    PrevalenceThreshold = prevalence_threshold,
    PrevalenceLabel = scales::percent(prevalence_threshold, accuracy = 1),
    CoreFeatures = c(
      sum(rarefied_one_read_prevalence >= prevalence_threshold),
      sum(rarefied_two_read_prevalence >= prevalence_threshold)
    ),
    stringsAsFactors = FALSE
  )
}
fixed_depth_sensitivity <- do.call(rbind, fixed_depth_rows)
threshold_plot_data <- rbind(
  data.frame(
    Panel = "Raw libraries",
    DetectionLabel = core_threshold_sensitivity$DetectionLabel,
    PrevalenceLabel = core_threshold_sensitivity$PrevalenceLabel,
    CoreFeatures = core_threshold_sensitivity$CoreFeatures,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "10,000-read branch",
    DetectionLabel = fixed_depth_sensitivity$DetectionLabel,
    PrevalenceLabel = fixed_depth_sensitivity$PrevalenceLabel,
    CoreFeatures = fixed_depth_sensitivity$CoreFeatures,
    stringsAsFactors = FALSE
  )
)
threshold_plot_data$DetectionLabel <- factor(
  threshold_plot_data$DetectionLabel,
  levels = c(detection_labels, "One read", "Two reads")
)
threshold_plot_data$PrevalenceLabel <- factor(
  threshold_plot_data$PrevalenceLabel,
  levels = scales::percent(prevalence_thresholds, accuracy = 1)
)
threshold_plot_data$Panel <- factor(
  threshold_plot_data$Panel,
  levels = c("Raw libraries", "10,000-read branch")
)
primary_tiles <- threshold_plot_data[
  (threshold_plot_data$Panel == "Raw libraries" &
    threshold_plot_data$DetectionLabel == "0.01%" &
    threshold_plot_data$PrevalenceLabel == "80%") |
    (threshold_plot_data$Panel == "10,000-read branch" &
      threshold_plot_data$DetectionLabel == "One read" &
      threshold_plot_data$PrevalenceLabel == "80%"),
  , drop = FALSE
]
threshold_plot <- ggplot2::ggplot(
  threshold_plot_data,
  ggplot2::aes(x = DetectionLabel, y = PrevalenceLabel, fill = CoreFeatures)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_tile(
    data = primary_tiles,
    fill = NA,
    colour = "#D55E00",
    linewidth = 1.15
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(CoreFeatures)),
    family = font_family,
    size = 3
  ) +
  ggplot2::facet_grid(
    . ~ Panel,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    trans = "sqrt",
    labels = scales::label_comma()
  ) +
  ggplot2::labs(
    title = "Core size changes with detection, prevalence and sequencing depth",
    subtitle = "Raw thresholds use all reads; the fixed-depth branch uses seed 20260728 and is a sensitivity analysis",
    x = "Detection rule",
    y = "Minimum prevalence",
    fill = "Core features",
    caption = paste0(
      "Orange outlines compare the primary raw core (97) with the rarefied one-read core (167) at 80%; ",
      "Jaccard = ", sprintf("%.3f", primary_rarefied_jaccard), "."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1),
    strip.text = ggplot2::element_text(size = 9)
  )

rare_plot_data <- rare_mass_all_thresholds
rare_plot_data$ThresholdLabel <- factor(
  rare_plot_data$ThresholdLabel,
  levels = rare_threshold_labels
)
rare_plot_data$Group <- factor(
  rare_plot_data$Group,
  levels = unname(group_display[group_codes])
)
rare_group_means <- rare_threshold_sensitivity[
  rare_threshold_sensitivity$GroupCode %in% group_codes,
  , drop = FALSE
]
rare_group_means$ThresholdLabel <- factor(
  rare_group_means$ThresholdLabel,
  levels = rare_threshold_labels
)
rare_group_means$Group <- factor(
  rare_group_means$Group,
  levels = unname(group_display[group_codes])
)
rare_plot <- ggplot2::ggplot(
  rare_plot_data,
  ggplot2::aes(x = ThresholdLabel, y = RareMass, fill = Group, colour = Group)
) +
  ggplot2::annotate(
    "rect",
    xmin = 2.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F0E442",
    alpha = 0.10
  ) +
  ggplot2::geom_boxplot(
    width = 0.68,
    outlier.shape = NA,
    alpha = 0.35,
    position = ggplot2::position_dodge(width = 0.78),
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    size = 0.85,
    alpha = 0.30,
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.11,
      dodge.width = 0.78,
      seed = primary_seed
    )
  ) +
  ggplot2::geom_line(
    data = rare_group_means,
    ggplot2::aes(y = MeanRareMass, group = Group),
    position = ggplot2::position_dodge(width = 0.78),
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = rare_group_means,
    ggplot2::aes(y = MeanRareMass),
    position = ggplot2::position_dodge(width = 0.78),
    size = 2.2,
    shape = 21,
    stroke = 0.7
  ) +
  ggplot2::scale_fill_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 0.93)) +
  ggplot2::labs(
    title = "The measured rare-biosphere mass is threshold dependent",
    subtitle = "Local rare mass sums features below the selected abundance threshold within each sample",
    x = "Local rare-abundance threshold",
    y = "Relative-abundance mass below threshold",
    fill = "Wetland group",
    colour = "Wetland group",
    caption = paste(
      "The highlighted 0.1% rule gives a 44.2% equal-sample mean (range 23.3–59.4%).",
      "Points are samples; lines connect descriptive group means without hypothesis tests."
    )
  ) +
  theme_pub(base_size = 8.8) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major.x = ggplot2::element_blank()
  )

figure_specification <- data.frame(
  stem = c(
    "28-occupancy-abundance",
    "28-group-core-membership",
    "28-threshold-depth-sensitivity",
    "28-rare-biosphere-mass"
  ),
  width_mm = c(183, 165, 170, 170),
  height_mm = c(132, 118, 112, 116),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  occupancy_plot,
  membership_plot,
  threshold_plot,
  rare_plot
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
  sum(format_audit$chinese_characters[format_audit$extension == "svg"], na.rm = TRUE),
  0,
  all(format_audit$chinese_characters[format_audit$extension == "svg"] == 0L)
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
write_tsv(feature_occupancy_abundance, "feature-occupancy-abundance.tsv")
write_tsv(core_members, "core-members.tsv")
write_tsv(global_core_summary, "global-core-summary.tsv")
write_tsv(group_core_summary, "group-core-summary.tsv")
write_tsv(core_membership_patterns, "core-membership-patterns.tsv")
write_tsv(core_threshold_sensitivity, "core-threshold-sensitivity.tsv")
write_tsv(detection_count_audit, "detection-count-audit.tsv")
write_tsv(rarefaction_ledger, "rarefaction-ledger.tsv")
write_tsv(depth_sensitivity, "depth-sensitivity.tsv")
write_tsv(feature_state_classification, "feature-state-classification.tsv")
write_tsv(feature_state_summary, "feature-state-summary.tsv")
write_tsv(rare_mass_by_sample, "rare-mass-by-sample.tsv")
write_tsv(rare_threshold_sensitivity, "rare-threshold-sensitivity.tsv")
write_tsv(group_state_mass, "group-state-mass.tsv")
write_tsv(core_feature_taxonomy, "core-feature-taxonomy.tsv")
write_tsv(closure_audit, "closure-audit.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 28L,
  title = "Core microbiome and the rare biosphere",
  dataset = list(
    source = "microeco 2.0.0 packaged Chinese wetland 16S tables",
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    groups = length(group_codes),
    samples_per_group = as.list(setNames(as.integer(table(metadata$Group)[group_codes]), group_codes)),
    library_depth_min = min(sample_depths),
    library_depth_max = max(sample_depths)
  ),
  primary_definition = list(
    feature_level = "OTU/feature",
    denominator = "All reads within each sample",
    sample_weighting = "Equal sample weight",
    detection_threshold = primary_detection,
    detection_boundary = "relative abundance >= threshold",
    prevalence_threshold = primary_prevalence,
    prevalence_boundary = "prevalence >= threshold",
    seed = primary_seed
  ),
  core = list(
    global_features = sum(global_core),
    global_equal_sample_abundance_mass = sum(feature_equal_mean[global_core]),
    global_features_below_point_one_percent_mean = global_core_summary$CoreFeaturesBelowPointOnePercentMean,
    group_features = as.list(setNames(group_core_summary$CoreFeatures, group_core_summary$GroupCode)),
    shared_three_group_features = shared_three,
    group_union_features = group_union,
    global_minus_shared_three = sum(global_core) - shared_three
  ),
  sensitivity = list(
    threshold_grid_rows = nrow(core_threshold_sensitivity),
    detection_thresholds = detection_thresholds,
    prevalence_thresholds = prevalence_thresholds,
    primary_minimum_reads_min = min(
      detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]
    ),
    primary_minimum_reads_max = max(
      detection_count_audit$MinimumReads[detection_count_audit$DetectionLabel == "0.01%"]
    ),
    rarefaction_depth = rarefaction_depth,
    rarefaction_seed = primary_seed,
    rarefaction_reproducible = rarefaction_reproducible,
    raw_presence_core_features = sum(raw_presence_core),
    rarefied_one_read_core_features = sum(rarefied_one_read_core),
    rarefied_two_read_core_features = sum(rarefied_two_read_core),
    primary_rarefied_one_read_jaccard = primary_rarefied_jaccard
  ),
  rare = list(
    local_definition = "0 < per-sample relative abundance < 0.001",
    local_rare_mass_mean = mean(rare_mass_by_sample$RareMass),
    local_rare_mass_minimum = min(rare_mass_by_sample$RareMass),
    local_rare_mass_maximum = max(rare_mass_by_sample$RareMass),
    classified_features = sum(feature_state_summary$Features),
    primary_core_features = feature_state_summary$Features[
      feature_state_summary$FeatureState == "Primary core"
    ],
    conditionally_rare_features = feature_state_summary$Features[
      feature_state_summary$FeatureState == "Conditionally rare"
    ],
    persistently_rare_features = feature_state_summary$Features[
      feature_state_summary$FeatureState == "Persistently rare"
    ],
    sampling_limited_features = feature_state_summary$Features[
      feature_state_summary$FeatureState == "Sampling-limited"
    ],
    intermediate_features = feature_state_summary$Features[
      feature_state_summary$FeatureState == "Intermediate"
    ],
    rare_threshold_means = as.list(setNames(
      rare_threshold_sensitivity$MeanRareMass[rare_threshold_sensitivity$GroupCode == "ALL"],
      rare_threshold_sensitivity$ThresholdLabel[rare_threshold_sensitivity$GroupCode == "ALL"]
    ))
  ),
  inference = list(
    group_hypothesis_test = FALSE,
    group_and_salinity_confounded = TRUE,
    absolute_abundance_claim = FALSE,
    activity_claim = FALSE,
    function_claim = FALSE
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
      format_audit$chinese_characters[format_audit$extension == "svg"],
      na.rm = TRUE
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
  file.path(output_dir, "core-rare-biosphere-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

validation_lines <- c(
  "Article 28 core/rare-biosphere validation",
  paste0("Status: ", ifelse(checks_failed == 0L, "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total),
  "Samples/features/reads: 90/13628/1619670",
  paste0("Primary core features: ", sum(global_core)),
  paste0("Primary core abundance mass: ", sprintf("%.6f", sum(feature_equal_mean[global_core]))),
  paste0("Group core IW/CW/TW: ", paste(group_core_summary$CoreFeatures, collapse = "/")),
  paste0("Shared three-group core: ", shared_three),
  paste0("Group-core union: ", group_union),
  paste0("Raw/rarefied-one-read core: ", sum(raw_presence_core), "/", sum(rarefied_one_read_core)),
  paste0("Primary/rarefied Jaccard: ", sprintf("%.6f", primary_rarefied_jaccard)),
  paste0("Primary local rare mass: ", sprintf("%.6f", mean(rare_mass_by_sample$RareMass))),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("R version: ", r_version)
)
writeLines(validation_lines, file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed_ids <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 28 validation failed: ",
    paste(failed_ids, collapse = ", "),
    call. = FALSE
  )
}

cat(
  sprintf(
    "Article 28 validation passed: %d/%d checks.\n",
    checks_passed,
    checks_total
  )
)
