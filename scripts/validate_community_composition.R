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
  "digest", "ggplot2", "jsonlite", "ragg", "readr", "scales",
  "svglite", "systemfonts"
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

aggregate_rank <- function(counts, taxonomy, rank_name) {
  labels <- clean_taxon(taxonomy[[rank_name]], rank_name)
  out <- rowsum(counts, group = labels, reorder = TRUE)
  out[order(rownames(out)), , drop = FALSE]
}

matrix_to_long <- function(x, value_name) {
  out <- data.frame(
    Taxon = rep(rownames(x), times = ncol(x)),
    SampleID = rep(colnames(x), each = nrow(x)),
    value = as.vector(x),
    stringsAsFactors = FALSE
  )
  names(out)[[3L]] <- value_name
  out
}

group_display <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan wetland"
)
type_display <- c(
  NE = "Northeast",
  NW = "Northwest",
  NC = "North China",
  YML = "Yangtze middle-lower",
  SC = "Southern coast",
  QTP = "Qinghai-Tibet Plateau"
)

primary_seed <- 20260726L
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
  # dd6310c added phylogeny export; re-exported count/taxonomy/metadata hashes
  # were checked unchanged before updating this preparation-script pin.
  prepare_script = "e7dc61c63e7118b719e69605bd1d2e9c199c1d77e6f18bde34c5d9dee5e9794e",
  theme_pub = "8d3821a485aeb529b12184bbe977c1ca952131a22d5eb69b802088ce91909beb"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 26 input(s): ",
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
  ragg = "1.3.2",
  readr = "2.1.5",
  scales = "1.3.0",
  svglite = "2.1.3",
  systemfonts = "1.1.0"
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
    "Community composition relative abundance prevalence phylum genus",
    "inland coastal Tibetan wetland unknown other 0123456789"
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
storage.mode(counts) <- "numeric"
feature_ids <- rownames(counts)
sample_ids <- colnames(counts)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]

metadata$SampleID <- rownames(metadata)
metadata$GroupDisplay <- unname(group_display[metadata$Group])
metadata$TypeDisplay <- unname(type_display[metadata$Type])
metadata$SampleNumber <- as.integer(sub("^S", "", metadata$SampleID))
metadata$GroupDisplay <- factor(
  metadata$GroupDisplay,
  levels = unname(group_display[c("IW", "CW", "TW")])
)
sample_order <- metadata$SampleID[order(
  metadata$GroupDisplay,
  metadata$TypeDisplay,
  metadata$SampleNumber
)]
metadata$SampleOrder <- match(metadata$SampleID, sample_order)

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_depths <- colSums(counts)
relative_feature <- sweep(counts, 2L, sample_depths, "/")

add_check("feature-count", "data", nrow(counts), 13628, nrow(counts) == 13628L)
add_check("sample-count", "data", ncol(counts), 90, ncol(counts) == 90L)
add_check("read-count", "data", sum(counts), 1619670, sum(counts) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628L)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90L)
add_check("taxonomy-seven-ranks", "data", ncol(taxonomy), 7, identical(colnames(taxonomy), rank_names))
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
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0))
add_check("counts-integer", "data", max(abs(counts - round(counts))), 0, all(counts == round(counts)))
add_check("nonempty-features", "data", sum(rowSums(counts) > 0), nrow(counts), all(rowSums(counts) > 0))
add_check("nonempty-samples", "data", sum(sample_depths > 0), ncol(counts), all(sample_depths > 0))
add_check("group-levels", "design", paste(sort(unique(metadata$Group)), collapse = ","), "CW,IW,TW", identical(sort(unique(metadata$Group)), c("CW", "IW", "TW")))
add_check("group-balanced", "design", paste(as.integer(table(metadata$Group)), collapse = ","), "30,30,30", all(table(metadata$Group) == 30L))
add_check("group-labels-english", "design", sum(grepl("[一-鿿]", metadata$GroupDisplay)), 0, !any(grepl("[一-鿿]", metadata$GroupDisplay)))
add_check("source-summary-samples", "provenance", source_summary$samples, 90, identical(source_summary$samples, 90L))
add_check("source-summary-features", "provenance", source_summary$features, 13628, identical(source_summary$features, 13628L))
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
  TypeCode = metadata$Type,
  Type = metadata$TypeDisplay,
  LibrarySize = as.numeric(sample_depths[metadata$SampleID]),
  DisplayOrder = metadata$SampleOrder,
  stringsAsFactors = FALSE
)

rank_matrices <- setNames(vector("list", length(rank_names)), rank_names)
rank_audit_rows <- list()
for (rank_name in rank_names) {
  rank_counts <- aggregate_rank(counts, taxonomy, rank_name)
  rank_rel <- sweep(rank_counts, 2L, sample_depths, "/")
  rank_matrices[[rank_name]] <- rank_rel
  unknown_label <- paste("Unknown", tolower(rank_name))
  known <- rownames(rank_rel) != unknown_label
  known_order <- rownames(rank_rel)[known][order(
    rowMeans(rank_rel[known, , drop = FALSE]),
    decreasing = TRUE
  )]
  top10 <- head(known_order, 10L)
  unknown_mean <- if (unknown_label %in% rownames(rank_rel)) {
    mean(rank_rel[unknown_label, ])
  } else {
    0
  }
  unknown_pooled <- if (unknown_label %in% rownames(rank_counts)) {
    sum(rank_counts[unknown_label, ]) / sum(rank_counts)
  } else {
    0
  }
  feature_labels <- clean_taxon(taxonomy[[rank_name]], rank_name)
  rank_audit_rows[[length(rank_audit_rows) + 1L]] <- data.frame(
    Rank = rank_name,
    Labels = nrow(rank_rel),
    KnownLabels = sum(known),
    UnknownFeatures = sum(feature_labels == unknown_label),
    UnknownFeatureFraction = mean(feature_labels == unknown_label),
    UnknownReadFractionPooled = unknown_pooled,
    UnknownEqualSampleMean = unknown_mean,
    Top10KnownCoverageMean = mean(colSums(rank_rel[top10, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-closure"), "taxonomy",
    max(abs(colSums(rank_rel) - 1)), "<=1e-12",
    max(abs(colSums(rank_rel) - 1)) <= 1e-12
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-labels-unique"), "taxonomy",
    anyDuplicated(rownames(rank_rel)), 0, !anyDuplicated(rownames(rank_rel))
  )
}
taxonomy_rank_audit <- do.call(rbind, rank_audit_rows)

phylum_rel <- rank_matrices[["Phylum"]]
genus_rel <- rank_matrices[["Genus"]]
unknown_phylum <- "Unknown phylum"
unknown_genus <- "Unknown genus"

known_phyla <- setdiff(rownames(phylum_rel), unknown_phylum)
equal_phylum_mean <- rowMeans(phylum_rel)
pooled_phylum_share <- rowSums(aggregate_rank(counts, taxonomy, "Phylum")) / sum(counts)
top10_phyla <- known_phyla[order(equal_phylum_mean[known_phyla], decreasing = TRUE)][1:10]
pooled_top10 <- known_phyla[order(pooled_phylum_share[known_phyla], decreasing = TRUE)][1:10]
top10_jaccard <- length(intersect(top10_phyla, pooled_top10)) /
  length(union(top10_phyla, pooled_top10))

phylum_display_labels <- c(top10_phyla, "Other known phyla", unknown_phylum)
phylum_assignment <- ifelse(
  rownames(phylum_rel) == unknown_phylum,
  unknown_phylum,
  ifelse(rownames(phylum_rel) %in% top10_phyla, rownames(phylum_rel), "Other known phyla")
)
phylum_display <- rowsum(
  phylum_rel,
  group = factor(phylum_assignment, levels = phylum_display_labels),
  reorder = FALSE
)
phylum_display <- phylum_display[phylum_display_labels, , drop = FALSE]

phylum_relative_abundance <- matrix_to_long(phylum_rel, "RelativeAbundance")
phylum_relative_abundance$Percent <- 100 * phylum_relative_abundance$RelativeAbundance
phylum_relative_abundance$Annotated <- phylum_relative_abundance$Taxon != unknown_phylum

phylum_display_composition <- matrix_to_long(phylum_display, "RelativeAbundance")
phylum_display_composition$Percent <- 100 * phylum_display_composition$RelativeAbundance
phylum_display_composition$GroupCode <- metadata[phylum_display_composition$SampleID, "Group"]
phylum_display_composition$Group <- as.character(
  metadata[phylum_display_composition$SampleID, "GroupDisplay"]
)
phylum_display_composition$Type <- metadata[phylum_display_composition$SampleID, "TypeDisplay"]
phylum_display_composition$SampleOrder <- metadata[
  phylum_display_composition$SampleID, "SampleOrder"
]
phylum_display_composition$DisplayClass <- ifelse(
  phylum_display_composition$Taxon == unknown_phylum,
  "Unknown",
  ifelse(phylum_display_composition$Taxon == "Other known phyla", "Other known", "Top 10 known")
)

group_codes <- c("IW", "CW", "TW")
phylum_group_rows <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  means <- rowMeans(phylum_display[, ids, drop = FALSE])
  pooled_counts <- aggregate_rank(counts[, ids, drop = FALSE], taxonomy, "Phylum")
  pooled_assignment <- ifelse(
    rownames(pooled_counts) == unknown_phylum,
    unknown_phylum,
    ifelse(rownames(pooled_counts) %in% top10_phyla, rownames(pooled_counts), "Other known phyla")
  )
  pooled_display <- rowsum(
    pooled_counts,
    group = factor(pooled_assignment, levels = phylum_display_labels),
    reorder = FALSE
  )
  pooled_share <- rowSums(pooled_display[phylum_display_labels, , drop = FALSE]) /
    sum(pooled_display)
  phylum_group_rows[[length(phylum_group_rows) + 1L]] <- data.frame(
    GroupCode = group_code,
    Group = unname(group_display[group_code]),
    Taxon = phylum_display_labels,
    EqualSampleMean = as.numeric(means[phylum_display_labels]),
    PooledReadShare = pooled_share,
    AbsoluteDifferencePP = 100 * abs(as.numeric(means[phylum_display_labels]) - pooled_share),
    stringsAsFactors = FALSE
  )
}
phylum_group_summary <- do.call(rbind, phylum_group_rows)

top_n_values <- c(5L, 8L, 10L, 15L, 20L)
top_n_rows <- lapply(top_n_values, function(n) {
  selected <- head(
    known_phyla[order(equal_phylum_mean[known_phyla], decreasing = TRUE)],
    n
  )
  coverage <- colSums(phylum_rel[selected, , drop = FALSE])
  data.frame(
    TopN = n,
    SelectedTaxa = paste(selected, collapse = "; "),
    MeanCoverage = mean(coverage),
    MinimumCoverage = min(coverage),
    MaximumCoverage = max(coverage),
    OtherKnownMean = mean(1 - coverage - phylum_rel[unknown_phylum, ]),
    UnknownMean = mean(phylum_rel[unknown_phylum, ]),
    stringsAsFactors = FALSE
  )
})
top_n_sensitivity <- do.call(rbind, top_n_rows)

ranking_taxa <- union(top10_phyla, pooled_top10)
ranking_sensitivity <- data.frame(
  Taxon = ranking_taxa,
  EqualSampleRank = match(ranking_taxa, known_phyla[
    order(equal_phylum_mean[known_phyla], decreasing = TRUE)
  ]),
  PooledReadRank = match(ranking_taxa, known_phyla[
    order(pooled_phylum_share[known_phyla], decreasing = TRUE)
  ]),
  EqualSampleMean = as.numeric(equal_phylum_mean[ranking_taxa]),
  PooledReadShare = as.numeric(pooled_phylum_share[ranking_taxa]),
  DifferencePP = 100 * (
    as.numeric(equal_phylum_mean[ranking_taxa]) -
      as.numeric(pooled_phylum_share[ranking_taxa])
  ),
  InEqualSampleTop10 = ranking_taxa %in% top10_phyla,
  InPooledReadTop10 = ranking_taxa %in% pooled_top10,
  stringsAsFactors = FALSE
)

known_genera <- setdiff(rownames(genus_rel), unknown_genus)
equal_genus_mean <- rowMeans(genus_rel)
top20_genera <- known_genera[order(equal_genus_mean[known_genera], decreasing = TRUE)][1:20]
top15_genera <- top20_genera[1:15]

genus_relative_abundance <- matrix_to_long(genus_rel, "RelativeAbundance")
genus_relative_abundance$Percent <- 100 * genus_relative_abundance$RelativeAbundance
genus_relative_abundance$Annotated <- genus_relative_abundance$Taxon != unknown_genus

genus_bubble_rows <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  for (taxon in c(top15_genera, unknown_genus)) {
    values <- genus_rel[taxon, ids]
    genus_bubble_rows[[length(genus_bubble_rows) + 1L]] <- data.frame(
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      Taxon = taxon,
      MeanRelativeAbundance = mean(values),
      MeanPercent = 100 * mean(values),
      Prevalence = mean(values > 0),
      stringsAsFactors = FALSE
    )
  }
}
genus_bubble_summary <- do.call(rbind, genus_bubble_rows)

denominator_rows <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  known_total <- 1 - genus_rel[unknown_genus, ids]
  for (taxon in top15_genera) {
    all_read <- mean(genus_rel[taxon, ids])
    known_only <- mean(genus_rel[taxon, ids] / known_total)
    denominator_rows[[length(denominator_rows) + 1L]] <- data.frame(
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      Taxon = taxon,
      AllReadMean = all_read,
      KnownOnlyMean = known_only,
      InflationPP = 100 * (known_only - all_read),
      UnknownMean = mean(genus_rel[unknown_genus, ids]),
      stringsAsFactors = FALSE
    )
  }
}
denominator_sensitivity <- do.call(rbind, denominator_rows)

aggregation_group_rows <- list()
for (group_code in group_codes) {
  ids <- metadata$SampleID[metadata$Group == group_code]
  equal_all <- rowMeans(phylum_rel[, ids, drop = FALSE])
  group_counts <- aggregate_rank(counts[, ids, drop = FALSE], taxonomy, "Phylum")
  pooled_all <- rowSums(group_counts) / sum(group_counts)
  taxa <- union(names(equal_all), names(pooled_all))
  equal_vector <- setNames(rep(0, length(taxa)), taxa)
  pooled_vector <- equal_vector
  equal_vector[names(equal_all)] <- equal_all
  pooled_vector[names(pooled_all)] <- pooled_all
  aggregation_group_rows[[length(aggregation_group_rows) + 1L]] <- data.frame(
    GroupCode = group_code,
    Group = unname(group_display[group_code]),
    TotalVariation = 0.5 * sum(abs(equal_vector - pooled_vector)),
    MaximumAbsoluteDifferencePP = 100 * max(abs(equal_vector - pooled_vector)),
    EqualSampleClosure = sum(equal_vector),
    PooledReadClosure = sum(pooled_vector),
    stringsAsFactors = FALSE
  )
}
aggregation_sensitivity <- do.call(rbind, aggregation_group_rows)

closure_audit <- data.frame(
  SampleID = sample_ids,
  FeatureClosure = colSums(relative_feature),
  PhylumClosure = colSums(phylum_rel),
  PhylumDisplayClosure = colSums(phylum_display),
  GenusClosure = colSums(genus_rel),
  MaximumAbsoluteError = pmax(
    abs(colSums(relative_feature) - 1),
    abs(colSums(phylum_rel) - 1),
    abs(colSums(phylum_display) - 1),
    abs(colSums(genus_rel) - 1)
  ),
  stringsAsFactors = FALSE
)

add_check("phylum-label-count", "composition", nrow(phylum_rel), 56, nrow(phylum_rel) == 56L)
add_check("genus-label-count", "composition", nrow(genus_rel), 772, nrow(genus_rel) == 772L)
add_check("phylum-top-n", "composition", length(top10_phyla), 10, length(top10_phyla) == 10L)
add_check("genus-bubble-known-top-n", "composition", length(top15_genera), 15, length(top15_genera) == 15L)
add_check("genus-heatmap-top-n", "composition", length(top20_genera), 20, length(top20_genera) == 20L)
add_check("unknown-phylum-explicit", "composition", unknown_phylum %in% rownames(phylum_display), TRUE, unknown_phylum %in% rownames(phylum_display))
add_check("other-known-explicit", "composition", "Other known phyla" %in% rownames(phylum_display), TRUE, "Other known phyla" %in% rownames(phylum_display))
add_check("display-closure", "composition", max(abs(colSums(phylum_display) - 1)), "<=1e-12", max(abs(colSums(phylum_display) - 1)) <= 1e-12)
add_check("top10-ranking-jaccard", "sensitivity", top10_jaccard, 1, near(top10_jaccard, 1))
add_check("top-n-branches", "sensitivity", nrow(top_n_sensitivity), 5, nrow(top_n_sensitivity) == 5L)
add_check("top10-coverage-mean", "composition", mean(colSums(phylum_rel[top10_phyla, , drop = FALSE])), 0.94803047094, near(mean(colSums(phylum_rel[top10_phyla, , drop = FALSE])), 0.94803047094, 1e-10))
add_check("unknown-phylum-mean", "composition", mean(phylum_rel[unknown_phylum, ]), 0.00346378239224, near(mean(phylum_rel[unknown_phylum, ]), 0.00346378239224, 1e-12))
add_check("unknown-genus-mean", "composition", mean(genus_rel[unknown_genus, ]), 0.63105611307, near(mean(genus_rel[unknown_genus, ]), 0.63105611307, 1e-10))
add_check("unknown-species-mean", "composition", taxonomy_rank_audit$UnknownEqualSampleMean[taxonomy_rank_audit$Rank == "Species"], 0.996345023736, near(taxonomy_rank_audit$UnknownEqualSampleMean[taxonomy_rank_audit$Rank == "Species"], 0.996345023736, 1e-10))
add_check("denominator-inflation-positive", "sensitivity", max(denominator_sensitivity$InflationPP), ">0", max(denominator_sensitivity$InflationPP) > 0)
add_check("aggregation-branches", "sensitivity", nrow(aggregation_sensitivity), 3, nrow(aggregation_sensitivity) == 3L)
add_check("all-closure-checks", "composition", max(closure_audit$MaximumAbsoluteError), "<=1e-12", max(closure_audit$MaximumAbsoluteError) <= 1e-12)
add_check("all-read-denominator-retained", "composition", all(genus_bubble_summary$MeanRelativeAbundance <= 1), TRUE, all(genus_bubble_summary$MeanRelativeAbundance <= 1))

taxon_palette <- c(
  Proteobacteria = "#0072B2",
  Chloroflexi = "#E69F00",
  Bacteroidetes = "#009E73",
  Acidobacteria = "#D55E00",
  Actinobacteria = "#CC79A7",
  Firmicutes = "#56B4E9",
  Verrucomicrobia = "#F0E442",
  Planctomycetes = "#8C564B",
  Gemmatimonadetes = "#17BECF",
  Nitrospirae = "#9467BD",
  `Other known phyla` = "#B8B8B8",
  `Unknown phylum` = "#2F2F2F"
)
taxon_palette <- taxon_palette[phylum_display_labels]

stack_data <- phylum_display_composition
stack_data$SampleID <- factor(stack_data$SampleID, levels = sample_order)
stack_data$Taxon <- factor(stack_data$Taxon, levels = rev(phylum_display_labels))
stack_data$Group <- factor(stack_data$Group, levels = unname(group_display[group_codes]))
stack_plot <- ggplot2::ggplot(
  stack_data,
  ggplot2::aes(x = SampleID, y = RelativeAbundance, fill = Taxon)
) +
  ggplot2::geom_col(width = 1, linewidth = 0) +
  ggplot2::facet_grid(~Group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_manual(values = rev(taxon_palette), drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Phylum-level composition across 90 wetland samples",
    subtitle = "Top 10 known phyla by equal-sample mean; unknown reads remain explicit",
    x = "Samples ordered by wetland group and region",
    y = "Relative abundance",
    fill = "Phylum",
    caption = "Each bar closes to 100%; Other contains known phyla only."
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.5, "mm"),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(3.3, "mm")
  )

bubble_data <- genus_bubble_summary
bubble_data$Group <- factor(bubble_data$Group, levels = unname(group_display[group_codes]))
bubble_order <- c(top15_genera, unknown_genus)
bubble_data$Taxon <- factor(bubble_data$Taxon, levels = rev(bubble_order))
bubble_plot <- ggplot2::ggplot(
  bubble_data,
  ggplot2::aes(
    x = Group,
    y = Taxon,
    size = MeanPercent,
    colour = Prevalence
  )
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_size_area(max_size = 13, breaks = c(1, 5, 20, 60)) +
  ggplot2::scale_colour_gradientn(
    colours = c("#D9EAF7", "#56B4E9", "#0072B2", "#003B5C"),
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      "Inland wetland" = "Inland\nwetland",
      "Coastal wetland" = "Coastal\nwetland",
      "Tibetan wetland" = "Tibetan\nwetland"
    )
  ) +
  ggplot2::labs(
    title = "Genus abundance and prevalence answer different questions",
    subtitle = "Top 15 known genera plus the explicit unclassified fraction",
    x = NULL,
    y = NULL,
    size = "Mean abundance (%)",
    colour = "Prevalence",
    caption = paste(
      "Means retain all reads in the denominator; prevalence is the fraction",
      "of non-zero samples."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )

make_flow_polygons <- function(summary_table, taxa, groups) {
  gap <- 0.12
  source_base <- setNames((seq_along(groups) - 1) * (1 + gap), groups)
  source_cursor <- source_base
  target_totals <- vapply(
    taxa,
    function(taxon) sum(summary_table$EqualSampleMean[summary_table$Taxon == taxon]),
    numeric(1)
  )
  target_gap <- 0.035
  target_base <- cumsum(c(0, head(target_totals, -1) + target_gap))
  names(target_base) <- taxa
  target_cursor <- target_base
  polygon_rows <- list()
  flow_id <- 0L
  for (group in groups) {
    for (taxon in taxa) {
      value <- summary_table$EqualSampleMean[
        summary_table$Group == group & summary_table$Taxon == taxon
      ]
      if (length(value) != 1L || value <= 0) next
      source_y0 <- source_cursor[[group]]
      source_y1 <- source_y0 + value
      target_y0 <- target_cursor[[taxon]]
      target_y1 <- target_y0 + value
      source_cursor[[group]] <- source_y1
      target_cursor[[taxon]] <- target_y1
      t <- seq(0, 1, length.out = 30)
      smooth <- 3 * t^2 - 2 * t^3
      x_forward <- 0.12 + 0.76 * t
      lower <- source_y0 + (target_y0 - source_y0) * smooth
      upper <- source_y1 + (target_y1 - source_y1) * smooth
      flow_id <- flow_id + 1L
      polygon_rows[[length(polygon_rows) + 1L]] <- data.frame(
        FlowID = flow_id,
        Group = group,
        Taxon = taxon,
        x = c(x_forward, rev(x_forward)),
        y = c(lower, rev(upper)),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    polygons = do.call(rbind, polygon_rows),
    source = data.frame(
      Group = groups,
      ymin = unname(source_base),
      ymax = unname(source_base) + 1,
      stringsAsFactors = FALSE
    ),
    target = data.frame(
      Taxon = taxa,
      ymin = unname(target_base),
      ymax = unname(target_base + target_totals),
      stringsAsFactors = FALSE
    )
  )
}

flow_parts <- make_flow_polygons(
  phylum_group_summary,
  phylum_display_labels,
  unname(group_display[group_codes])
)
alluvial_plot <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = flow_parts$polygons,
    ggplot2::aes(x = x, y = y, group = FlowID, fill = Taxon),
    alpha = 0.78,
    colour = NA
  ) +
  ggplot2::geom_rect(
    data = flow_parts$source,
    ggplot2::aes(xmin = 0.07, xmax = 0.12, ymin = ymin, ymax = ymax),
    fill = "#4B5563",
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_rect(
    data = flow_parts$target,
    ggplot2::aes(xmin = 0.88, xmax = 0.93, ymin = ymin, ymax = ymax, fill = Taxon),
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = flow_parts$source,
    ggplot2::aes(x = 0.055, y = (ymin + ymax) / 2, label = Group),
    hjust = 1,
    size = 2.7,
    family = font_family
  ) +
  ggplot2::geom_text(
    data = flow_parts$target,
    ggplot2::aes(x = 0.945, y = (ymin + ymax) / 2, label = Taxon),
    hjust = 0,
    size = 2.55,
    family = font_family
  ) +
  ggplot2::scale_fill_manual(values = taxon_palette, drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(-0.12, 1.23), clip = "off") +
  ggplot2::labs(
    title = "Group means reveal how wetland composition is redistributed",
    subtitle = "Each group contributes one equal-sample composition; ribbon width is mean relative abundance",
    caption = "The flow diagram is descriptive: it does not represent migration, succession, or causality."
  ) +
  ggplot2::theme_void(base_family = font_family, base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 10.5),
    plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
    plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
    legend.position = "none",
    plot.margin = ggplot2::margin(5.5, 65, 5.5, 78)
  )

heatmap_matrix <- log1p(10000 * genus_rel[top20_genera, sample_order, drop = FALSE])
row_order <- rownames(heatmap_matrix)[
  hclust(dist(heatmap_matrix), method = "complete")$order
]
genus_heatmap_matrix <- data.frame(
  Genus = rownames(heatmap_matrix),
  heatmap_matrix,
  check.names = FALSE
)
heatmap_data <- matrix_to_long(heatmap_matrix, "Log1pAbundance")
heatmap_data$SampleID <- factor(heatmap_data$SampleID, levels = sample_order)
heatmap_data$Taxon <- factor(heatmap_data$Taxon, levels = rev(row_order))
heatmap_data$Group <- factor(
  as.character(metadata[as.character(heatmap_data$SampleID), "GroupDisplay"]),
  levels = unname(group_display[group_codes])
)
heatmap_plot <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = SampleID, y = Taxon, fill = Log1pAbundance)
) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(~Group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
  ) +
  ggplot2::labs(
    title = "Sample-level heatmap preserves within-group heterogeneity",
    subtitle = "Top 20 known genera; rows clustered, columns ordered by metadata rather than outcome",
    x = "Samples ordered by wetland group and region",
    y = NULL,
    fill = "log1p(relative\nabundance x 10,000)",
    caption = "Display transform only; the exported abundance table remains on the original relative scale."
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.5, "mm"),
    legend.position = "right"
  )

figure_specification <- data.frame(
  stem = c(
    "26-phylum-stacked",
    "26-genus-bubble",
    "26-group-phylum-alluvial",
    "26-genus-heatmap"
  ),
  width_mm = c(183, 140, 183, 183),
  height_mm = c(112, 118, 112, 122),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(stack_plot, bubble_plot, alluvial_plot, heatmap_plot)
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
write_tsv(taxonomy_rank_audit, "taxonomy-rank-audit.tsv")
write_tsv(phylum_relative_abundance, "phylum-relative-abundance.tsv")
write_tsv(phylum_display_composition, "phylum-display-composition.tsv")
write_tsv(phylum_group_summary, "phylum-group-summary.tsv")
write_tsv(genus_relative_abundance, "genus-relative-abundance.tsv")
write_tsv(genus_bubble_summary, "genus-bubble-summary.tsv")
write_tsv(genus_heatmap_matrix, "genus-heatmap-matrix.tsv")
write_tsv(data.frame(RowOrder = seq_along(row_order), Genus = row_order), "genus-heatmap-row-order.tsv")
write_tsv(top_n_sensitivity, "top-n-sensitivity.tsv")
write_tsv(ranking_sensitivity, "ranking-sensitivity.tsv")
write_tsv(denominator_sensitivity, "denominator-sensitivity.tsv")
write_tsv(aggregation_sensitivity, "aggregation-sensitivity.tsv")
write_tsv(closure_audit, "closure-audit.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

top10_coverage <- colSums(phylum_rel[top10_phyla, , drop = FALSE])
summary_payload <- list(
  article = 26L,
  title = "Community composition: stacked bars, bubbles, flows, and heatmaps",
  dataset = list(
    source = "microeco 2.0.0 packaged Chinese wetland 16S tables",
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    groups = nlevels(metadata$GroupDisplay),
    library_depth_min = min(sample_depths),
    library_depth_max = max(sample_depths)
  ),
  taxonomy = list(
    ranks = length(rank_names),
    phylum_labels = nrow(phylum_rel),
    genus_labels = nrow(genus_rel),
    phylum_unknown_equal_sample_mean = mean(phylum_rel[unknown_phylum, ]),
    genus_unknown_equal_sample_mean = mean(genus_rel[unknown_genus, ]),
    species_unknown_equal_sample_mean = taxonomy_rank_audit$UnknownEqualSampleMean[
      taxonomy_rank_audit$Rank == "Species"
    ]
  ),
  primary = list(
    rank = "Phylum",
    top_n = length(top10_phyla),
    top_taxa = top10_phyla,
    top10_coverage_mean = mean(top10_coverage),
    top10_coverage_min = min(top10_coverage),
    top10_coverage_max = max(top10_coverage),
    unknown_separate_from_other = TRUE,
    aggregation = "Equal-sample mean after per-sample closure",
    inferential_claim = FALSE
  ),
  sensitivity = list(
    top_n_branches = nrow(top_n_sensitivity),
    equal_vs_pooled_top10_jaccard = top10_jaccard,
    maximum_known_only_inflation_pp = max(denominator_sensitivity$InflationPP),
    maximum_equal_vs_pooled_group_taxon_difference_pp = max(
      phylum_group_summary$AbsoluteDifferencePP
    ),
    maximum_group_total_variation = max(aggregation_sensitivity$TotalVariation)
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
  file.path(output_dir, "community-composition-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

validation_lines <- c(
  "Article 26 community-composition validation",
  paste0("Status: ", ifelse(checks_failed == 0L, "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total),
  paste0("Samples/features/reads: 90/13628/1619670"),
  paste0("Phylum Top 10 mean coverage: ", sprintf("%.4f", mean(top10_coverage))),
  paste0("Unknown phylum mean: ", sprintf("%.4f", mean(phylum_rel[unknown_phylum, ]))),
  paste0("Unknown genus mean: ", sprintf("%.4f", mean(genus_rel[unknown_genus, ]))),
  paste0("Equal-sample vs pooled Top 10 Jaccard: ", sprintf("%.3f", top10_jaccard)),
  paste0("Maximum known-only inflation: ", sprintf("%.3f pp", max(denominator_sensitivity$InflationPP))),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("R version: ", r_version)
)
writeLines(validation_lines, file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed_ids <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 26 validation failed: ",
    paste(failed_ids, collapse = ", "),
    call. = FALSE
  )
}

cat(
  sprintf(
    "Article 26 validation passed: %d/%d checks.\n",
    checks_passed,
    checks_total
  )
)
