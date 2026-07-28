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

matrix_to_long <- function(x, rank_name) {
  data.frame(
    Rank = rank_name,
    Taxon = rep(rownames(x), times = ncol(x)),
    RankQualifiedTaxon = paste0(
      substr(rank_name, 1L, 1L), " · ",
      rep(rownames(x), times = ncol(x))
    ),
    SampleID = rep(colnames(x), each = nrow(x)),
    RelativeAbundance = as.vector(x),
    Percent = 100 * as.vector(x),
    stringsAsFactors = FALSE
  )
}

wrap_label <- function(x, width = 16L) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

group_display <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan wetland"
)
group_codes <- c("IW", "CW", "TW")
rank_names <- c("Phylum", "Class", "Order", "Family", "Genus")
rank_prefix <- c(Phylum = "P", Class = "C", Order = "O", Family = "F", Genus = "G")
primary_seed <- 20260727L
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
    "Missing Article 27 input(s): ",
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
    "Phylum class order family genus lineage known unknown",
    "inland coastal Tibetan wetland prevalence abundance 0123456789"
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
metadata$GroupDisplay <- factor(
  metadata$GroupDisplay,
  levels = unname(group_display[group_codes])
)

sample_depths <- colSums(counts)
relative_feature <- sweep(counts, 2L, sample_depths, "/")
feature_equal_mean <- rowMeans(relative_feature)
all_rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

add_check("feature-count", "data", nrow(counts), 13628, nrow(counts) == 13628L)
add_check("sample-count", "data", ncol(counts), 90, ncol(counts) == 90L)
add_check("read-count", "data", sum(counts), 1619670, sum(counts) == 1619670)
add_check("taxonomy-feature-count", "data", nrow(taxonomy), 13628, nrow(taxonomy) == 13628L)
add_check("metadata-sample-count", "data", nrow(metadata), 90, nrow(metadata) == 90L)
add_check(
  "taxonomy-seven-ranks", "data", paste(colnames(taxonomy), collapse = ","),
  paste(all_rank_names, collapse = ","), identical(colnames(taxonomy), all_rank_names)
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
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0))
add_check("counts-integer", "data", max(abs(counts - round(counts))), 0, all(counts == round(counts)))
add_check("nonempty-features", "data", sum(rowSums(counts) > 0), nrow(counts), all(rowSums(counts) > 0))
add_check("nonempty-samples", "data", sum(sample_depths > 0), ncol(counts), all(sample_depths > 0))
add_check(
  "group-levels", "design", paste(sort(unique(metadata$Group)), collapse = ","),
  "CW,IW,TW", identical(sort(unique(metadata$Group)), c("CW", "IW", "TW"))
)
add_check(
  "group-balanced", "design", paste(as.integer(table(metadata$Group)), collapse = ","),
  "30,30,30", all(table(metadata$Group) == 30L)
)
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
unknown_labels <- setNames(paste("Unknown", tolower(rank_names)), rank_names)
known_matrix <- vapply(
  rank_names,
  function(rank_name) taxonomy_clean[[rank_name]] != unknown_labels[[rank_name]],
  logical(nrow(taxonomy_clean))
)
colnames(known_matrix) <- rank_names

complete_depth <- apply(known_matrix, 1L, function(x) {
  first_unknown <- which(!x)[1L]
  if (is.na(first_unknown)) length(x) else first_unknown - 1L
})
lineage_key <- apply(taxonomy_clean[, rank_names, drop = FALSE], 1L, paste, collapse = " > ")
taxonomy_lineage_map <- data.frame(
  FeatureID = feature_ids,
  taxonomy_clean[, rank_names, drop = FALSE],
  FullLineageKey = lineage_key,
  CompleteDepth = complete_depth,
  DeepestContiguousRank = c("Below phylum", rank_names)[complete_depth + 1L],
  EqualSampleMean = feature_equal_mean,
  PooledReadFraction = rowSums(counts) / sum(counts),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rank_relative <- list()
rank_counts <- list()
rank_resolution_rows <- list()
multirank_long_rows <- list()
expected_label_counts <- c(Phylum = 56L, Class = 128L, Order = 217L, Family = 394L, Genus = 772L)
expected_known_means <- c(
  Phylum = 0.996536217607761,
  Class = 0.967978544271532,
  Order = 0.855727280469946,
  Family = 0.679627539342075,
  Genus = 0.368943886930081
)

for (rank_name in rank_names) {
  labels <- taxonomy_clean[[rank_name]]
  counts_rank <- rowsum(counts, group = labels, reorder = TRUE)
  rel_rank <- sweep(counts_rank, 2L, sample_depths, "/")
  rank_counts[[rank_name]] <- counts_rank
  rank_relative[[rank_name]] <- rel_rank
  multirank_long_rows[[rank_name]] <- matrix_to_long(rel_rank, rank_name)
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_order <- known_taxa[order(rowMeans(rel_rank[known_taxa, , drop = FALSE]), decreasing = TRUE)]
  top5 <- head(known_order, 5L)
  top10 <- head(known_order, 10L)
  known_sample <- 1 - rel_rank[unknown_label, ]
  incertae_features <- grepl("Incertae Sedis", labels, fixed = TRUE)
  rank_resolution_rows[[rank_name]] <- data.frame(
    Rank = rank_name,
    RankIndex = match(rank_name, rank_names),
    Labels = nrow(rel_rank),
    KnownLabels = length(known_taxa),
    UnknownFeatures = sum(labels == unknown_label),
    UnknownFeatureFraction = mean(labels == unknown_label),
    UnknownReadFractionPooled = sum(counts_rank[unknown_label, ]) / sum(counts_rank),
    UnknownEqualSampleMean = mean(rel_rank[unknown_label, ]),
    KnownEqualSampleMean = mean(known_sample),
    KnownMinimum = min(known_sample),
    KnownMaximum = max(known_sample),
    Top5KnownCoverageMean = mean(colSums(rel_rank[top5, , drop = FALSE])),
    Top10KnownCoverageMean = mean(colSums(rel_rank[top10, , drop = FALSE])),
    IncertaeSedisFeatures = sum(incertae_features),
    IncertaeSedisEqualSampleMean = sum(feature_equal_mean[incertae_features]),
    stringsAsFactors = FALSE
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-closure"), "taxonomy",
    max(abs(colSums(rel_rank) - 1)), "<=1e-12",
    max(abs(colSums(rel_rank) - 1)) <= 1e-12
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-labels-unique"), "taxonomy",
    anyDuplicated(rownames(rel_rank)), 0, !anyDuplicated(rownames(rel_rank))
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-label-count"), "taxonomy",
    nrow(rel_rank), expected_label_counts[[rank_name]],
    nrow(rel_rank) == expected_label_counts[[rank_name]]
  )
  add_check(
    paste0("rank-", tolower(rank_name), "-known-mean"), "taxonomy",
    mean(known_sample), expected_known_means[[rank_name]],
    near(mean(known_sample), expected_known_means[[rank_name]], 1e-12)
  )
}

rank_resolution_audit <- do.call(rbind, rank_resolution_rows)
rank_resolution_audit$ResolutionLossFromParentPP <- c(
  NA_real_,
  100 * head(rank_resolution_audit$KnownEqualSampleMean, -1L) -
    100 * tail(rank_resolution_audit$KnownEqualSampleMean, -1L)
)
multirank_relative_abundance <- do.call(rbind, multirank_long_rows)
multirank_relative_abundance$GroupCode <- metadata[
  multirank_relative_abundance$SampleID, "Group"
]
multirank_relative_abundance$Group <- as.character(metadata[
  multirank_relative_abundance$SampleID, "GroupDisplay"
])

group_resolution_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_sample <- 1 - rel_rank[unknown_label, ]
  group_resolution_rows[[length(group_resolution_rows) + 1L]] <- data.frame(
    GroupCode = "ALL",
    Group = "All samples",
    Rank = rank_name,
    RankIndex = match(rank_name, rank_names),
    Samples = length(sample_ids),
    KnownMean = mean(known_sample),
    KnownMinimum = min(known_sample),
    KnownMaximum = max(known_sample),
    UnknownMean = mean(1 - known_sample),
    stringsAsFactors = FALSE
  )
  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    values <- known_sample[ids]
    group_resolution_rows[[length(group_resolution_rows) + 1L]] <- data.frame(
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      Rank = rank_name,
      RankIndex = match(rank_name, rank_names),
      Samples = length(ids),
      KnownMean = mean(values),
      KnownMinimum = min(values),
      KnownMaximum = max(values),
      UnknownMean = mean(1 - values),
      stringsAsFactors = FALSE
    )
  }
}
group_rank_resolution <- do.call(rbind, group_resolution_rows)

lineage_gap_rows <- list()
for (i in 2:length(rank_names)) {
  parent_rank <- rank_names[[i - 1L]]
  child_rank <- rank_names[[i]]
  parent_unknown <- taxonomy_clean[[parent_rank]] == unknown_labels[[parent_rank]]
  child_unknown <- taxonomy_clean[[child_rank]] == unknown_labels[[child_rank]]
  gap_mask <- parent_unknown & !child_unknown
  loss_mask <- !parent_unknown & child_unknown
  lineage_gap_rows[[length(lineage_gap_rows) + 1L]] <- data.frame(
    ParentRank = parent_rank,
    ChildRank = child_rank,
    ChildKnownParentUnknownFeatures = sum(gap_mask),
    ChildKnownParentUnknownTaxa = length(unique(taxonomy_clean[[child_rank]][gap_mask])),
    ChildKnownParentUnknownEqualSampleMean = sum(feature_equal_mean[gap_mask]),
    ChildKnownParentUnknownPooledReadFraction = sum(counts[gap_mask, , drop = FALSE]) / sum(counts),
    ParentKnownChildUnknownFeatures = sum(loss_mask),
    ParentKnownChildUnknownEqualSampleMean = sum(feature_equal_mean[loss_mask]),
    stringsAsFactors = FALSE
  )
}
lineage_gap_audit <- do.call(rbind, lineage_gap_rows)

depth_labels <- c("Below phylum", rank_names)
complete_depth_rows <- lapply(0:length(rank_names), function(depth) {
  mask <- complete_depth == depth
  data.frame(
    CompleteDepth = depth,
    DeepestContiguousRank = depth_labels[[depth + 1L]],
    FeatureCount = sum(mask),
    FeatureFraction = mean(mask),
    PooledReadFraction = sum(counts[mask, , drop = FALSE]) / sum(counts),
    EqualSampleMean = sum(feature_equal_mean[mask]),
    stringsAsFactors = FALSE
  )
})
complete_depth_audit <- do.call(rbind, complete_depth_rows)

collision_rows <- list()
for (i in 2:length(rank_names)) {
  rank_name <- rank_names[[i]]
  unknown_label <- unknown_labels[[rank_name]]
  labels <- taxonomy_clean[[rank_name]]
  for (taxon in sort(setdiff(unique(labels), unknown_label))) {
    mask <- labels == taxon
    parent_paths <- apply(
      taxonomy_clean[mask, rank_names[seq_len(i - 1L)], drop = FALSE],
      1L,
      paste,
      collapse = " > "
    )
    parent_paths <- sort(unique(parent_paths))
    if (length(parent_paths) > 1L) {
      collision_rows[[length(collision_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        Taxon = taxon,
        ParentLineages = length(parent_paths),
        ParentPaths = paste(parent_paths, collapse = " | "),
        FeatureCount = sum(mask),
        EqualSampleMean = sum(feature_equal_mean[mask]),
        PooledReadFraction = sum(counts[mask, , drop = FALSE]) / sum(counts),
        stringsAsFactors = FALSE
      )
    }
  }
}
label_collision_audit <- if (length(collision_rows)) {
  do.call(rbind, collision_rows)
} else {
  data.frame(
    Rank = character(), Taxon = character(), ParentLineages = integer(),
    ParentPaths = character(), FeatureCount = integer(), EqualSampleMean = numeric(),
    PooledReadFraction = numeric(), stringsAsFactors = FALSE
  )
}

rank_label_pairs <- do.call(rbind, lapply(rank_names, function(rank_name) {
  labels <- setdiff(unique(taxonomy_clean[[rank_name]]), unknown_labels[[rank_name]])
  data.frame(Label = labels, Rank = rank_name, stringsAsFactors = FALSE)
}))
label_rank_split <- split(rank_label_pairs$Rank, rank_label_pairs$Label)
cross_rank_rows <- lapply(names(label_rank_split), function(label) {
  ranks <- unique(label_rank_split[[label]])
  if (length(ranks) <= 1L) return(NULL)
  data.frame(
    Label = label,
    RankCount = length(ranks),
    Ranks = paste(rank_names[rank_names %in% ranks], collapse = "; "),
    stringsAsFactors = FALSE
  )
})
cross_rank_rows <- Filter(Negate(is.null), cross_rank_rows)
cross_rank_name_collisions <- do.call(rbind, cross_rank_rows)
cross_rank_name_collisions <- cross_rank_name_collisions[
  order(-cross_rank_name_collisions$RankCount, cross_rank_name_collisions$Label),
  , drop = FALSE
]

rank_top_rows <- list()
rank_group_rows <- list()
top_taxa_by_rank <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_means <- rowMeans(rel_rank[known_taxa, , drop = FALSE])
  top5 <- names(sort(known_means, decreasing = TRUE))[1:5]
  top_taxa_by_rank[[rank_name]] <- top5
  for (j in seq_along(top5)) {
    taxon <- top5[[j]]
    rank_top_rows[[length(rank_top_rows) + 1L]] <- data.frame(
      Rank = rank_name,
      RankIndex = match(rank_name, rank_names),
      TopRank = j,
      Taxon = taxon,
      RankQualifiedTaxon = paste0(rank_prefix[[rank_name]], " · ", taxon),
      OverallMean = mean(rel_rank[taxon, ]),
      OverallPercent = 100 * mean(rel_rank[taxon, ]),
      stringsAsFactors = FALSE
    )
  }
  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    display_taxa <- c(top5, unknown_label)
    for (taxon in display_taxa) {
      values <- rel_rank[taxon, ids]
      rank_group_rows[[length(rank_group_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        RankIndex = match(rank_name, rank_names),
        GroupCode = group_code,
        Group = unname(group_display[group_code]),
        Taxon = taxon,
        RankQualifiedTaxon = paste0(rank_prefix[[rank_name]], " · ", taxon),
        Selection = ifelse(taxon == unknown_label, "Unknown", "Global Top 5 known"),
        MeanRelativeAbundance = mean(values),
        MeanPercent = 100 * mean(values),
        Prevalence = mean(values > 0),
        stringsAsFactors = FALSE
      )
    }
  }
}
rank_top_taxa <- do.call(rbind, rank_top_rows)
rank_taxon_group_summary <- do.call(rbind, rank_group_rows)

top_n_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_means <- rowMeans(rel_rank[known_taxa, , drop = FALSE])
  known_order <- names(sort(known_means, decreasing = TRUE))
  branches <- list(
    `Top 5` = head(known_order, 5L),
    `Top 10` = head(known_order, 10L),
    `Top 20` = head(known_order, 20L),
    `Top 50` = head(known_order, 50L),
    `All known` = known_order
  )
  for (label in names(branches)) {
    selected <- branches[[label]]
    coverage <- colSums(rel_rank[selected, , drop = FALSE])
    top_n_rows[[length(top_n_rows) + 1L]] <- data.frame(
      Rank = rank_name,
      RankIndex = match(rank_name, rank_names),
      NLabel = label,
      SelectedN = length(selected),
      SelectedTaxa = ifelse(label == "All known", "ALL_KNOWN", paste(selected, collapse = "; ")),
      MeanCoverage = mean(coverage),
      MinimumCoverage = min(coverage),
      MaximumCoverage = max(coverage),
      OtherKnownMean = mean(1 - rel_rank[unknown_label, ] - coverage),
      UnknownMean = mean(rel_rank[unknown_label, ]),
      stringsAsFactors = FALSE
    )
  }
}
rank_topn_sensitivity <- do.call(rbind, top_n_rows)

complete_mask <- complete_depth == length(rank_names)
complete_rel <- rowsum(relative_feature[complete_mask, , drop = FALSE], lineage_key[complete_mask], reorder = TRUE)
complete_counts <- rowsum(counts[complete_mask, , drop = FALSE], lineage_key[complete_mask], reorder = TRUE)
complete_means <- rowMeans(complete_rel)
complete_order <- names(sort(complete_means, decreasing = TRUE))
lineage_first_index <- match(complete_order, lineage_key)
complete_lineage_abundance <- data.frame(
  LineageRank = seq_along(complete_order),
  LineageKey = complete_order,
  taxonomy_clean[lineage_first_index, rank_names, drop = FALSE],
  FeatureCount = as.integer(table(lineage_key[complete_mask])[complete_order]),
  EqualSampleMean = as.numeric(complete_means[complete_order]),
  PooledReadFraction = rowSums(complete_counts[complete_order, , drop = FALSE]) / sum(counts),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
top12_lineages <- head(complete_order, 12L)

dominant_rows <- list()
ladder_rows <- list()
for (lineage_rank in seq_along(top12_lineages)) {
  key <- top12_lineages[[lineage_rank]]
  path_index <- match(key, lineage_key)
  path <- unlist(taxonomy_clean[path_index, rank_names, drop = FALSE], use.names = TRUE)
  leaf_values <- complete_rel[key, ]
  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    dominant_rows[[length(dominant_rows) + 1L]] <- data.frame(
      LineageRank = lineage_rank,
      LineageKey = key,
      as.list(path),
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      MeanRelativeAbundance = mean(leaf_values[ids]),
      MeanPercent = 100 * mean(leaf_values[ids]),
      Prevalence = mean(leaf_values[ids] > 0),
      OverallMean = mean(leaf_values),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  for (rank_i in seq_along(rank_names)) {
    rank_name <- rank_names[[rank_i]]
    prefix_match <- rep(TRUE, nrow(taxonomy_clean))
    for (prefix_i in seq_len(rank_i)) {
      prefix_rank <- rank_names[[prefix_i]]
      prefix_match <- prefix_match & taxonomy_clean[[prefix_rank]] == path[[prefix_rank]]
    }
    node_mean <- sum(feature_equal_mean[prefix_match])
    ladder_rows[[length(ladder_rows) + 1L]] <- data.frame(
      LineageRank = lineage_rank,
      LineageKey = key,
      LeafGenus = path[["Genus"]],
      Rank = rank_name,
      RankIndex = rank_i,
      Taxon = path[[rank_name]],
      RankQualifiedTaxon = paste0(rank_prefix[[rank_name]], " · ", path[[rank_name]]),
      NodeEqualSampleMean = node_mean,
      NodePercent = 100 * node_mean,
      LeafEqualSampleMean = mean(leaf_values),
      LeafPercent = 100 * mean(leaf_values),
      stringsAsFactors = FALSE
    )
  }
}
dominant_complete_lineages <- do.call(rbind, dominant_rows)
lineage_ladder <- do.call(rbind, ladder_rows)

aggregation_rows <- list()
denominator_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  counts_rank <- rank_counts[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  top5 <- top_taxa_by_rank[[rank_name]]
  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    equal_vector <- rowMeans(rel_rank[, ids, drop = FALSE])
    pooled_vector <- rowSums(counts_rank[, ids, drop = FALSE]) / sum(counts_rank[, ids, drop = FALSE])
    aggregation_rows[[length(aggregation_rows) + 1L]] <- data.frame(
      Rank = rank_name,
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      TotalVariation = 0.5 * sum(abs(equal_vector - pooled_vector)),
      MaximumAbsoluteTaxonDifferencePP = 100 * max(abs(equal_vector - pooled_vector)),
      EqualSampleKnownFraction = 1 - equal_vector[[unknown_label]],
      PooledReadKnownFraction = 1 - pooled_vector[[unknown_label]],
      KnownFractionDifferencePP = 100 * (
        (1 - equal_vector[[unknown_label]]) - (1 - pooled_vector[[unknown_label]])
      ),
      stringsAsFactors = FALSE
    )
    known_total <- 1 - rel_rank[unknown_label, ids]
    for (taxon in top5) {
      all_read_mean <- mean(rel_rank[taxon, ids])
      known_only_mean <- mean(rel_rank[taxon, ids] / known_total)
      denominator_rows[[length(denominator_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        GroupCode = group_code,
        Group = unname(group_display[group_code]),
        Taxon = taxon,
        AllReadMean = all_read_mean,
        KnownOnlyMean = known_only_mean,
        InflationPP = 100 * (known_only_mean - all_read_mean),
        GroupUnknownMean = mean(rel_rank[unknown_label, ids]),
        stringsAsFactors = FALSE
      )
    }
  }
}
aggregation_sensitivity <- do.call(rbind, aggregation_rows)
denominator_sensitivity <- do.call(rbind, denominator_rows)

closure_audit <- data.frame(
  SampleID = sample_ids,
  FeatureClosure = colSums(relative_feature),
  stringsAsFactors = FALSE
)
for (rank_name in rank_names) {
  closure_audit[[paste0(rank_name, "Closure")]] <- colSums(rank_relative[[rank_name]])
}
closure_columns <- grep("Closure$", names(closure_audit), value = TRUE)
closure_audit$MaximumAbsoluteError <- apply(
  abs(as.matrix(closure_audit[, closure_columns, drop = FALSE]) - 1),
  1L,
  max
)

family_genus_gap <- lineage_gap_audit[
  lineage_gap_audit$ParentRank == "Family" & lineage_gap_audit$ChildRank == "Genus",
  , drop = FALSE
]
complete_genus_mean <- sum(feature_equal_mean[complete_mask])
top12_complete_mean <- sum(complete_means[top12_lineages])
largest_resolution_loss_pp <- max(rank_resolution_audit$ResolutionLossFromParentPP, na.rm = TRUE)

add_check("rank-count", "taxonomy", length(rank_names), 5, length(rank_names) == 5L)
add_check(
  "known-means-monotone", "taxonomy",
  paste(round(rank_resolution_audit$KnownEqualSampleMean, 12), collapse = ","),
  "strictly decreasing", all(diff(rank_resolution_audit$KnownEqualSampleMean) < 0)
)
add_check(
  "largest-resolution-loss-family-genus", "taxonomy",
  rank_resolution_audit$Rank[which.max(rank_resolution_audit$ResolutionLossFromParentPP)],
  "Genus",
  identical(rank_resolution_audit$Rank[which.max(rank_resolution_audit$ResolutionLossFromParentPP)], "Genus")
)
add_check(
  "largest-resolution-loss-value", "taxonomy", largest_resolution_loss_pp,
  31.0683652412, near(largest_resolution_loss_pp, 31.0683652412, 1e-9)
)
for (group_code in group_codes) {
  observed <- group_rank_resolution$KnownMean[
    group_rank_resolution$GroupCode == group_code & group_rank_resolution$Rank == "Genus"
  ]
  expected <- c(IW = 0.289913056827080, CW = 0.383846784456552, TW = 0.433071819506610)[[group_code]]
  add_check(
    paste0("genus-known-", tolower(group_code)), "taxonomy",
    observed, expected, near(observed, expected, 1e-12)
  )
}
add_check(
  "early-parent-child-gaps-zero", "hierarchy",
  sum(head(lineage_gap_audit$ChildKnownParentUnknownFeatures, 3L)), 0,
  sum(head(lineage_gap_audit$ChildKnownParentUnknownFeatures, 3L)) == 0L
)
add_check(
  "family-genus-gap-features", "hierarchy",
  family_genus_gap$ChildKnownParentUnknownFeatures, 324,
  family_genus_gap$ChildKnownParentUnknownFeatures == 324L
)
add_check(
  "family-genus-gap-taxa", "hierarchy",
  family_genus_gap$ChildKnownParentUnknownTaxa, 17,
  family_genus_gap$ChildKnownParentUnknownTaxa == 17L
)
add_check(
  "family-genus-gap-mean", "hierarchy",
  family_genus_gap$ChildKnownParentUnknownEqualSampleMean, 0.0251827731587,
  near(family_genus_gap$ChildKnownParentUnknownEqualSampleMean, 0.0251827731587, 1e-12)
)
add_check(
  "complete-genus-lineages", "hierarchy", nrow(complete_lineage_abundance), 754,
  nrow(complete_lineage_abundance) == 754L
)
add_check(
  "complete-genus-mean", "hierarchy", complete_genus_mean, 0.343761113771,
  near(complete_genus_mean, 0.343761113771, 1e-12)
)
add_check(
  "genus-known-decomposes", "hierarchy",
  complete_genus_mean + family_genus_gap$ChildKnownParentUnknownEqualSampleMean,
  expected_known_means[["Genus"]],
  near(
    complete_genus_mean + family_genus_gap$ChildKnownParentUnknownEqualSampleMean,
    expected_known_means[["Genus"]],
    1e-12
  )
)
add_check("top12-lineages", "hierarchy", length(top12_lineages), 12, length(top12_lineages) == 12L)
add_check(
  "top12-lineage-coverage", "hierarchy", top12_complete_mean, 0.101876289589,
  near(top12_complete_mean, 0.101876289589, 1e-12)
)
add_check(
  "parent-specific-label-collisions", "hierarchy", nrow(label_collision_audit), 2,
  nrow(label_collision_audit) == 2L
)
add_check(
  "parent-specific-collision-mass", "hierarchy", sum(label_collision_audit$EqualSampleMean),
  0.00630608162952,
  near(sum(label_collision_audit$EqualSampleMean), 0.00630608162952, 1e-12)
)
add_check(
  "cross-rank-name-collisions", "hierarchy", nrow(cross_rank_name_collisions), 13,
  nrow(cross_rank_name_collisions) == 13L
)
add_check("rank-top5-rows", "composition", nrow(rank_top_taxa), 25, nrow(rank_top_taxa) == 25L)
add_check(
  "rank-group-bubble-rows", "composition", nrow(rank_taxon_group_summary), 90,
  nrow(rank_taxon_group_summary) == 90L
)
add_check(
  "rank-topn-rows", "sensitivity", nrow(rank_topn_sensitivity), 25,
  nrow(rank_topn_sensitivity) == 25L
)
expected_top10 <- c(
  Phylum = 0.948030470939761,
  Class = 0.590045317291106,
  Order = 0.347652784912770,
  Family = 0.236710656674214,
  Genus = 0.0952180453933754
)
observed_top10 <- setNames(
  rank_topn_sensitivity$MeanCoverage[rank_topn_sensitivity$NLabel == "Top 10"],
  rank_topn_sensitivity$Rank[rank_topn_sensitivity$NLabel == "Top 10"]
)[rank_names]
add_check(
  "top10-rank-coverage", "sensitivity", observed_top10, expected_top10,
  near(observed_top10, expected_top10, 1e-12)
)
add_check(
  "all-rank-closures", "composition", max(closure_audit$MaximumAbsoluteError),
  "<=1e-12", max(closure_audit$MaximumAbsoluteError) <= 1e-12
)
add_check(
  "known-only-inflation-positive", "sensitivity", max(denominator_sensitivity$InflationPP),
  ">0", max(denominator_sensitivity$InflationPP) > 0
)
add_check(
  "known-only-inflation-maximum", "sensitivity", max(denominator_sensitivity$InflationPP),
  7.27296856377, near(max(denominator_sensitivity$InflationPP), 7.27296856377, 1e-9)
)
add_check(
  "aggregation-sensitivity-rows", "sensitivity", nrow(aggregation_sensitivity), 15,
  nrow(aggregation_sensitivity) == 15L
)

rank_colours <- c(
  Phylum = "#0072B2",
  Class = "#E69F00",
  Order = "#009E73",
  Family = "#CC79A7",
  Genus = "#D55E00"
)
group_colours <- c(
  `All samples` = "#1F2937",
  `Inland wetland` = "#E69F00",
  `Coastal wetland` = "#0072B2",
  `Tibetan wetland` = "#009E73"
)

cascade_data <- group_rank_resolution
cascade_data$Rank <- factor(cascade_data$Rank, levels = rank_names)
cascade_data$Group <- factor(
  cascade_data$Group,
  levels = c("All samples", unname(group_display[group_codes]))
)
cascade_plot <- ggplot2::ggplot(
  cascade_data,
  ggplot2::aes(
    x = Rank,
    y = KnownMean,
    colour = Group,
    linetype = Group,
    group = Group
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(
    values = c(
      `All samples` = "22",
      `Inland wetland` = "solid",
      `Coastal wetland` = "solid",
      `Tibetan wetland` = "solid"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Taxonomic resolution falls as the requested rank becomes finer",
    subtitle = "Known fraction uses all reads and equal sample weight within each wetland group",
    x = NULL,
    y = "Reads assigned at rank",
    colour = NULL,
    linetype = NULL,
    caption = paste(
      "Genus-level resolution: inland 29.0%, coastal 38.4%, Tibetan 43.3%;\n",
      "differences in annotation completeness remain visible."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

bubble_data <- rank_taxon_group_summary
bubble_data$Rank <- factor(bubble_data$Rank, levels = rank_names)
bubble_data$Group <- factor(
  bubble_data$Group,
  levels = unname(group_display[group_codes])
)
bubble_levels <- unlist(lapply(rank_names, function(rank_name) {
  c(
    paste0(rank_prefix[[rank_name]], " · ", top_taxa_by_rank[[rank_name]]),
    paste0(rank_prefix[[rank_name]], " · ", unknown_labels[[rank_name]])
  )
}), use.names = FALSE)
bubble_data$RankQualifiedTaxon <- factor(
  bubble_data$RankQualifiedTaxon,
  levels = rev(bubble_levels)
)
bubble_plot <- ggplot2::ggplot(
  bubble_data,
  ggplot2::aes(
    x = Group,
    y = RankQualifiedTaxon,
    size = MeanPercent,
    colour = Prevalence
  )
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::facet_wrap(~Rank, ncol = 2, scales = "free_y") +
  ggplot2::scale_size_area(max_size = 11, breaks = c(1, 5, 20, 60)) +
  ggplot2::scale_colour_gradientn(
    colours = c("#D9EAF7", "#56B4E9", "#0072B2", "#003B5C"),
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      `Inland wetland` = "IW",
      `Coastal wetland` = "CW",
      `Tibetan wetland` = "TW"
    )
  ) +
  ggplot2::labs(
    title = "The same Top-5 rule yields different views at each rank",
    subtitle = "Global Top 5 named taxa plus rank-specific Unknown; one fixed list is reused across groups",
    x = "Wetland group (IW / CW / TW)",
    y = NULL,
    size = "Mean abundance (%)",
    colour = "Prevalence",
    caption = paste(
      "IW, inland; CW, coastal; TW, Tibetan wetland. P/C/O/F/G prefixes prevent",
      "equal-looking names at different ranks from being merged."
    )
  ) +
  theme_pub(base_size = 7.7) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.25),
    panel.spacing.x = grid::unit(1.2, "mm"),
    strip.text = ggplot2::element_text(size = 8)
  )

ladder_data <- lineage_ladder
ladder_data$Rank <- factor(ladder_data$Rank, levels = rank_names)
leaf_order <- unique(ladder_data$LeafGenus[order(ladder_data$LineageRank)])
ladder_data$LeafGenus <- factor(ladder_data$LeafGenus, levels = rev(leaf_order))
ladder_data$CellLabel <- wrap_label(ladder_data$Taxon, width = 15L)
ladder_dark <- ladder_data$NodePercent >= 9
ladder_plot <- ggplot2::ggplot(
  ladder_data,
  ggplot2::aes(x = Rank, y = LeafGenus, fill = NodePercent)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(
    data = ladder_data[!ladder_dark, , drop = FALSE],
    ggplot2::aes(label = CellLabel),
    colour = "#17202A",
    family = font_family,
    size = 2.05,
    lineheight = 0.88
  ) +
  ggplot2::geom_text(
    data = ladder_data[ladder_dark, , drop = FALSE],
    ggplot2::aes(label = CellLabel),
    colour = "white",
    family = font_family,
    size = 2.05,
    lineheight = 0.88
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    trans = "sqrt",
    labels = scales::label_number(accuracy = 0.1, suffix = "%")
  ) +
  ggplot2::scale_y_discrete(labels = function(x) wrap_label(x, width = 21L)) +
  ggplot2::labs(
    title = "A genus is a complete path, not an isolated display name",
    subtitle = "Top 12 complete genus lineages cover 10.19% of all reads; tile colour is prefix-node abundance",
    x = NULL,
    y = "Leaf genus",
    fill = "Node abundance",
    caption = paste(
      "Each row traces one Phylum-to-Genus path. Prefix-node percentages repeat",
      "context and must not be added across columns."
    )
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.ticks = ggplot2::element_blank()
  )

topn_plot_data <- rank_topn_sensitivity
topn_plot_data$Rank <- factor(topn_plot_data$Rank, levels = rank_names)
topn_plot_data$NLabel <- factor(
  topn_plot_data$NLabel,
  levels = c("Top 5", "Top 10", "Top 20", "Top 50", "All known")
)
top10_labels <- topn_plot_data[topn_plot_data$NLabel == "Top 10", , drop = FALSE]
top10_labels$Label <- scales::percent(top10_labels$MeanCoverage, accuracy = 0.1)
topn_plot <- ggplot2::ggplot(
  topn_plot_data,
  ggplot2::aes(x = NLabel, y = MeanCoverage, colour = Rank, group = Rank)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::geom_text(
    data = top10_labels,
    ggplot2::aes(label = Label),
    nudge_x = 0.14,
    hjust = 0,
    size = 2.4,
    show.legend = FALSE,
    family = font_family
  ) +
  ggplot2::scale_colour_manual(values = rank_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.03),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Top 10 does not mean the same coverage at every rank",
    subtitle = "Coverage is measured against all reads after per-sample closure",
    x = "Display branch",
    y = "Mean covered abundance",
    colour = "Rank",
    caption = paste(
      "Top-10 coverage falls from 94.8% at phylum to 9.5% at genus;\n",
      "All known retains the rank-specific unclassified fraction."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

figure_specification <- data.frame(
  stem = c(
    "27-rank-resolution-cascade",
    "27-multirank-bubble",
    "27-lineage-ladder",
    "27-topn-coverage"
  ),
  width_mm = c(150, 183, 183, 160),
  height_mm = c(96, 198, 168, 102),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(cascade_plot, bubble_plot, ladder_plot, topn_plot)
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
write_tsv(taxonomy_lineage_map, "taxonomy-lineage-map.tsv")
write_tsv(multirank_relative_abundance, "multirank-relative-abundance.tsv")
write_tsv(rank_resolution_audit, "rank-resolution-audit.tsv")
write_tsv(group_rank_resolution, "group-rank-resolution.tsv")
write_tsv(lineage_gap_audit, "lineage-gap-audit.tsv")
write_tsv(complete_depth_audit, "complete-depth-audit.tsv")
write_tsv(label_collision_audit, "label-collision-audit.tsv")
write_tsv(cross_rank_name_collisions, "cross-rank-name-collisions.tsv")
write_tsv(rank_top_taxa, "rank-top-taxa.tsv")
write_tsv(rank_taxon_group_summary, "rank-taxon-group-summary.tsv")
write_tsv(rank_topn_sensitivity, "rank-topn-sensitivity.tsv")
write_tsv(complete_lineage_abundance, "complete-lineage-abundance.tsv")
write_tsv(dominant_complete_lineages, "dominant-complete-lineages.tsv")
write_tsv(lineage_ladder, "lineage-ladder.tsv")
write_tsv(aggregation_sensitivity, "aggregation-sensitivity.tsv")
write_tsv(denominator_sensitivity, "denominator-sensitivity.tsv")
write_tsv(closure_audit, "closure-audit.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

genus_group_known <- setNames(
  group_rank_resolution$KnownMean[
    group_rank_resolution$Rank == "Genus" &
      group_rank_resolution$GroupCode %in% group_codes
  ],
  group_rank_resolution$GroupCode[
    group_rank_resolution$Rank == "Genus" &
      group_rank_resolution$GroupCode %in% group_codes
  ]
)[group_codes]

summary_payload <- list(
  article = 27L,
  title = "Multirank composition from phylum to genus",
  dataset = list(
    source = "microeco 2.0.0 packaged Chinese wetland 16S tables",
    features = nrow(counts),
    samples = ncol(counts),
    reads = sum(counts),
    groups = length(group_codes),
    library_depth_min = min(sample_depths),
    library_depth_max = max(sample_depths)
  ),
  taxonomy = list(
    ranks = length(rank_names),
    rank_names = rank_names,
    labels_by_rank = as.list(setNames(rank_resolution_audit$Labels, rank_names)),
    known_equal_sample_mean = as.list(setNames(
      rank_resolution_audit$KnownEqualSampleMean,
      rank_names
    )),
    genus_known_equal_sample_mean = expected_known_means[["Genus"]],
    genus_known_by_group = as.list(genus_group_known),
    largest_resolution_loss_transition = "Family to Genus",
    largest_resolution_loss_pp = largest_resolution_loss_pp,
    incertae_sedis_features = sum(rank_resolution_audit$IncertaeSedisFeatures)
  ),
  hierarchy = list(
    complete_genus_lineages = nrow(complete_lineage_abundance),
    complete_genus_lineage_mean = complete_genus_mean,
    family_unknown_genus_known_features = family_genus_gap$ChildKnownParentUnknownFeatures,
    family_unknown_genus_known_taxa = family_genus_gap$ChildKnownParentUnknownTaxa,
    family_unknown_genus_known_mean = family_genus_gap$ChildKnownParentUnknownEqualSampleMean,
    parent_specific_label_collisions = nrow(label_collision_audit),
    parent_specific_collision_mean = sum(label_collision_audit$EqualSampleMean),
    cross_rank_name_collisions = nrow(cross_rank_name_collisions),
    dominant_complete_lineages = length(top12_lineages),
    dominant_complete_lineage_mean = top12_complete_mean,
    rank_qualified_labels = TRUE,
    missing_parent_backfill = FALSE
  ),
  primary = list(
    ranks = rank_names,
    top_known_per_rank = 5L,
    aggregation = "Equal-sample mean after per-sample closure",
    denominator = "All reads",
    unknown_separate = TRUE,
    inferential_claim = FALSE
  ),
  sensitivity = list(
    top_n_rows = nrow(rank_topn_sensitivity),
    top_n_branches_per_rank = 5L,
    top10_coverage_by_rank = as.list(observed_top10),
    maximum_known_only_inflation_pp = max(denominator_sensitivity$InflationPP),
    maximum_equal_vs_pooled_total_variation = max(aggregation_sensitivity$TotalVariation),
    maximum_equal_vs_pooled_taxon_difference_pp = max(
      aggregation_sensitivity$MaximumAbsoluteTaxonDifferencePP
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
  file.path(output_dir, "multirank-composition-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

session_lines <- sanitize_text(capture.output(sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))

validation_lines <- c(
  "Article 27 multirank-composition validation",
  paste0("Status: ", ifelse(checks_failed == 0L, "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total),
  "Samples/features/reads: 90/13628/1619670",
  paste0(
    "Known fractions Phylum-to-Genus: ",
    paste(sprintf("%.4f", rank_resolution_audit$KnownEqualSampleMean), collapse = "/")
  ),
  paste0("Largest resolution loss: ", sprintf("%.3f pp", largest_resolution_loss_pp)),
  paste0("Complete genus lineage mean: ", sprintf("%.4f", complete_genus_mean)),
  paste0("Family-unknown/genus-known features: ", family_genus_gap$ChildKnownParentUnknownFeatures),
  paste0("Cross-rank name collisions: ", nrow(cross_rank_name_collisions)),
  paste0("Top-12 complete lineage coverage: ", sprintf("%.4f", top12_complete_mean)),
  paste0("Maximum known-only inflation: ", sprintf("%.3f pp", max(denominator_sensitivity$InflationPP))),
  paste0("Primary exports: ", nrow(format_audit), " files"),
  paste0("R version: ", r_version)
)
writeLines(validation_lines, file.path(output_dir, "validation.log"))

if (checks_failed > 0L) {
  failed_ids <- validation_checks$check_id[validation_checks$status == "FAIL"]
  stop(
    "Article 27 validation failed: ",
    paste(failed_ids, collapse = ", "),
    call. = FALSE
  )
}

cat(
  sprintf(
    "Article 27 validation passed: %d/%d checks.\n",
    checks_passed,
    checks_total
  )
)
