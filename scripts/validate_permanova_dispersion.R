#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260723)

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
project_root <- normalizePath(
  as_absolute(args[["project-root"]]),
  mustWork = TRUE
)
input_dir <- normalizePath(
  as_absolute(args[["input-dir"]]),
  mustWork = TRUE
)
output_dir <- as_absolute(args[["output-dir"]])
figure_dir <- as_absolute(args[["figure-dir"]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "digest", "ggplot2", "jsonlite", "permute", "ragg", "readr",
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
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  source_summary = file.path(input_dir, "source-summary.json"),
  prepare_script = file.path(
    project_root,
    "scripts",
    "prepare_permanova_dispersion_data.R"
  ),
  theme_pub = file.path(project_root, "R", "theme_pub.R")
)
expected_sha256 <- c(
  otutab = "de52dc7457f15f07776c216561803afd92782af7d6723ec93deb8dc0e932b900",
  taxonomy = "599ed73e3d1f59bc2867782c58113100def184c3ebd2cab36286b2fc96dbaee0",
  metadata = "5ebf2fd60582e1a0ef80fe20d251594bec3607df11d8bb74d0981236e03f982c",
  source_summary = "0b7f7c5a6db0c3d9c09d9ae147ff926cb1ce8cbbed804290c074c088100e076d",
  prepare_script = "dc6ff8650872d9693afc3d8f5965dd5d21446929e58ef97b0ddb703156815cf4",
  theme_pub = "8d3821a485aeb529b12184bbe977c1ca952131a22d5eb69b802088ce91909beb"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 23 input(s): ",
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
    "PERMANOVA PERMDISP Bray Curtis split plot Warming Clipping",
    "Ambient Block spatial median pseudo F Holm 0123456789"
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

otu_frame <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
source_summary <- jsonlite::read_json(input_paths[["source_summary"]])

counts_character <- as.matrix(otu_frame)
counts <- matrix(
  suppressWarnings(as.numeric(counts_character)),
  nrow = nrow(counts_character),
  ncol = ncol(counts_character),
  dimnames = dimnames(counts_character)
)
feature_ids <- rownames(counts)
library_ids <- colnames(counts)
rank_columns <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
required_metadata <- c(
  "SourceTreatment", "BiologicalSampleID", "Block", "MainPlot",
  "Warming", "Clipping", "Treatment", "TechnicalLibraryIndex",
  "TechnicalLibraryCount"
)

add_check("feature-count", "data", nrow(counts), 16825, nrow(counts) == 16825L)
add_check("library-count", "data", ncol(counts), 56, ncol(counts) == 56L)
add_check("read-count", "data", sum(counts), 98022, sum(counts) == 98022)
add_check(
  "taxonomy-feature-count",
  "data",
  nrow(taxonomy),
  16825,
  nrow(taxonomy) == 16825L
)
add_check(
  "metadata-library-count",
  "data",
  nrow(metadata),
  56,
  nrow(metadata) == 56L
)
add_check("feature-ids-unique", "data", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("library-ids-unique", "data", anyDuplicated(library_ids), 0, !anyDuplicated(library_ids))
add_check(
  "taxonomy-feature-ids",
  "data",
  length(intersect(feature_ids, rownames(taxonomy))),
  length(feature_ids),
  setequal(feature_ids, rownames(taxonomy))
)
add_check(
  "metadata-library-ids",
  "data",
  length(intersect(library_ids, rownames(metadata))),
  length(library_ids),
  setequal(library_ids, rownames(metadata))
)
add_check("counts-finite", "data", sum(is.finite(counts)), length(counts), all(is.finite(counts)))
add_check("counts-nonnegative", "data", min(counts), ">=0", all(counts >= 0))
add_check(
  "counts-integer",
  "data",
  max(abs(counts - round(counts))),
  0,
  all(abs(counts - round(counts)) < .Machine$double.eps^0.5)
)
add_check("nonempty-features", "data", sum(rowSums(counts) > 0), nrow(counts), all(rowSums(counts) > 0))
add_check("nonempty-libraries", "data", sum(colSums(counts) > 0), ncol(counts), all(colSums(counts) > 0))
add_check(
  "taxonomy-rank-columns",
  "data",
  intersect(rank_columns, colnames(taxonomy)),
  rank_columns,
  all(rank_columns %in% colnames(taxonomy))
)
rank_values_empty <- all(vapply(
  taxonomy[rank_columns],
  function(x) all(is.na(x) | as.character(x) == ""),
  logical(1)
))
add_check(
  "taxonomy-ranks-explicitly-empty",
  "data",
  rank_values_empty,
  TRUE,
  rank_values_empty
)
taxonomy_boundary <- all(
  taxonomy$TaxonomyStatus ==
    "Not distributed with phyloseq::soilrep; not used in this analysis"
)
add_check(
  "taxonomy-boundary-declared",
  "data",
  taxonomy_boundary,
  TRUE,
  taxonomy_boundary
)
add_check(
  "metadata-columns",
  "data",
  intersect(required_metadata, colnames(metadata)),
  required_metadata,
  all(required_metadata %in% colnames(metadata))
)
add_check(
  "source-object-checksum",
  "provenance",
  source_summary$source_sha256,
  "47b25a0a033b794fca1db30ae6fdb68c55a44e6c4cd0d1be76788d1d77b51f36",
  identical(
    source_summary$source_sha256,
    "47b25a0a033b794fca1db30ae6fdb68c55a44e6c4cd0d1be76788d1d77b51f36"
  )
)
add_check(
  "source-study-doi",
  "provenance",
  source_summary$primary_study$doi,
  "10.1038/ismej.2011.11",
  identical(source_summary$primary_study$doi, "10.1038/ismej.2011.11")
)

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[library_ids, , drop = FALSE]
treatment_levels <- c(
  "Ambient + Unclipped",
  "Ambient + Clipped",
  "Warmed + Unclipped",
  "Warmed + Clipped"
)
metadata$Block <- factor(metadata$Block, levels = paste("Block", 1:6))
metadata$Warming <- factor(
  metadata$Warming,
  levels = c("Ambient", "Warmed")
)
metadata$Clipping <- factor(
  metadata$Clipping,
  levels = c("Unclipped", "Clipped")
)
metadata$Treatment <- factor(
  metadata$Treatment,
  levels = treatment_levels
)

add_check(
  "biological-sample-count-source",
  "design",
  length(unique(metadata$BiologicalSampleID)),
  24,
  length(unique(metadata$BiologicalSampleID)) == 24L
)
add_check(
  "block-count",
  "design",
  length(unique(metadata$Block)),
  6,
  length(unique(metadata$Block)) == 6L
)
add_check(
  "main-plot-count",
  "design",
  length(unique(metadata$MainPlot)),
  12,
  length(unique(metadata$MainPlot)) == 12L
)
add_check(
  "treatment-library-balance",
  "design",
  as.integer(table(metadata$Treatment)),
  rep(14L, 4L),
  identical(as.integer(table(metadata$Treatment)), rep(14L, 4L))
)

biological_ids <- unique(as.character(metadata$BiologicalSampleID))
biological_meta <- do.call(
  rbind,
  lapply(biological_ids, function(id) {
    rows <- metadata$BiologicalSampleID == id
    fields <- c("Block", "MainPlot", "Warming", "Clipping", "Treatment")
    stopifnot(all(vapply(metadata[rows, fields, drop = FALSE], function(x) {
      length(unique(as.character(x))) == 1L
    }, logical(1))))
    data.frame(
      BiologicalSampleID = id,
      Block = as.character(metadata$Block[rows][[1L]]),
      MainPlot = as.character(metadata$MainPlot[rows][[1L]]),
      Warming = as.character(metadata$Warming[rows][[1L]]),
      Clipping = as.character(metadata$Clipping[rows][[1L]]),
      Treatment = as.character(metadata$Treatment[rows][[1L]]),
      TechnicalLibraries = sum(rows),
      AggregateReads = sum(counts[, rows, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  })
)
biological_meta$Block <- factor(
  biological_meta$Block,
  levels = paste("Block", 1:6)
)
biological_meta$Warming <- factor(
  biological_meta$Warming,
  levels = c("Ambient", "Warmed")
)
biological_meta$Clipping <- factor(
  biological_meta$Clipping,
  levels = c("Unclipped", "Clipped")
)
biological_meta$Treatment <- factor(
  biological_meta$Treatment,
  levels = treatment_levels
)
biological_meta <- biological_meta[
  order(
    biological_meta$Block,
    biological_meta$Warming,
    biological_meta$Clipping
  ),
  ,
  drop = FALSE
]
rownames(biological_meta) <- biological_meta$BiologicalSampleID
biological_ids <- rownames(biological_meta)

aggregate_counts <- vapply(
  biological_ids,
  function(id) {
    rowSums(counts[, metadata$BiologicalSampleID == id, drop = FALSE])
  },
  numeric(nrow(counts))
)
rownames(aggregate_counts) <- feature_ids
colnames(aggregate_counts) <- biological_ids

design_table <- with(
  biological_meta,
  table(Block, Warming, Clipping)
)
add_check("aggregated-biological-samples", "aggregation", ncol(aggregate_counts), 24, ncol(aggregate_counts) == 24L)
add_check("aggregation-read-conservation", "aggregation", sum(aggregate_counts), sum(counts), sum(aggregate_counts) == sum(counts))
add_check(
  "split-plot-cell-balance",
  "design",
  as.integer(design_table),
  rep(1L, 24L),
  all(design_table == 1L)
)
add_check(
  "technical-library-count-range",
  "aggregation",
  range(biological_meta$TechnicalLibraries),
  c(1, 3),
  identical(as.integer(range(biological_meta$TechnicalLibraries)), c(1L, 3L))
)
add_check(
  "biological-depth-positive",
  "aggregation",
  min(colSums(aggregate_counts)),
  ">0",
  all(colSums(aggregate_counts) > 0)
)

aggregation_audit <- data.frame(
  LibraryID = library_ids,
  BiologicalSampleID = metadata$BiologicalSampleID,
  Block = as.character(metadata$Block),
  MainPlot = metadata$MainPlot,
  Treatment = as.character(metadata$Treatment),
  LibraryReads = as.numeric(colSums(counts)),
  TechnicalLibraryIndex = as.integer(metadata$TechnicalLibraryIndex),
  TechnicalLibraryCount = as.integer(metadata$TechnicalLibraryCount),
  InferenceUnit = "BiologicalSampleID after count aggregation",
  stringsAsFactors = FALSE
)
design_audit <- data.frame(
  BiologicalSampleID = biological_meta$BiologicalSampleID,
  Block = as.character(biological_meta$Block),
  MainPlot = biological_meta$MainPlot,
  Warming = as.character(biological_meta$Warming),
  Clipping = as.character(biological_meta$Clipping),
  Treatment = as.character(biological_meta$Treatment),
  TechnicalLibraries = biological_meta$TechnicalLibraries,
  AggregateReads = biological_meta$AggregateReads,
  InferentialContribution = 1L,
  stringsAsFactors = FALSE
)

# The primary prevalence rule mirrors the study's Figure 2 preprocessing:
# retain a feature detected in at least four tagged libraries in any treatment.
presence_by_treatment <- vapply(
  treatment_levels,
  function(group) {
    rowSums(counts[, metadata$Treatment == group, drop = FALSE] > 0)
  },
  numeric(nrow(counts))
)
colnames(presence_by_treatment) <- treatment_levels
primary_keep <- apply(presence_by_treatment >= 4L, 1L, any)
overall_three_keep <- rowSums(aggregate_counts > 0) >= 3L

feature_filter_audit <- data.frame(
  Branch = c(
    "Primary: >=4 tagged libraries in any treatment",
    "Sensitivity: no prevalence filter",
    "Descriptive: >=3 biological samples overall"
  ),
  InputFeatures = nrow(counts),
  RetainedFeatures = c(
    sum(primary_keep),
    nrow(counts),
    sum(overall_three_keep)
  ),
  ThresholdDefinedBeforeTesting = c(TRUE, TRUE, TRUE),
  UsedForPrimaryInference = c(TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
add_check("primary-filter-features", "standardization", sum(primary_keep), 1824, sum(primary_keep) == 1824L)
add_check(
  "primary-filter-nonempty",
  "standardization",
  sum(rowSums(aggregate_counts[primary_keep, , drop = FALSE]) > 0),
  sum(primary_keep),
  all(rowSums(aggregate_counts[primary_keep, , drop = FALSE]) > 0)
)

primary_counts <- aggregate_counts[primary_keep, , drop = FALSE]
primary_community <- t(primary_counts)
primary_relative <- sweep(
  primary_community,
  1L,
  rowSums(primary_community),
  "/"
)
add_check(
  "relative-abundance-row-sums",
  "standardization",
  max(abs(rowSums(primary_relative) - 1)),
  "<1e-12",
  max(abs(rowSums(primary_relative) - 1)) < 1e-12
)
primary_distance <- vegan::vegdist(primary_relative, method = "bray")
distance_matrix <- as.matrix(primary_distance)
add_check("distance-finite", "distance", sum(is.finite(distance_matrix)), length(distance_matrix), all(is.finite(distance_matrix)))
add_check("distance-symmetric", "distance", max(abs(distance_matrix - t(distance_matrix))), "<1e-12", max(abs(distance_matrix - t(distance_matrix))) < 1e-12)
add_check("distance-diagonal-zero", "distance", max(abs(diag(distance_matrix))), "<1e-12", max(abs(diag(distance_matrix))) < 1e-12)
add_check("distance-sample-order", "distance", rownames(distance_matrix), biological_ids, identical(rownames(distance_matrix), biological_ids))

make_splitplot_permutations <- function(meta, nset, seed) {
  block_levels <- levels(meta$Block)
  total_possible <- 8^length(block_levels)
  if (nset > total_possible - 1L) {
    stop("Requested more unique non-identity permutations than available.")
  }
  set.seed(seed)
  codes <- sample.int(total_possible - 1L, nset, replace = FALSE)
  permutations <- matrix(
    rep(seq_len(nrow(meta)), each = nset),
    nrow = nset,
    byrow = FALSE
  )
  for (row_i in seq_len(nset)) {
    code <- codes[[row_i]]
    for (block_i in seq_along(block_levels)) {
      state <- code %% 8L
      code <- code %/% 8L
      destination <- which(meta$Block == block_levels[[block_i]])
      ambient <- destination[meta$Warming[destination] == "Ambient"]
      warmed <- destination[meta$Warming[destination] == "Warmed"]
      ambient <- ambient[order(meta$Clipping[ambient])]
      warmed <- warmed[order(meta$Clipping[warmed])]
      if (bitwAnd(state, 2L) != 0L) ambient <- rev(ambient)
      if (bitwAnd(state, 4L) != 0L) warmed <- rev(warmed)
      source <- if (bitwAnd(state, 1L) != 0L) {
        c(warmed, ambient)
      } else {
        c(ambient, warmed)
      }
      permutations[row_i, destination] <- source
    }
  }
  list(
    matrix = permutations,
    codes = codes,
    possible = total_possible
  )
}

is_legal_splitplot <- function(permutation, meta) {
  for (block_level in levels(meta$Block)) {
    destination <- which(meta$Block == block_level)
    source <- permutation[destination]
    if (!all(meta$Block[source] == block_level)) return(FALSE)
    destination_ambient <- destination[meta$Warming[destination] == "Ambient"]
    destination_warmed <- destination[meta$Warming[destination] == "Warmed"]
    source_plot_a <- unique(as.character(meta$MainPlot[permutation[destination_ambient]]))
    source_plot_w <- unique(as.character(meta$MainPlot[permutation[destination_warmed]]))
    if (length(source_plot_a) != 1L || length(source_plot_w) != 1L) return(FALSE)
    if (identical(source_plot_a, source_plot_w)) return(FALSE)
    if (!setequal(meta$Clipping[permutation[destination_ambient]], levels(meta$Clipping))) {
      return(FALSE)
    }
    if (!setequal(meta$Clipping[permutation[destination_warmed]], levels(meta$Clipping))) {
      return(FALSE)
    }
  }
  TRUE
}

primary_seed <- 20260723L
primary_nperm <- 9999L
splitplot <- make_splitplot_permutations(
  biological_meta,
  nset = primary_nperm,
  seed = primary_seed
)
permutation_matrix <- splitplot$matrix
control_design <- permute::how(
  blocks = biological_meta$Block,
  plots = permute::Plots(
    strata = factor(biological_meta$MainPlot),
    type = "free"
  ),
  within = permute::Within(type = "free"),
  nperm = primary_nperm
)
possible_from_permute <- permute::numPerms(
  nrow(biological_meta),
  control = control_design
)
identity_rows <- sum(apply(
  permutation_matrix,
  1L,
  function(x) identical(as.integer(x), seq_len(nrow(biological_meta)))
))
legal_rows <- vapply(
  seq_len(nrow(permutation_matrix)),
  function(i) is_legal_splitplot(permutation_matrix[i, ], biological_meta),
  logical(1)
)
minimum_primary_p <- 1 / (primary_nperm + 1)

add_check("permutation-space-manual", "permutation", splitplot$possible, 262144, splitplot$possible == 262144)
add_check("permutation-space-permute", "permutation", possible_from_permute, 262144, possible_from_permute == 262144)
add_check("permutation-count", "permutation", nrow(permutation_matrix), 9999, nrow(permutation_matrix) == 9999L)
add_check("permutation-width", "permutation", ncol(permutation_matrix), 24, ncol(permutation_matrix) == 24L)
add_check("permutation-unique", "permutation", nrow(unique(permutation_matrix)), 9999, nrow(unique(permutation_matrix)) == 9999L)
add_check("permutation-no-identity", "permutation", identity_rows, 0, identity_rows == 0L)
add_check("permutation-all-legal", "permutation", sum(legal_rows), 9999, all(legal_rows))
add_check("permutation-minimum-p", "permutation", minimum_primary_p, 0.0001, identical(minimum_primary_p, 0.0001))

permutation_audit <- data.frame(
  Scheme = "Primary split-plot Monte Carlo",
  Seed = primary_seed,
  Units = nrow(biological_meta),
  Blocks = nlevels(biological_meta$Block),
  MainPlots = length(unique(biological_meta$MainPlot)),
  PossibleIncludingIdentity = splitplot$possible,
  NonIdentityAvailable = splitplot$possible - 1L,
  PermutationsUsed = nrow(permutation_matrix),
  UniqueRows = nrow(unique(permutation_matrix)),
  IdentityRows = identity_rows,
  LegalRows = sum(legal_rows),
  MinimumP = minimum_primary_p,
  MatrixSHA256 = digest::digest(
    permutation_matrix,
    algo = "sha256",
    serialize = TRUE
  ),
  stringsAsFactors = FALSE
)

biological_meta$Treatment <- droplevels(biological_meta$Treatment)
set.seed(primary_seed)
omnibus_fit <- vegan::adonis2(
  primary_distance ~ Block + Treatment,
  data = biological_meta,
  permutations = permutation_matrix,
  by = "terms",
  parallel = 1
)
omnibus_residual_ss <- omnibus_fit["Residual", "SumOfSqs"]
omnibus_ss <- omnibus_fit["Treatment", "SumOfSqs"]
permanova_omnibus <- data.frame(
  Term = "Treatment (4 levels)",
  Df = omnibus_fit["Treatment", "Df"],
  SumOfSquares = omnibus_ss,
  ConditioningBlockR2 = omnibus_fit["Block", "R2"],
  R2Total = omnibus_fit["Treatment", "R2"],
  R2Partial = omnibus_ss / (omnibus_ss + omnibus_residual_ss),
  PseudoF = omnibus_fit["Treatment", "F"],
  PValue = omnibus_fit["Treatment", "Pr(>F)"],
  Permutations = primary_nperm,
  MinimumP = minimum_primary_p,
  PermutationScheme = "Block + whole-plot + within-main-plot",
  stringsAsFactors = FALSE
)

set.seed(primary_seed)
factorial_fit <- vegan::adonis2(
  primary_distance ~ Block + Warming * Clipping,
  data = biological_meta,
  permutations = permutation_matrix,
  by = "terms",
  parallel = 1
)
factorial_terms <- c("Warming", "Clipping", "Warming:Clipping")
factorial_residual_ss <- factorial_fit["Residual", "SumOfSqs"]
permanova_factorial <- data.frame(
  Term = factorial_terms,
  Df = factorial_fit[factorial_terms, "Df"],
  SumOfSquares = factorial_fit[factorial_terms, "SumOfSqs"],
  R2Total = factorial_fit[factorial_terms, "R2"],
  R2Partial = factorial_fit[factorial_terms, "SumOfSqs"] /
    (factorial_fit[factorial_terms, "SumOfSqs"] + factorial_residual_ss),
  PseudoF = factorial_fit[factorial_terms, "F"],
  PValue = factorial_fit[factorial_terms, "Pr(>F)"],
  Permutations = primary_nperm,
  MinimumP = minimum_primary_p,
  stringsAsFactors = FALSE
)
permanova_factorial$PAdjustedHolm <- stats::p.adjust(
  permanova_factorial$PValue,
  method = "holm"
)
permanova_factorial$RejectHolm05 <-
  permanova_factorial$PAdjustedHolm < 0.05

swapped_fit <- vegan::adonis2(
  primary_distance ~ Block + Clipping * Warming,
  data = biological_meta,
  permutations = permutation_matrix[seq_len(2L), , drop = FALSE],
  by = "terms",
  parallel = 1
)
order_parity <- c(
  Warming = factorial_fit["Warming", "SumOfSqs"] -
    swapped_fit["Warming", "SumOfSqs"],
  Clipping = factorial_fit["Clipping", "SumOfSqs"] -
    swapped_fit["Clipping", "SumOfSqs"],
  Interaction = factorial_fit["Warming:Clipping", "SumOfSqs"] -
    swapped_fit["Clipping:Warming", "SumOfSqs"]
)
add_check(
  "factorial-omnibus-ss-parity",
  "permanova",
  sum(permanova_factorial$SumOfSquares) - permanova_omnibus$SumOfSquares,
  "<1e-10",
  abs(sum(permanova_factorial$SumOfSquares) -
    permanova_omnibus$SumOfSquares) < 1e-10
)
add_check(
  "orthogonal-term-order-parity",
  "permanova",
  max(abs(order_parity)),
  "<1e-10",
  max(abs(order_parity)) < 1e-10
)
add_check(
  "omnibus-r2-range",
  "permanova",
  permanova_omnibus$R2Total,
  "0..1",
  permanova_omnibus$R2Total >= 0 && permanova_omnibus$R2Total <= 1
)
add_check(
  "conditioning-block-r2-range",
  "permanova",
  permanova_omnibus$ConditioningBlockR2,
  "0..1",
  permanova_omnibus$ConditioningBlockR2 >= 0 &&
    permanova_omnibus$ConditioningBlockR2 <= 1
)
add_check(
  "omnibus-p-grid",
  "permanova",
  permanova_omnibus$PValue * (primary_nperm + 1),
  "integer",
  abs(permanova_omnibus$PValue * (primary_nperm + 1) -
    round(permanova_omnibus$PValue * (primary_nperm + 1))) < 1e-8
)
add_check(
  "factorial-statistics-finite",
  "permanova",
  sum(is.finite(as.matrix(permanova_factorial[c(
    "SumOfSquares", "R2Total", "R2Partial", "PseudoF", "PValue"
  )]))),
  15,
  all(is.finite(as.matrix(permanova_factorial[c(
    "SumOfSquares", "R2Total", "R2Partial", "PseudoF", "PValue"
  )])))
)
add_check(
  "factorial-holm-monotonic",
  "permanova",
  all(permanova_factorial$PAdjustedHolm >= permanova_factorial$PValue),
  TRUE,
  all(permanova_factorial$PAdjustedHolm >= permanova_factorial$PValue)
)
add_check(
  "factorial-holm-rejections",
  "permanova",
  sum(permanova_factorial$RejectHolm05),
  0,
  sum(permanova_factorial$RejectHolm05) == 0L
)

make_paired_permutations <- function(meta) {
  block_levels <- levels(meta$Block)
  permutations <- matrix(
    rep(seq_len(nrow(meta)), each = 2^length(block_levels) - 1L),
    nrow = 2^length(block_levels) - 1L,
    byrow = FALSE
  )
  for (code in seq_len(2^length(block_levels) - 1L)) {
    for (block_i in seq_along(block_levels)) {
      if (bitwAnd(code, bitwShiftL(1L, block_i - 1L)) != 0L) {
        destination <- which(meta$Block == block_levels[[block_i]])
        permutations[code, destination] <- rev(destination)
      }
    }
  }
  permutations
}

planned_contrasts <- list(
  list(
    label = "Warming | Unclipped",
    group_a = "Ambient + Unclipped",
    group_b = "Warmed + Unclipped"
  ),
  list(
    label = "Warming | Clipped",
    group_a = "Ambient + Clipped",
    group_b = "Warmed + Clipped"
  ),
  list(
    label = "Clipping | Ambient",
    group_a = "Ambient + Unclipped",
    group_b = "Ambient + Clipped"
  ),
  list(
    label = "Clipping | Warmed",
    group_a = "Warmed + Unclipped",
    group_b = "Warmed + Clipped"
  )
)

pairwise_rows <- vector("list", length(planned_contrasts))
pair_audit_rows <- vector("list", length(planned_contrasts))
for (contrast_i in seq_along(planned_contrasts)) {
  contrast <- planned_contrasts[[contrast_i]]
  keep_samples <- biological_meta$Treatment %in%
    c(contrast$group_a, contrast$group_b)
  pair_meta <- biological_meta[keep_samples, , drop = FALSE]
  pair_meta$ContrastGroup <- factor(
    as.character(pair_meta$Treatment),
    levels = c(contrast$group_a, contrast$group_b)
  )
  pair_meta <- pair_meta[
    order(pair_meta$Block, pair_meta$ContrastGroup),
    ,
    drop = FALSE
  ]
  pair_distance <- stats::as.dist(
    distance_matrix[rownames(pair_meta), rownames(pair_meta), drop = FALSE]
  )
  pair_permutations <- make_paired_permutations(pair_meta)
  set.seed(primary_seed + contrast_i)
  pair_fit <- vegan::adonis2(
    pair_distance ~ Block + ContrastGroup,
    data = pair_meta,
    permutations = pair_permutations,
    by = "terms",
    parallel = 1
  )
  pair_ss <- pair_fit["ContrastGroup", "SumOfSqs"]
  pair_residual_ss <- pair_fit["Residual", "SumOfSqs"]
  pairwise_rows[[contrast_i]] <- data.frame(
    Contrast = contrast$label,
    GroupA = contrast$group_a,
    GroupB = contrast$group_b,
    Samples = nrow(pair_meta),
    Blocks = nlevels(droplevels(pair_meta$Block)),
    Df = pair_fit["ContrastGroup", "Df"],
    SumOfSquares = pair_ss,
    R2Total = pair_fit["ContrastGroup", "R2"],
    R2Partial = pair_ss / (pair_ss + pair_residual_ss),
    PseudoF = pair_fit["ContrastGroup", "F"],
    PValue = pair_fit["ContrastGroup", "Pr(>F)"],
    Permutations = nrow(pair_permutations),
    MinimumP = 1 / (nrow(pair_permutations) + 1),
    stringsAsFactors = FALSE
  )
  pair_audit_rows[[contrast_i]] <- data.frame(
    Scheme = paste("Exact paired:", contrast$label),
    Seed = NA_integer_,
    Units = nrow(pair_meta),
    Blocks = nlevels(droplevels(pair_meta$Block)),
    MainPlots = NA_integer_,
    PossibleIncludingIdentity = 64L,
    NonIdentityAvailable = 63L,
    PermutationsUsed = nrow(pair_permutations),
    UniqueRows = nrow(unique(pair_permutations)),
    IdentityRows = sum(apply(
      pair_permutations,
      1L,
      function(x) identical(as.integer(x), seq_len(nrow(pair_meta)))
    )),
    LegalRows = sum(vapply(
      seq_len(nrow(pair_permutations)),
      function(i) {
        all(pair_meta$Block[pair_permutations[i, ]] == pair_meta$Block)
      },
      logical(1)
    )),
    MinimumP = 1 / 64,
    MatrixSHA256 = digest::digest(
      pair_permutations,
      algo = "sha256",
      serialize = TRUE
    ),
    stringsAsFactors = FALSE
  )
}
pairwise_permanova <- do.call(rbind, pairwise_rows)
pairwise_permanova$PAdjustedHolm <- stats::p.adjust(
  pairwise_permanova$PValue,
  method = "holm"
)
pairwise_permanova$RejectHolm05 <-
  pairwise_permanova$PAdjustedHolm < 0.05
permutation_audit <- rbind(
  permutation_audit,
  do.call(rbind, pair_audit_rows)
)
add_check("pairwise-test-count", "pairwise", nrow(pairwise_permanova), 4, nrow(pairwise_permanova) == 4L)
add_check("pairwise-permutation-count", "pairwise", unique(pairwise_permanova$Permutations), 63, all(pairwise_permanova$Permutations == 63L))
add_check("pairwise-minimum-p", "pairwise", unique(pairwise_permanova$MinimumP), 1 / 64, all(pairwise_permanova$MinimumP == 1 / 64))
add_check("pairwise-holm-monotonic", "pairwise", all(pairwise_permanova$PAdjustedHolm >= pairwise_permanova$PValue), TRUE, all(pairwise_permanova$PAdjustedHolm >= pairwise_permanova$PValue))
add_check(
  "pairwise-p-grid",
  "pairwise",
  pairwise_permanova$PValue * 64,
  "integers",
  all(abs(pairwise_permanova$PValue * 64 -
    round(pairwise_permanova$PValue * 64)) < 1e-8)
)

dispersion_model <- vegan::betadisper(
  primary_distance,
  group = biological_meta$Treatment,
  type = "median",
  bias.adjust = TRUE
)
set.seed(primary_seed)
dispersion_test <- vegan::permutest(
  dispersion_model,
  permutations = permutation_matrix,
  pairwise = FALSE,
  parallel = 1
)
dispersion_global <- data.frame(
  Grouping = "Treatment (4 levels)",
  Center = "Spatial median",
  BiasAdjustment = TRUE,
  Df = dispersion_test$tab["Groups", "Df"],
  SumOfSquares = dispersion_test$tab["Groups", "Sum Sq"],
  MeanSquare = dispersion_test$tab["Groups", "Mean Sq"],
  FValue = dispersion_test$tab["Groups", "F"],
  PValue = dispersion_test$tab["Groups", "Pr(>F)"],
  Permutations = dispersion_test$tab["Groups", "N.Perm"],
  MinimumP = minimum_primary_p,
  stringsAsFactors = FALSE
)
dispersion_distances <- data.frame(
  BiologicalSampleID = names(dispersion_model$distances),
  Treatment = as.character(biological_meta[
    names(dispersion_model$distances),
    "Treatment"
  ]),
  DistanceToSpatialMedian = as.numeric(dispersion_model$distances),
  stringsAsFactors = FALSE
)
add_check("dispersion-distance-count", "dispersion", nrow(dispersion_distances), 24, nrow(dispersion_distances) == 24L)
add_check("dispersion-distances-finite", "dispersion", sum(is.finite(dispersion_distances$DistanceToSpatialMedian)), 24, all(is.finite(dispersion_distances$DistanceToSpatialMedian)))
add_check("dispersion-permutations", "dispersion", dispersion_global$Permutations, 9999, dispersion_global$Permutations == 9999L)
add_check(
  "dispersion-p-grid",
  "dispersion",
  dispersion_global$PValue * 10000,
  "integer",
  abs(dispersion_global$PValue * 10000 -
    round(dispersion_global$PValue * 10000)) < 1e-8
)

set.seed(primary_seed)
anosim_fit <- vegan::anosim(
  primary_distance,
  grouping = biological_meta$Treatment,
  permutations = permutation_matrix,
  parallel = 1
)
set.seed(primary_seed)
mrpp_fit <- vegan::mrpp(
  primary_distance,
  grouping = biological_meta$Treatment,
  permutations = permutation_matrix,
  weight.type = 1,
  parallel = 1
)
alternative_tests <- data.frame(
  Method = c("PERMANOVA", "ANOSIM", "MRPP"),
  Statistic = c("Pseudo-F", "R", "A"),
  StatisticValue = c(
    permanova_omnibus$PseudoF,
    unname(anosim_fit$statistic),
    unname(mrpp_fit$A)
  ),
  PValue = c(
    permanova_omnibus$PValue,
    anosim_fit$signif,
    mrpp_fit$Pvalue
  ),
  Permutations = primary_nperm,
  EffectSizeInterpretation = c(
    "R2 is reported in permanova-omnibus.tsv",
    "Rank separation; not comparable with R2",
    "Chance-corrected agreement; not comparable with R2"
  ),
  Role = c("Primary", "Secondary sensitivity", "Secondary sensitivity"),
  stringsAsFactors = FALSE
)
add_check("alternative-test-count", "sensitivity", nrow(alternative_tests), 3, nrow(alternative_tests) == 3L)
add_check("alternative-test-permutations", "sensitivity", unique(alternative_tests$Permutations), 9999, all(alternative_tests$Permutations == 9999L))
add_check("alternative-test-p-finite", "sensitivity", sum(is.finite(alternative_tests$PValue)), 3, all(is.finite(alternative_tests$PValue)))

make_free_permutations <- function(n, nset, seed) {
  set.seed(seed)
  out <- t(replicate(nset, sample.int(n), simplify = "matrix"))
  if (nrow(unique(out)) != nset) {
    stop("Free permutation matrix unexpectedly contains duplicates.")
  }
  out
}

run_omnibus <- function(distance, meta, permutations, include_block = TRUE) {
  if (include_block) {
    fit <- vegan::adonis2(
      distance ~ Block + Treatment,
      data = meta,
      permutations = permutations,
      by = "terms",
      parallel = 1
    )
  } else {
    fit <- vegan::adonis2(
      distance ~ Treatment,
      data = meta,
      permutations = permutations,
      by = "terms",
      parallel = 1
    )
  }
  ss <- fit["Treatment", "SumOfSqs"]
  residual_ss <- fit["Residual", "SumOfSqs"]
  c(
    R2Total = fit["Treatment", "R2"],
    R2Partial = ss / (ss + residual_ss),
    PseudoF = fit["Treatment", "F"],
    PValue = fit["Treatment", "Pr(>F)"]
  )
}

unfiltered_relative <- sweep(
  t(aggregate_counts),
  1L,
  rowSums(t(aggregate_counts)),
  "/"
)
unfiltered_distance <- vegan::vegdist(
  unfiltered_relative,
  method = "bray"
)
unfiltered_result <- run_omnibus(
  unfiltered_distance,
  biological_meta,
  permutation_matrix,
  include_block = TRUE
)
free_biological_permutations <- make_free_permutations(
  nrow(biological_meta),
  primary_nperm,
  primary_seed + 1L
)
free_biological_result <- run_omnibus(
  primary_distance,
  biological_meta,
  free_biological_permutations,
  include_block = TRUE
)

library_relative <- sweep(
  t(counts[primary_keep, , drop = FALSE]),
  1L,
  colSums(counts[primary_keep, , drop = FALSE]),
  "/"
)
library_distance <- vegan::vegdist(library_relative, method = "bray")
library_meta <- metadata
free_library_permutations <- make_free_permutations(
  nrow(library_meta),
  primary_nperm,
  primary_seed + 2L
)
technical_result <- run_omnibus(
  library_distance,
  library_meta,
  free_library_permutations,
  include_block = FALSE
)

within_biological_levels <- vapply(
  split(as.character(metadata$Treatment), metadata$BiologicalSampleID),
  function(x) length(unique(x)),
  integer(1)
)
add_check(
  "technical-within-unit-label-immobility",
  "design",
  max(within_biological_levels),
  1,
  all(within_biological_levels == 1L)
)

sensitivity_analysis <- data.frame(
  Branch = c(
    "Primary paper filter + split-plot restriction",
    "No prevalence filter + split-plot restriction",
    "Paper filter + unrestricted biological units",
    "Paper filter + unrestricted technical libraries"
  ),
  Features = c(sum(primary_keep), nrow(counts), sum(primary_keep), sum(primary_keep)),
  Units = c(24L, 24L, 24L, 56L),
  UnitDefinition = c(
    "Biological soil sample",
    "Biological soil sample",
    "Biological soil sample",
    "PCR/tag technical library"
  ),
  PermutationScheme = c(
    "Split-plot restricted",
    "Split-plot restricted",
    "Free; violates paired nested design",
    "Free; pseudoreplication"
  ),
  ValidForPrimaryInference = c(TRUE, TRUE, FALSE, FALSE),
  R2Total = c(
    permanova_omnibus$R2Total,
    unfiltered_result[["R2Total"]],
    free_biological_result[["R2Total"]],
    technical_result[["R2Total"]]
  ),
  R2Partial = c(
    permanova_omnibus$R2Partial,
    unfiltered_result[["R2Partial"]],
    free_biological_result[["R2Partial"]],
    technical_result[["R2Partial"]]
  ),
  PseudoF = c(
    permanova_omnibus$PseudoF,
    unfiltered_result[["PseudoF"]],
    free_biological_result[["PseudoF"]],
    technical_result[["PseudoF"]]
  ),
  PValue = c(
    permanova_omnibus$PValue,
    unfiltered_result[["PValue"]],
    free_biological_result[["PValue"]],
    technical_result[["PValue"]]
  ),
  Permutations = primary_nperm,
  stringsAsFactors = FALSE
)
add_check("sensitivity-branch-count", "sensitivity", nrow(sensitivity_analysis), 4, nrow(sensitivity_analysis) == 4L)
add_check("sensitivity-valid-branches", "sensitivity", sum(sensitivity_analysis$ValidForPrimaryInference), 2, sum(sensitivity_analysis$ValidForPrimaryInference) == 2L)
add_check("technical-pseudorep-marked-invalid", "sensitivity", sensitivity_analysis$ValidForPrimaryInference[[4L]], FALSE, !sensitivity_analysis$ValidForPrimaryInference[[4L]])
add_check("sensitivity-statistics-finite", "sensitivity", sum(is.finite(as.matrix(sensitivity_analysis[c("R2Total", "R2Partial", "PseudoF", "PValue")]))), 16, all(is.finite(as.matrix(sensitivity_analysis[c("R2Total", "R2Partial", "PseudoF", "PValue")]))))

# Use the same full-space principal-coordinate representation that supports
# distances to spatial medians in betadisper. The first two coordinates are
# only a visual projection; all tests above use the complete distance matrix.
site_vectors <- as.data.frame(
  dispersion_model$vectors[, seq_len(2L), drop = FALSE]
)
colnames(site_vectors) <- c("PCoA1", "PCoA2")
site_vectors$BiologicalSampleID <- rownames(site_vectors)
site_vectors$Treatment <- as.character(biological_meta[
  site_vectors$BiologicalSampleID,
  "Treatment"
])
site_vectors$Warming <- as.character(biological_meta[
  site_vectors$BiologicalSampleID,
  "Warming"
])
site_vectors$Clipping <- as.character(biological_meta[
  site_vectors$BiologicalSampleID,
  "Clipping"
])
site_vectors$Block <- as.character(biological_meta[
  site_vectors$BiologicalSampleID,
  "Block"
])

center_vectors <- as.data.frame(
  dispersion_model$centroids[, seq_len(2L), drop = FALSE]
)
colnames(center_vectors) <- c("Center1", "Center2")
center_vectors$Treatment <- rownames(center_vectors)
site_vectors$Center1 <- center_vectors$Center1[
  match(site_vectors$Treatment, center_vectors$Treatment)
]
site_vectors$Center2 <- center_vectors$Center2[
  match(site_vectors$Treatment, center_vectors$Treatment)
]

positive_eigenvalues <- dispersion_model$eig[dispersion_model$eig > 0]
axis_percent <- 100 * dispersion_model$eig[seq_len(2L)] /
  sum(positive_eigenvalues)
pcoa_scores <- data.frame(
  BiologicalSampleID = site_vectors$BiologicalSampleID,
  Block = site_vectors$Block,
  Warming = site_vectors$Warming,
  Clipping = site_vectors$Clipping,
  Treatment = site_vectors$Treatment,
  PCoA1 = site_vectors$PCoA1,
  PCoA2 = site_vectors$PCoA2,
  SpatialMedianPCoA1 = site_vectors$Center1,
  SpatialMedianPCoA2 = site_vectors$Center2,
  DistanceToSpatialMedian = dispersion_model$distances[
    site_vectors$BiologicalSampleID
  ],
  stringsAsFactors = FALSE
)
add_check("pcoa-score-count", "ordination", nrow(pcoa_scores), 24, nrow(pcoa_scores) == 24L)
add_check("pcoa-scores-finite", "ordination", sum(is.finite(as.matrix(pcoa_scores[c("PCoA1", "PCoA2", "SpatialMedianPCoA1", "SpatialMedianPCoA2", "DistanceToSpatialMedian")]))), 120, all(is.finite(as.matrix(pcoa_scores[c("PCoA1", "PCoA2", "SpatialMedianPCoA1", "SpatialMedianPCoA2", "DistanceToSpatialMedian")]))))
add_check("pcoa-axis-percent-positive", "ordination", axis_percent, ">0", all(axis_percent > 0))

treatment_palette <- c(
  "Ambient + Unclipped" = pal_pub[["blue"]],
  "Ambient + Clipped" = pal_pub[["sky"]],
  "Warmed + Unclipped" = pal_pub[["vermillion"]],
  "Warmed + Clipped" = pal_pub[["orange"]]
)
warming_palette <- c(
  "Ambient" = pal_pub[["blue"]],
  "Warmed" = pal_pub[["vermillion"]]
)
clipping_shapes <- c("Unclipped" = 21, "Clipped" = 24)

design_plot_data <- design_audit
design_plot_data$Block <- factor(
  design_plot_data$Block,
  levels = rev(paste("Block", 1:6))
)
design_plot_data$Treatment <- factor(
  design_plot_data$Treatment,
  levels = treatment_levels
)
design_plot_data$Warming <- factor(
  design_plot_data$Warming,
  levels = c("Ambient", "Warmed")
)
design_plot_data$CellLabel <- paste0(
  design_plot_data$TechnicalLibraries,
  ifelse(design_plot_data$TechnicalLibraries == 1L, " library", " libraries"),
  "\n1 analysis unit"
)
design_plot <- ggplot2::ggplot(
  design_plot_data,
  ggplot2::aes(Treatment, Block, fill = Warming)
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 1.1,
    width = 0.96,
    height = 0.92
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = CellLabel),
    colour = "white",
    size = 2.65,
    lineheight = 0.9,
    family = font_family
  ) +
  ggplot2::scale_fill_manual(values = warming_palette, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = c(
    "Ambient + Unclipped" = "Ambient\nUnclipped",
    "Ambient + Clipped" = "Ambient\nClipped",
    "Warmed + Unclipped" = "Warmed\nUnclipped",
    "Warmed + Clipped" = "Warmed\nClipped"
  )) +
  ggplot2::labs(
    title = "Technical libraries collapse to 24 biological units",
    subtitle = paste(
      "One row is one block; whole plots may swap within blocks,",
      "and clipping may swap only within each main plot"
    ),
    x = "Field treatment",
    y = NULL,
    fill = "Whole-plot factor",
    caption = paste(
      "Labels show 1–3 PCR/tag libraries aggregated before inference.",
      "The resulting split-plot space contains 8^6 = 262,144 assignments."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8.3, lineheight = 0.9),
    axis.text.y = ggplot2::element_text(size = 8.3),
    legend.position = "top"
  )

pcoa_plot <- ggplot2::ggplot(
  site_vectors,
  ggplot2::aes(PCoA1, PCoA2)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      xend = Center1,
      yend = Center2,
      colour = Treatment
    ),
    linewidth = 0.42,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = Treatment,
      shape = Clipping
    ),
    colour = "#1A1A1A",
    size = 3.0,
    stroke = 0.55
  ) +
  ggplot2::geom_point(
    data = center_vectors,
    ggplot2::aes(Center1, Center2, colour = Treatment),
    inherit.aes = FALSE,
    shape = 4,
    size = 4.2,
    stroke = 1.2,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = treatment_palette, drop = FALSE) +
  ggplot2::scale_colour_manual(values = treatment_palette, drop = FALSE) +
  ggplot2::scale_shape_manual(values = clipping_shapes, drop = FALSE) +
  ggplot2::labs(
    title = "Location and spread must be read together",
    subtitle = sprintf(
      paste0(
        "Blocked PERMANOVA: R² = %.3f, P = %s; ",
        "PERMDISP: P = %s"
      ),
      permanova_omnibus$R2Total,
      format_p(permanova_omnibus$PValue),
      format_p(dispersion_global$PValue)
    ),
    x = sprintf("PCoA 1 (%.1f%%)", axis_percent[[1L]]),
    y = sprintf("PCoA 2 (%.1f%%)", axis_percent[[2L]]),
    fill = "Treatment",
    shape = "Clipping",
    caption = paste(
      "Crosses are treatment spatial medians; segments show the two-dimensional projection.",
      "PERMDISP uses full-space distances to spatial medians.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    legend.box = "vertical"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 3.0,
        colour = "#1A1A1A"
      )
    ),
    shape = ggplot2::guide_legend(nrow = 1)
  )

pairwise_plot_data <- pairwise_permanova
pairwise_plot_data$Contrast <- factor(
  pairwise_plot_data$Contrast,
  levels = rev(vapply(planned_contrasts, `[[`, character(1), "label"))
)
pairwise_plot_data$Decision <- ifelse(
  pairwise_plot_data$RejectHolm05,
  "Holm-adjusted P < 0.05",
  "Holm-adjusted P >= 0.05"
)
pairwise_plot_data$Label <- sprintf(
  "Partial R² = %.3f  |  Holm P %s",
  pairwise_plot_data$R2Partial,
  format_p(pairwise_plot_data$PAdjustedHolm)
)
pairwise_decision_palette <- c(
  "Holm-adjusted P < 0.05" = pal_pub[["vermillion"]],
  "Holm-adjusted P >= 0.05" = pal_pub[["grey"]]
)
pairwise_x_max <- max(pairwise_plot_data$R2Partial) * 1.75
pairwise_plot <- ggplot2::ggplot(
  pairwise_plot_data,
  ggplot2::aes(R2Partial, Contrast, colour = Decision)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = R2Partial, yend = Contrast),
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    hjust = -0.08,
    colour = "#333333",
    size = 2.65,
    family = font_family
  ) +
  ggplot2::scale_colour_manual(values = pairwise_decision_palette) +
  ggplot2::scale_x_continuous(
    limits = c(0, pairwise_x_max),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Planned simple effects use exact paired permutations",
    subtitle = "Each contrast enumerates 63 non-identity swaps across six blocks",
    x = "Partial variance explained",
    y = NULL,
    colour = "Multiplicity decision",
    caption = paste(
      "The four prespecified P values form one Holm family.",
      "With 64 assignments including the observed order, the minimum P is 0.015625.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "top"
  )

dispersion_plot_data <- dispersion_distances
dispersion_plot_data$Treatment <- factor(
  dispersion_plot_data$Treatment,
  levels = treatment_levels
)
dispersion_plot <- ggplot2::ggplot(
  dispersion_plot_data,
  ggplot2::aes(Treatment, DistanceToSpatialMedian, fill = Treatment)
) +
  ggplot2::geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.35,
    colour = "#333333",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    shape = 21,
    colour = "#1A1A1A",
    size = 2.4,
    stroke = 0.45,
    position = ggplot2::position_jitter(
      width = 0.09,
      height = 0,
      seed = primary_seed
    )
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 4,
    size = 3.4,
    stroke = 1.0,
    colour = "#1A1A1A"
  ) +
  ggplot2::scale_fill_manual(values = treatment_palette, guide = "none") +
  ggplot2::scale_x_discrete(labels = c(
    "Ambient + Unclipped" = "Ambient\nUnclipped",
    "Ambient + Clipped" = "Ambient\nClipped",
    "Warmed + Unclipped" = "Warmed\nUnclipped",
    "Warmed + Clipped" = "Warmed\nClipped"
  )) +
  ggplot2::labs(
    title = "PERMDISP tests spread, not group location",
    subtitle = sprintf(
      "Spatial median, bias-adjusted; F = %.3f, restricted P = %s",
      dispersion_global$FValue,
      format_p(dispersion_global$PValue)
    ),
    x = "Treatment",
    y = "Distance to spatial median",
    caption = paste(
      "Points are 24 biological soil samples; crosses are group means of the distances.",
      "A non-significant result is not proof that dispersions are identical.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8.2, lineheight = 0.9)
  )

figure_specification <- data.frame(
  stem = c(
    "23-exchangeability-design",
    "23-pcoa-centroid-dispersion",
    "23-pairwise-permanova",
    "23-dispersion-by-treatment"
  ),
  width_mm = rep(183, 4L),
  height_mm = c(104, 120, 96, 104),
  raster_ppi = rep(600L, 4L),
  stringsAsFactors = FALSE
)
plot_objects <- list(
  design_plot,
  pcoa_plot,
  pairwise_plot,
  dispersion_plot
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
    text_nodes = lengths(regmatches(
      svg_text,
      gregexpr("<text\\b", svg_text, perl = TRUE)
    )),
    font_family_present = grepl(
      'font-family: "DejaVu Sans"',
      svg_text,
      fixed = TRUE
    ),
    chinese_characters = lengths(regmatches(
      svg_text,
      gregexpr("[\u4e00-\u9fff]", svg_text, perl = TRUE)
    )),
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
      chinese_characters <- svg_info$chinese_characters
      font_status <- ifelse(
        svg_info$font_family_present,
        "declared",
        "missing"
      )
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
      add_check(
        paste0(spec$stem, "-svg-english-only"),
        "vector",
        chinese_characters,
        0,
        chinese_characters == 0L
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
      expected_width_px <- floor(
        spec$width_mm / 25.4 * spec$raster_ppi
      )
      expected_height_px <- floor(
        spec$height_mm / 25.4 * spec$raster_ppi
      )
      add_check(
        paste0(spec$stem, "-", extension, "-pixel-dimensions"),
        "raster",
        c(pixel_width, pixel_height),
        c(expected_width_px, expected_height_px),
        all(c(pixel_width, pixel_height) ==
          c(expected_width_px, expected_height_px))
      )
      add_check(
        paste0(spec$stem, "-", extension, "-effective-ppi"),
        "raster",
        round(c(effective_ppi_x, effective_ppi_y), 3),
        paste0(spec$raster_ppi, " +/- 0.5"),
        all(abs(c(effective_ppi_x, effective_ppi_y) -
          spec$raster_ppi) <= 0.5)
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
      chinese_characters = chinese_characters,
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
add_check(
  "all-svg-english-only",
  "export",
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
write_tsv(design_audit, "design-audit.tsv")
write_tsv(aggregation_audit, "aggregation-audit.tsv")
write_tsv(feature_filter_audit, "feature-filter-audit.tsv")
write_tsv(permutation_audit, "permutation-audit.tsv")
write_tsv(permanova_omnibus, "permanova-omnibus.tsv")
write_tsv(permanova_factorial, "permanova-factorial.tsv")
write_tsv(pairwise_permanova, "pairwise-permanova.tsv")
write_tsv(dispersion_global, "dispersion-global.tsv")
write_tsv(dispersion_distances, "dispersion-distances.tsv")
write_tsv(alternative_tests, "alternative-tests.tsv")
write_tsv(sensitivity_analysis, "sensitivity-analysis.tsv")
write_tsv(pcoa_scores, "pcoa-scores.tsv")
write_tsv(figure_specification, "figure-specification.tsv")
write_tsv(format_audit, "format-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

summary_payload <- list(
  article = 23L,
  title = "PERMANOVA, restricted permutations, and multivariate dispersion",
  dataset = list(
    source = "phyloseq 1.48.0 soilrep",
    features = nrow(counts),
    technical_libraries = ncol(counts),
    biological_samples = ncol(aggregate_counts),
    blocks = nlevels(biological_meta$Block),
    main_plots = length(unique(biological_meta$MainPlot)),
    treatments = nlevels(biological_meta$Treatment),
    reads = sum(counts),
    minimum_biological_depth = min(colSums(aggregate_counts)),
    median_biological_depth = unname(stats::median(colSums(aggregate_counts))),
    maximum_biological_depth = max(colSums(aggregate_counts))
  ),
  standardization = list(
    primary = paste(
      "Relative abundance Bray-Curtis after retaining features detected",
      "in >=4 tagged libraries in any treatment"
    ),
    features_retained = sum(primary_keep),
    distance = "Bray-Curtis"
  ),
  permutations = list(
    seed = primary_seed,
    possible = splitplot$possible,
    primary = primary_nperm,
    primary_unique = nrow(unique(permutation_matrix)),
    primary_identity_rows = identity_rows,
    primary_legal_rows = sum(legal_rows),
    primary_minimum_p = minimum_primary_p,
    pairwise_each = 63L,
    pairwise_minimum_p = 1 / 64
  ),
  analysis = list(
    conditioning_block_r2 = permanova_omnibus$ConditioningBlockR2,
    omnibus_r2_total = permanova_omnibus$R2Total,
    omnibus_r2_partial = permanova_omnibus$R2Partial,
    omnibus_p = permanova_omnibus$PValue,
    factorial_terms = nrow(permanova_factorial),
    factorial_rejections_holm = sum(permanova_factorial$RejectHolm05),
    pairwise_tests = nrow(pairwise_permanova),
    pairwise_rejections_holm = sum(pairwise_permanova$RejectHolm05),
    dispersion_p = dispersion_global$PValue,
    anosim_r = unname(anosim_fit$statistic),
    anosim_p = anosim_fit$signif,
    mrpp_a = unname(mrpp_fit$A),
    mrpp_p = mrpp_fit$Pvalue,
    sensitivity_branches = nrow(sensitivity_analysis)
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
    dimension_checks_passed = sum(format_audit$dimension_status == "PASS")
  ),
  versions = as.list(observed_versions),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "permanova-dispersion-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 10
)

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"))
log_lines <- c(
  "Article 23 PERMANOVA and dispersion validation",
  paste0("R version: ", r_version),
  paste0(
    "Dataset: ", nrow(counts), " features x ", ncol(counts),
    " technical libraries; ", ncol(aggregate_counts),
    " biological samples"
  ),
  paste0("Reads conserved after aggregation: ", sum(aggregate_counts)),
  paste0("Primary features retained: ", sum(primary_keep)),
  paste0(
    "Primary permutation space / used: ", splitplot$possible,
    " / ", primary_nperm
  ),
  sprintf(
    "Omnibus PERMANOVA: R2=%.6f; partial R2=%.6f; P=%.6f",
    permanova_omnibus$R2Total,
    permanova_omnibus$R2Partial,
    permanova_omnibus$PValue
  ),
  sprintf(
    "PERMDISP: F=%.6f; P=%.6f",
    dispersion_global$FValue,
    dispersion_global$PValue
  ),
  paste0(
    "Pairwise Holm rejections: ",
    sum(pairwise_permanova$RejectHolm05),
    "/",
    nrow(pairwise_permanova)
  ),
  sprintf(
    "ANOSIM R/P: %.6f / %.6f; MRPP A/P: %.6f / %.6f",
    unname(anosim_fit$statistic),
    anosim_fit$signif,
    unname(mrpp_fit$A),
    mrpp_fit$Pvalue
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
writeLines(
  sanitize_text(log_lines),
  file.path(output_dir, "validation.log")
)

if (checks_failed > 0L) {
  failed <- validation_checks$check_id[
    validation_checks$status == "FAIL"
  ]
  stop(
    "Article 23 validation failed: ",
    paste(failed, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Article 23 validation passed: ",
  checks_passed,
  "/",
  checks_total,
  " checks."
)
