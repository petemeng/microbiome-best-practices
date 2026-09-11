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
  "ape", "Biostrings", "digest", "dplyr", "ggplot2", "jsonlite",
  "phyloseq", "readr", "scales", "tibble", "tidyr"
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

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
  x <- gsub(input_dir, "<INPUT_DIR>", x, fixed = TRUE)
  x <- gsub(output_dir, "<OUTPUT_DIR>", x, fixed = TRUE)
  x <- gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
  x
}

write_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(output_dir, filename), na = "")
}

save_publication_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 350,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 350,
    compression = "lzw",
    bg = "white"
  )
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "#4D4D4D"),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", color = "#B3B3B3"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
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

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  taxonomy_confidence = file.path(input_dir, "taxonomy-confidence.tsv"),
  representative_sequences = file.path(input_dir, "representative-sequences.fasta.gz"),
  rooted_insertion_tree = file.path(input_dir, "rooted-insertion-tree.nwk.gz"),
  source_summary = file.path(input_dir, "source-summary.json")
)
expected_sha256 <- c(
  otutab = "6c43bb8f929f935b60e42ed32737632d3be18ddd491f5298ba3a6abd6e3f18a9",
  taxonomy = "8ef235ea19923307eef4d8e1905c2f34b78cf25cf1a83f08aa00353f2c5c0d38",
  metadata = "6628dd51051c53884e82c540ceaff8cc34f69900a78d372239a2ed41ff5b7263",
  taxonomy_confidence = "26fd4c66301f97c6feac8f3302ab19f269e9429c9ce65dea5c66ae277eca55b7",
  representative_sequences = "a19e083cfb0ab757d557569bdf76aa3de78576ae35033a7e69a647976a370bc2",
  rooted_insertion_tree = "7e541e2ed9e3f81d9b6c6bf98999335e7f684a37ac5dd5f2b56e61b2f4c4c47c",
  source_summary = "358710a72413e1bc1ff0e13cb1a3722a848583d37321e1c2442523f09d3226fd"
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing Article 17 input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
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
  phyloseq = "1.48.0",
  ape = "5.8",
  Biostrings = "2.72.1",
  ggplot2 = "3.5.2",
  readr = "2.1.5",
  digest = "0.6.36",
  jsonlite = "1.8.8"
)
observed_versions <- vapply(
  names(expected_versions),
  function(package) as.character(utils::packageVersion(package)),
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

source_summary <- jsonlite::read_json(
  input_paths[["source_summary"]],
  simplifyVector = TRUE
)
otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])
taxonomy_confidence <- readr::read_tsv(
  input_paths[["taxonomy_confidence"]],
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal",
  na = character()
)

otu_matrix <- as.matrix(otutab)
storage.mode(otu_matrix) <- "numeric"
tax_matrix <- as.matrix(taxonomy)
storage.mode(tax_matrix) <- "character"
feature_ids <- rownames(otu_matrix)
sample_ids <- colnames(otu_matrix)

sequences <- Biostrings::readDNAStringSet(
  input_paths[["representative_sequences"]]
)
full_tree <- ape::read.tree(
  gzfile(input_paths[["rooted_insertion_tree"]])
)

expected_sample_reads <- c("1" = 1144, "1a" = 1251, "2" = 1503, "2a" = 1315)
expected_sample_asvs <- c("1" = 82, "1a" = 133, "2" = 106, "2a" = 138)
sample_reads <- colSums(otu_matrix)
sample_asvs <- colSums(otu_matrix > 0)

add_check("feature-count", "shape", nrow(otu_matrix), 366, nrow(otu_matrix) == 366)
add_check("sample-count", "shape", ncol(otu_matrix), 4, ncol(otu_matrix) == 4)
add_check("read-count", "shape", sum(otu_matrix), 5213, sum(otu_matrix) == 5213)
add_check("taxonomy-ranks", "shape", ncol(tax_matrix), 7, ncol(tax_matrix) == 7)
add_check(
  "taxonomy-rank-names",
  "shape",
  colnames(tax_matrix),
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  identical(
    colnames(tax_matrix),
    c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
)
add_check("feature-ids-unique", "ids", anyDuplicated(feature_ids), 0, !anyDuplicated(feature_ids))
add_check("sample-ids-unique", "ids", anyDuplicated(sample_ids), 0, !anyDuplicated(sample_ids))
add_check("counts-finite", "counts", all(is.finite(otu_matrix)), TRUE, all(is.finite(otu_matrix)))
add_check("counts-nonnegative", "counts", min(otu_matrix), ">=0", min(otu_matrix) >= 0)
add_check(
  "counts-integer",
  "counts",
  max(abs(otu_matrix - round(otu_matrix))),
  0,
  all(otu_matrix == round(otu_matrix))
)
add_check(
  "sample-library-sizes",
  "counts",
  sample_reads[names(expected_sample_reads)],
  expected_sample_reads,
  identical(as.numeric(sample_reads[names(expected_sample_reads)]), as.numeric(expected_sample_reads))
)
add_check(
  "sample-observed-asvs",
  "counts",
  sample_asvs[names(expected_sample_asvs)],
  expected_sample_asvs,
  identical(as.numeric(sample_asvs[names(expected_sample_asvs)]), as.numeric(expected_sample_asvs))
)
add_check(
  "taxonomy-feature-ids",
  "ids",
  length(intersect(feature_ids, rownames(tax_matrix))),
  366,
  setequal(feature_ids, rownames(tax_matrix))
)
add_check(
  "metadata-sample-ids",
  "ids",
  length(intersect(sample_ids, rownames(metadata))),
  4,
  setequal(sample_ids, rownames(metadata))
)
add_check(
  "confidence-feature-ids",
  "ids",
  length(intersect(feature_ids, taxonomy_confidence[["FeatureID"]])),
  366,
  setequal(feature_ids, taxonomy_confidence[["FeatureID"]])
)
add_check(
  "sequence-feature-ids",
  "ids",
  length(intersect(feature_ids, names(sequences))),
  366,
  setequal(feature_ids, names(sequences))
)

sequence_lengths <- Biostrings::width(sequences)
sequence_strings <- as.character(sequences)
valid_dna <- !grepl("[^ACGTN]", sequence_strings)
add_check("sequence-count", "sequences", length(sequences), 366, length(sequences) == 366)
add_check("sequence-min-length", "sequences", min(sequence_lengths), 252, min(sequence_lengths) == 252)
add_check("sequence-mode-length", "sequences", names(which.max(table(sequence_lengths))), 253, names(which.max(table(sequence_lengths))) == "253")
add_check("sequence-max-length", "sequences", max(sequence_lengths), 303, max(sequence_lengths) == 303)
add_check("sequence-alphabet", "sequences", sum(valid_dna), 366, all(valid_dna))

full_tree_tips <- full_tree$tip.label
tree_query_tips <- intersect(full_tree_tips, feature_ids)
tree_reference_tips <- setdiff(full_tree_tips, feature_ids)
add_check("full-tree-tip-count", "tree", ape::Ntip(full_tree), 203818, ape::Ntip(full_tree) == 203818)
add_check("full-tree-query-tips", "tree", length(tree_query_tips), 366, setequal(tree_query_tips, feature_ids))
add_check("full-tree-reference-tips", "tree", length(tree_reference_tips), 203452, length(tree_reference_tips) == 203452)
add_check("full-tree-rooted", "tree", ape::is.rooted(full_tree), TRUE, ape::is.rooted(full_tree))
add_check("full-tree-missing-branches", "tree", sum(is.na(full_tree$edge.length)), 0, !anyNA(full_tree$edge.length))
add_check("full-tree-negative-branches", "tree", sum(full_tree$edge.length < 0), 0, all(full_tree$edge.length >= 0))

query_tree <- ape::keep.tip(full_tree, feature_ids)
add_check("query-tree-tip-count", "tree", ape::Ntip(query_tree), 366, ape::Ntip(query_tree) == 366)
add_check("query-tree-tip-ids", "tree", length(intersect(query_tree$tip.label, feature_ids)), 366, setequal(query_tree$tip.label, feature_ids))
add_check("query-tree-rooted", "tree", ape::is.rooted(query_tree), TRUE, ape::is.rooted(query_tree))
add_check("query-tree-missing-branches", "tree", sum(is.na(query_tree$edge.length)), 0, !anyNA(query_tree$edge.length))
add_check("query-tree-negative-branches", "tree", sum(query_tree$edge.length < 0), 0, all(query_tree$edge.length >= 0))

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
sequences <- sequences[feature_ids]

phyloseq_object <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(taxonomy)),
  phyloseq::sample_data(metadata),
  phyloseq::phy_tree(query_tree),
  sequences
)

component_names <- slotNames(phyloseq_object)[
  !vapply(
    slotNames(phyloseq_object),
    function(slot_name) is.null(methods::slot(phyloseq_object, slot_name)),
    logical(1)
  )
]
object_reads <- sum(phyloseq::sample_sums(phyloseq_object))
add_check("object-class", "object", class(phyloseq_object)[[1L]], "phyloseq", inherits(phyloseq_object, "phyloseq"))
add_check("object-components", "object", length(component_names), 5, length(component_names) == 5)
add_check("object-component-names", "object", sort(component_names), sort(c("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")), setequal(component_names, c("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")))
add_check("object-taxa", "object", phyloseq::ntaxa(phyloseq_object), 366, phyloseq::ntaxa(phyloseq_object) == 366)
add_check("object-samples", "object", phyloseq::nsamples(phyloseq_object), 4, phyloseq::nsamples(phyloseq_object) == 4)
add_check("object-reads", "object", object_reads, 5213, object_reads == 5213)
add_check("object-taxa-orientation", "object", phyloseq::taxa_are_rows(phyloseq_object), TRUE, phyloseq::taxa_are_rows(phyloseq_object))
add_check("object-taxa-ids", "object", length(intersect(phyloseq::taxa_names(phyloseq_object), feature_ids)), 366, setequal(phyloseq::taxa_names(phyloseq_object), feature_ids))
add_check("object-sample-ids", "object", length(intersect(phyloseq::sample_names(phyloseq_object), sample_ids)), 4, setequal(phyloseq::sample_names(phyloseq_object), sample_ids))
add_check("object-tree-tips", "object", ape::Ntip(phyloseq::phy_tree(phyloseq_object)), 366, setequal(phyloseq::phy_tree(phyloseq_object)$tip.label, feature_ids))
add_check("object-refseq-ids", "object", length(phyloseq::refseq(phyloseq_object)), 366, setequal(names(phyloseq::refseq(phyloseq_object)), feature_ids))
add_check("object-sample-sums", "object", phyloseq::sample_sums(phyloseq_object)[sample_ids], sample_reads[sample_ids], identical(as.numeric(phyloseq::sample_sums(phyloseq_object)[sample_ids]), as.numeric(sample_reads[sample_ids])))

rds_path <- file.path(output_dir, "atacama-v4-phyloseq.rds")
saveRDS(phyloseq_object, rds_path, version = 3, compress = "xz")
roundtrip_object <- readRDS(rds_path)
roundtrip_audit <- data.frame(
  metric = c("class", "components", "features", "samples", "reads", "taxa_ids", "sample_ids"),
  before = c(
    class(phyloseq_object)[[1L]],
    length(component_names),
    phyloseq::ntaxa(phyloseq_object),
    phyloseq::nsamples(phyloseq_object),
    object_reads,
    length(phyloseq::taxa_names(phyloseq_object)),
    length(phyloseq::sample_names(phyloseq_object))
  ),
  after = c(
    class(roundtrip_object)[[1L]],
    sum(!vapply(slotNames(roundtrip_object), function(slot_name) is.null(methods::slot(roundtrip_object, slot_name)), logical(1))),
    phyloseq::ntaxa(roundtrip_object),
    phyloseq::nsamples(roundtrip_object),
    sum(phyloseq::sample_sums(roundtrip_object)),
    length(phyloseq::taxa_names(roundtrip_object)),
    length(phyloseq::sample_names(roundtrip_object))
  ),
  status = "PASS",
  check.names = FALSE
)
roundtrip_audit$status[roundtrip_audit$before != roundtrip_audit$after] <- "FAIL"
add_check("rds-roundtrip-class", "roundtrip", class(roundtrip_object)[[1L]], "phyloseq", inherits(roundtrip_object, "phyloseq"))
add_check("rds-roundtrip-components", "roundtrip", roundtrip_audit$after[roundtrip_audit$metric == "components"], 5, roundtrip_audit$after[roundtrip_audit$metric == "components"] == "5")
add_check("rds-roundtrip-taxa", "roundtrip", phyloseq::ntaxa(roundtrip_object), 366, phyloseq::ntaxa(roundtrip_object) == 366)
add_check("rds-roundtrip-samples", "roundtrip", phyloseq::nsamples(roundtrip_object), 4, phyloseq::nsamples(roundtrip_object) == 4)
add_check("rds-roundtrip-reads", "roundtrip", sum(phyloseq::sample_sums(roundtrip_object)), 5213, sum(phyloseq::sample_sums(roundtrip_object)) == 5213)
add_check("rds-roundtrip-taxa-ids", "roundtrip", length(intersect(phyloseq::taxa_names(roundtrip_object), feature_ids)), 366, setequal(phyloseq::taxa_names(roundtrip_object), feature_ids))
add_check("rds-roundtrip-sample-ids", "roundtrip", length(intersect(phyloseq::sample_names(roundtrip_object), sample_ids)), 4, setequal(phyloseq::sample_names(roundtrip_object), sample_ids))

rank_names <- colnames(taxonomy)
taxonomy_coverage <- do.call(
  rbind,
  lapply(rank_names, function(rank_name) {
    assigned <- !is.na(taxonomy[[rank_name]]) &
      nzchar(taxonomy[[rank_name]]) &
      taxonomy[[rank_name]] != "Unassigned"
    data.frame(
      rank = rank_name,
      assigned_features = sum(assigned),
      total_features = length(assigned),
      feature_percent = 100 * mean(assigned),
      assigned_reads = sum(rowSums(otu_matrix)[assigned]),
      total_reads = sum(otu_matrix),
      read_percent = 100 * sum(rowSums(otu_matrix)[assigned]) / sum(otu_matrix),
      check.names = FALSE
    )
  })
)
genus_row <- taxonomy_coverage[taxonomy_coverage$rank == "Genus", , drop = FALSE]
species_row <- taxonomy_coverage[taxonomy_coverage$rank == "Species", , drop = FALSE]
add_check("genus-assigned-features", "taxonomy", genus_row$assigned_features, 318, genus_row$assigned_features == 318)
add_check("genus-assigned-reads", "taxonomy", genus_row$assigned_reads, 4217, genus_row$assigned_reads == 4217)
add_check("species-assigned-features", "taxonomy", species_row$assigned_features, 0, species_row$assigned_features == 0)
add_check("species-assigned-reads", "taxonomy", species_row$assigned_reads, 0, species_row$assigned_reads == 0)

sample_library_audit <- data.frame(
  sample_id = sample_ids,
  library_label = metadata[sample_ids, "LibraryLabel"],
  filename_set = metadata[sample_ids, "FilenameSet"],
  reads = as.numeric(sample_reads[sample_ids]),
  observed_asvs = as.numeric(sample_asvs[sample_ids]),
  library_percent = 100 * as.numeric(sample_reads[sample_ids]) / sum(sample_reads),
  status = ifelse(
    sample_reads[sample_ids] == expected_sample_reads[sample_ids] &
      sample_asvs[sample_ids] == expected_sample_asvs[sample_ids],
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)

id_reconciliation <- data.frame(
  component = c(
    "Taxonomy", "Representative sequences", "Full-tree query tips",
    "Pruned tree", "phyloseq taxa", "Metadata", "phyloseq samples"
  ),
  id_domain = c(rep("Feature IDs", 5), rep("Sample IDs", 2)),
  expected_ids = c(rep(366, 5), rep(4, 2)),
  matched_ids = c(
    length(intersect(feature_ids, rownames(taxonomy))),
    length(intersect(feature_ids, names(sequences))),
    length(tree_query_tips),
    length(intersect(feature_ids, query_tree$tip.label)),
    length(intersect(feature_ids, phyloseq::taxa_names(phyloseq_object))),
    length(intersect(sample_ids, rownames(metadata))),
    length(intersect(sample_ids, phyloseq::sample_names(phyloseq_object)))
  ),
  missing_ids = c(
    length(setdiff(feature_ids, rownames(taxonomy))),
    length(setdiff(feature_ids, names(sequences))),
    length(setdiff(feature_ids, full_tree_tips)),
    length(setdiff(feature_ids, query_tree$tip.label)),
    length(setdiff(feature_ids, phyloseq::taxa_names(phyloseq_object))),
    length(setdiff(sample_ids, rownames(metadata))),
    length(setdiff(sample_ids, phyloseq::sample_names(phyloseq_object)))
  ),
  unexpected_ids = c(
    length(setdiff(rownames(taxonomy), feature_ids)),
    length(setdiff(names(sequences), feature_ids)),
    0,
    length(setdiff(query_tree$tip.label, feature_ids)),
    length(setdiff(phyloseq::taxa_names(phyloseq_object), feature_ids)),
    length(setdiff(rownames(metadata), sample_ids)),
    length(setdiff(phyloseq::sample_names(phyloseq_object), sample_ids))
  ),
  allowed_reference_only_ids = c(0, 0, length(tree_reference_tips), 0, 0, 0, 0),
  check.names = FALSE
)
id_reconciliation$match_percent <- 100 * id_reconciliation$matched_ids / id_reconciliation$expected_ids
id_reconciliation$status <- ifelse(
  id_reconciliation$matched_ids == id_reconciliation$expected_ids &
    id_reconciliation$missing_ids == 0 &
    id_reconciliation$unexpected_ids == 0,
  "PASS",
  "FAIL"
)

component_audit <- data.frame(
  component = c("otu_table", "tax_table", "sample_data", "phy_tree", "refseq"),
  class = c(
    class(phyloseq::otu_table(phyloseq_object))[[1L]],
    class(phyloseq::tax_table(phyloseq_object))[[1L]],
    class(phyloseq::sample_data(phyloseq_object))[[1L]],
    class(phyloseq::phy_tree(phyloseq_object))[[1L]],
    class(phyloseq::refseq(phyloseq_object))[[1L]]
  ),
  primary_dimension = c(366, 366, 4, 366, 366),
  secondary_dimension = c(4, 7, ncol(metadata), ape::Nnode(query_tree), NA_integer_),
  orientation_or_role = c(
    "features x samples", "features x ranks", "samples x variables",
    "ASV-pruned rooted tree", "ASV reference sequences"
  ),
  status = "PASS",
  check.names = FALSE
)

tree_audit <- data.frame(
  tree = c("Complete SEPP insertion tree", "ASV-pruned object tree"),
  tips = c(ape::Ntip(full_tree), ape::Ntip(query_tree)),
  query_tips = c(length(tree_query_tips), ape::Ntip(query_tree)),
  reference_only_tips = c(length(tree_reference_tips), 0),
  internal_nodes = c(ape::Nnode(full_tree), ape::Nnode(query_tree)),
  rooted = c(ape::is.rooted(full_tree), ape::is.rooted(query_tree)),
  missing_branch_lengths = c(sum(is.na(full_tree$edge.length)), sum(is.na(query_tree$edge.length))),
  negative_branch_lengths = c(sum(full_tree$edge.length < 0), sum(query_tree$edge.length < 0)),
  status = "PASS",
  check.names = FALSE
)

sequence_audit <- data.frame(
  metric = c("records", "minimum_length", "median_length", "modal_length", "maximum_length", "invalid_alphabet_records"),
  observed = c(
    length(sequences),
    min(sequence_lengths),
    stats::median(sequence_lengths),
    as.numeric(names(which.max(table(sequence_lengths)))),
    max(sequence_lengths),
    sum(!valid_dna)
  ),
  expected = c(366, 252, 253, 253, 303, 0),
  status = "PASS",
  check.names = FALSE
)
sequence_audit$status[sequence_audit$observed != sequence_audit$expected] <- "FAIL"

workflow_nodes <- data.frame(
  id = c("counts", "taxonomy", "metadata", "sequences", "tree", "gate", "object"),
  x = c(1.0, 1.0, 1.0, 1.0, 1.0, 3.4, 5.8),
  y = c(4.2, 3.2, 2.2, 1.2, 0.2, 2.2, 2.2),
  width = c(rep(1.7, 5), 1.8, 2.3),
  height = c(rep(0.65, 5), 1.05, 2.85),
  label = c(
    "otutab\n366 x 4",
    "taxonomy\n366 x 7",
    "metadata\n4 x 5",
    "refseq\n366 ASVs",
    "full rooted tree\n203,818 tips",
    "ID + orientation\ngate",
    "phyloseq\notu_table\ntax_table\nsample_data\nphy_tree\nrefseq"
  ),
  group = c(rep("Input", 5), "Gate", "Object"),
  stringsAsFactors = FALSE
)
workflow_edges <- data.frame(
  x = c(rep(1.85, 5), 4.3),
  y = c(4.2, 3.2, 2.2, 1.2, 0.2, 2.2),
  xend = c(rep(2.5, 5), 4.65),
  yend = c(rep(2.2, 5), 2.2)
)
workflow_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = workflow_edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    color = "#7F7F7F",
    linewidth = 0.6,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = workflow_nodes,
    ggplot2::aes(
      xmin = x - width / 2,
      xmax = x + width / 2,
      ymin = y - height / 2,
      ymax = y + height / 2,
      fill = group
    ),
    color = "#333333",
    linewidth = 0.45
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = x, y = y, label = label),
    size = 3.2,
    lineheight = 0.93
  ) +
  ggplot2::scale_fill_manual(
    values = c(Input = "#D7EBF7", Gate = "#FCE5CD", Object = "#D9EAD3")
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 7.15), ylim = c(-0.45, 4.75), clip = "off") +
  ggplot2::labs(
    title = "A phyloseq object is an audited intersection",
    subtitle = "Every component is keyed before the five-slot object is serialized",
    fill = NULL
  ) +
  ggplot2::theme_void(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "#4D4D4D"),
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )
save_publication_plot(workflow_plot, "17-phyloseq-object-contract", 8.0, 5.0)

id_plot <- ggplot2::ggplot(
  id_reconciliation,
  ggplot2::aes(x = match_percent, y = stats::reorder(component, match_percent), color = id_domain)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = match_percent, yend = stats::reorder(component, match_percent)),
    color = "#D9D9D9",
    linewidth = 1.7
  ) +
  ggplot2::geom_point(size = 3.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(matched_ids, " / ", expected_ids)),
    hjust = 1.12,
    color = "#222222",
    size = 3.1,
    fontface = "bold"
  ) +
  ggplot2::scale_color_manual(values = c("Feature IDs" = "#0072B2", "Sample IDs" = "#D55E00")) +
  ggplot2::scale_x_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Identifier reconciliation before analysis",
    subtitle = "No feature or sample is silently lost; 203,452 reference-only tips are expected before pruning",
    x = "Matched identifiers",
    y = NULL,
    color = "ID domain"
  ) +
  theme_pub(10.5)
save_publication_plot(id_plot, "17-id-reconciliation-audit", 8.0, 5.1)

library_long <- tidyr::pivot_longer(
  sample_library_audit,
  cols = c("reads", "observed_asvs"),
  names_to = "metric",
  values_to = "value"
)
library_long$metric <- factor(
  library_long$metric,
  levels = c("reads", "observed_asvs"),
  labels = c("Reads", "Observed ASVs")
)
library_long$display_label <- factor(
  gsub("_", " ", library_long$library_label, fixed = TRUE),
  levels = gsub("_", " ", sample_library_audit$library_label, fixed = TRUE)
)
library_plot <- ggplot2::ggplot(
  library_long,
  ggplot2::aes(x = display_label, y = value, fill = filename_set)
) +
  ggplot2::geom_col(width = 0.72, color = "#333333", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(value)),
    vjust = -0.35,
    size = 3.1
  ) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
  ggplot2::scale_fill_manual(values = c(S103 = "#0072B2", S115 = "#E69F00")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
  ggplot2::labs(
    title = "The imported object preserves each library",
    subtitle = "Filename sets are technical labels, not biological treatment groups",
    x = NULL,
    y = NULL,
    fill = "Filename set"
  ) +
  theme_pub(10.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
save_publication_plot(library_plot, "17-library-feature-audit", 7.2, 6.3)

coverage_long <- tidyr::pivot_longer(
  taxonomy_coverage,
  cols = c("feature_percent", "read_percent"),
  names_to = "coverage_type",
  values_to = "percent"
)
coverage_long$coverage_type <- factor(
  coverage_long$coverage_type,
  levels = c("feature_percent", "read_percent"),
  labels = c("Feature coverage", "Read-weighted coverage")
)
coverage_long$rank <- factor(coverage_long$rank, levels = rank_names)
coverage_plot <- ggplot2::ggplot(
  coverage_long,
  ggplot2::aes(x = rank, y = percent, fill = coverage_type)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.66) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", percent)),
    position = ggplot2::position_dodge(width = 0.72),
    vjust = ifelse(coverage_long$percent > 12, 1.3, -0.35),
    color = ifelse(coverage_long$percent > 12, "white", "black"),
    size = 2.8
  ) +
  ggplot2::scale_fill_manual(values = c("Feature coverage" = "#0072B2", "Read-weighted coverage" = "#009E73")) +
  ggplot2::scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Taxonomic depth is part of the object contract",
    subtitle = "SILVA 138.2 V4 classification at confidence 0.70; species calls are intentionally absent",
    x = "Taxonomic rank",
    y = "Assigned fraction",
    fill = NULL
  ) +
  theme_pub(10.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_publication_plot(coverage_plot, "17-taxonomy-coverage-audit", 8.0, 5.0)

figure_stems <- c(
  "17-phyloseq-object-contract",
  "17-id-reconciliation-audit",
  "17-library-feature-audit",
  "17-taxonomy-coverage-audit"
)
for (stem in figure_stems) {
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    add_check(
      paste0("figure-", stem, "-", extension),
      "figures",
      file.exists(path),
      TRUE,
      file.exists(path) && file.info(path)$size > 0
    )
  }
}

validation_checks <- do.call(rbind, check_rows)
checks_passed <- sum(validation_checks$status == "PASS")
checks_total <- nrow(validation_checks)

write_tsv(input_audit, "input-audit.tsv")
write_tsv(package_audit, "package-audit.tsv")
write_tsv(id_reconciliation, "id-reconciliation.tsv")
write_tsv(component_audit, "component-audit.tsv")
write_tsv(sample_library_audit, "sample-library-audit.tsv")
write_tsv(taxonomy_coverage, "taxonomy-coverage.tsv")
write_tsv(tree_audit, "tree-audit.tsv")
write_tsv(sequence_audit, "sequence-audit.tsv")
write_tsv(roundtrip_audit, "rds-roundtrip-audit.tsv")
write_tsv(validation_checks, "validation-checks.tsv")

session_lines <- sanitize_text(capture.output(utils::sessionInfo()))
writeLines(session_lines, file.path(output_dir, "r-session-info.txt"), useBytes = TRUE)

summary <- list(
  schema_version = 1,
  bundle_id = source_summary$bundle_id,
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_total - checks_passed,
  environment = list(
    r = paste(R.version$major, R.version$minor, sep = "."),
    packages = as.list(observed_versions)
  ),
  input = list(
    features = nrow(otu_matrix),
    samples = ncol(otu_matrix),
    reads = as.integer(sum(otu_matrix)),
    taxonomy_ranks = ncol(taxonomy),
    representative_sequences = length(sequences)
  ),
  object = list(
    class = class(phyloseq_object)[[1L]],
    components = length(component_names),
    component_names = unname(component_names),
    taxa_are_rows = phyloseq::taxa_are_rows(phyloseq_object),
    features = phyloseq::ntaxa(phyloseq_object),
    samples = phyloseq::nsamples(phyloseq_object),
    reads = as.integer(object_reads)
  ),
  tree = list(
    full_tips = ape::Ntip(full_tree),
    full_reference_tips = length(tree_reference_tips),
    full_query_tips = length(tree_query_tips),
    full_rooted = ape::is.rooted(full_tree),
    object_tips = ape::Ntip(query_tree),
    object_rooted = ape::is.rooted(query_tree),
    missing_branch_lengths = sum(is.na(query_tree$edge.length)),
    negative_branch_lengths = sum(query_tree$edge.length < 0)
  ),
  taxonomy = list(
    genus_assigned_features = genus_row$assigned_features[[1L]],
    genus_assigned_reads = genus_row$assigned_reads[[1L]],
    species_assigned_features = species_row$assigned_features[[1L]],
    species_assigned_reads = species_row$assigned_reads[[1L]]
  ),
  sequences = list(
    minimum_length = min(sequence_lengths),
    modal_length = as.integer(names(which.max(table(sequence_lengths)))),
    maximum_length = max(sequence_lengths)
  ),
  rds = list(
    file = basename(rds_path),
    bytes = unname(file.info(rds_path)$size),
    sha256 = sha256_file(rds_path),
    roundtrip_passed = all(roundtrip_audit$status == "PASS")
  )
)
jsonlite::write_json(
  summary,
  file.path(output_dir, "phyloseq-import-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

log_lines <- c(
  "Article 17 phyloseq import validation",
  paste0("Input bundle: ", source_summary$bundle_id),
  paste0("Input: ", nrow(otu_matrix), " features x ", ncol(otu_matrix), " samples; ", sum(otu_matrix), " reads"),
  paste0("Complete tree: ", ape::Ntip(full_tree), " tips; object tree: ", ape::Ntip(query_tree), " tips"),
  paste0("Object components: ", paste(component_names, collapse = ", ")),
  paste0("Genus coverage: ", genus_row$assigned_features, "/366 features; ", genus_row$assigned_reads, "/5213 reads"),
  paste0("RDS round-trip: ", ifelse(all(roundtrip_audit$status == "PASS"), "PASS", "FAIL")),
  paste0("Checks: ", checks_passed, "/", checks_total, " PASS")
)
writeLines(sanitize_text(log_lines), file.path(output_dir, "validation.log"), useBytes = TRUE)

cat(paste(log_lines, collapse = "\n"), "\n", sep = "")
if (checks_passed != checks_total) {
  failed <- validation_checks[validation_checks$status != "PASS", , drop = FALSE]
  print(failed, row.names = FALSE)
  stop("Article 17 validation failed.", call. = FALSE)
}
