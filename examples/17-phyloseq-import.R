# 导出到 R：构建 phyloseq 对象、数据结构详解
# Run sequentially in a new working directory.
# Required packages: Biostrings, ape, digest, dplyr, ggplot2, phyloseq, readr, tibble, tidyr.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/validate_phyloseq_import.R",
  "data/small/phyloseq-import/metadata.tsv",
  "data/small/phyloseq-import/otutab.tsv",
  "data/small/phyloseq-import/representative-sequences.fasta.gz",
  "data/small/phyloseq-import/rooted-insertion-tree.nwk.gz",
  "data/small/phyloseq-import/source-summary.json",
  "data/small/phyloseq-import/taxonomy-confidence.tsv",
  "data/small/phyloseq-import/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

set.seed(20260721)

required_packages <- c(
  "ape", "Biostrings", "digest", "dplyr", "ggplot2",
  "phyloseq", "readr", "scales", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
stopifnot(length(missing_packages) == 0L)

pal_pub <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  sky = "#56B4E9",
  yellow = "#F0E442",
  grey = "#6B7280"
)

scale_color_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_color_manual(..., values = values)
}

scale_fill_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_fill_manual(..., values = values)
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
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", color = "#B3B3B3"
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(plot, stem, width = 125, height = 95, dpi = 350) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot,
    width = width, height = height, units = "mm",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot,
    width = width, height = height, units = "mm",
    dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), plot,
    width = width, height = height, units = "mm",
    dpi = dpi, compression = "lzw", bg = "white"
  )
  invisible(plot)
}

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/validate_phyloseq_import.R", "--project-root", ".", "--input-dir", "data/small/phyloseq-import", "--output-dir", "results/17-phyloseq-import", "--figure-dir", "figures"))

input_dir <- "data/small/phyloseq-import"

input_paths <- c(
  otutab = file.path(input_dir, "otutab.tsv"),
  taxonomy = file.path(input_dir, "taxonomy.tsv"),
  metadata = file.path(input_dir, "metadata.tsv"),
  representative_sequences = file.path(
    input_dir, "representative-sequences.fasta.gz"
  ),
  rooted_insertion_tree = file.path(
    input_dir, "rooted-insertion-tree.nwk.gz"
  )
)

expected_sha256 <- c(
  otutab = "6c43bb8f929f935b60e42ed32737632d3be18ddd491f5298ba3a6abd6e3f18a9",
  taxonomy = "8ef235ea19923307eef4d8e1905c2f34b78cf25cf1a83f08aa00353f2c5c0d38",
  metadata = "6628dd51051c53884e82c540ceaff8cc34f69900a78d372239a2ed41ff5b7263",
  representative_sequences = "a19e083cfb0ab757d557569bdf76aa3de78576ae35033a7e69a647976a370bc2",
  rooted_insertion_tree = "7e541e2ed9e3f81d9b6c6bf98999335e7f684a37ac5dd5f2b56e61b2f4c4c47c"
)

observed_sha256 <- vapply(
  input_paths,
  digest::digest,
  character(1),
  algo = "sha256",
  serialize = FALSE,
  file = TRUE
)

hash_audit <- tibble::tibble(
  asset = names(input_paths),
  observed_sha256 = unname(observed_sha256),
  expected_sha256 = unname(expected_sha256),
  status = ifelse(observed_sha256 == expected_sha256, "PASS", "FAIL")
)
stopifnot(all(hash_audit$status == "PASS"))
hash_audit

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(length(ids) > 0L, !anyDuplicated(ids), !anyNA(ids), nzchar(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

otutab <- read_keyed_tsv("data/small/phyloseq-import/otutab.tsv")
taxonomy <- read_keyed_tsv("data/small/phyloseq-import/taxonomy.tsv")
metadata <- read_keyed_tsv("data/small/phyloseq-import/metadata.tsv")

otu_matrix <- as.matrix(otutab)
storage.mode(otu_matrix) <- "numeric"
tax_matrix <- as.matrix(taxonomy)
storage.mode(tax_matrix) <- "character"

feature_ids <- rownames(otu_matrix)
sample_ids <- colnames(otu_matrix)

stopifnot(
  dim(otu_matrix)[[1L]] == 366L,
  dim(otu_matrix)[[2L]] == 4L,
  sum(otu_matrix) == 5213,
  all(is.finite(otu_matrix)),
  all(otu_matrix >= 0),
  all(otu_matrix == round(otu_matrix)),
  setequal(feature_ids, rownames(tax_matrix)),
  setequal(sample_ids, rownames(metadata))
)

list(
  otutab_dimension = dim(otu_matrix),
  taxonomy_dimension = dim(tax_matrix),
  metadata_dimension = dim(metadata),
  total_reads = sum(otu_matrix)
)

utils::head(
  data.frame(FeatureID = rownames(otutab), otutab, check.names = FALSE),
  5
)
utils::head(
  data.frame(FeatureID = rownames(taxonomy), taxonomy, check.names = FALSE),
  5
)
utils::head(
  data.frame(SampleID = rownames(metadata), metadata, check.names = FALSE),
  4
)

sample_audit <- tibble::tibble(
  SampleID = sample_ids,
  Library = metadata[sample_ids, "LibraryLabel"],
  Reads = as.numeric(colSums(otu_matrix)),
  ObservedASVs = as.numeric(colSums(otu_matrix > 0))
)

stopifnot(
  identical(sample_audit$Reads, c(1144, 1251, 1503, 1315)),
  identical(sample_audit$ObservedASVs, c(82, 133, 106, 138))
)
sample_audit

sequences <- Biostrings::readDNAStringSet(
  "data/small/phyloseq-import/representative-sequences.fasta.gz"
)
full_tree <- ape::read.tree(
  gzfile("data/small/phyloseq-import/rooted-insertion-tree.nwk.gz")
)

full_query_tips <- intersect(full_tree$tip.label, feature_ids)
reference_only_tips <- setdiff(full_tree$tip.label, feature_ids)

stopifnot(
  length(sequences) == 366L,
  setequal(names(sequences), feature_ids),
  ape::Ntip(full_tree) == 203818L,
  length(full_query_tips) == 366L,
  length(reference_only_tips) == 203452L,
  ape::is.rooted(full_tree),
  !anyNA(full_tree$edge.length),
  all(full_tree$edge.length >= 0)
)

query_tree <- ape::keep.tip(full_tree, feature_ids)

stopifnot(
  ape::Ntip(query_tree) == 366L,
  setequal(query_tree$tip.label, feature_ids),
  ape::is.rooted(query_tree),
  !anyNA(query_tree$edge.length),
  all(query_tree$edge.length >= 0)
)

tibble::tibble(
  Tree = c("Complete SEPP insertion tree", "ASV-pruned object tree"),
  Tips = c(ape::Ntip(full_tree), ape::Ntip(query_tree)),
  QueryTips = c(length(full_query_tips), ape::Ntip(query_tree)),
  ReferenceOnlyTips = c(length(reference_only_tips), 0L),
  Rooted = c(ape::is.rooted(full_tree), ape::is.rooted(query_tree))
)

id_gate <- tibble::tribble(
  ~Component, ~Domain, ~Expected, ~Matched, ~Missing, ~Unexpected,
  "Taxonomy", "Feature IDs", length(feature_ids),
    length(intersect(feature_ids, rownames(taxonomy))),
    length(setdiff(feature_ids, rownames(taxonomy))),
    length(setdiff(rownames(taxonomy), feature_ids)),
  "Representative sequences", "Feature IDs", length(feature_ids),
    length(intersect(feature_ids, names(sequences))),
    length(setdiff(feature_ids, names(sequences))),
    length(setdiff(names(sequences), feature_ids)),
  "Pruned tree", "Feature IDs", length(feature_ids),
    length(intersect(feature_ids, query_tree$tip.label)),
    length(setdiff(feature_ids, query_tree$tip.label)),
    length(setdiff(query_tree$tip.label, feature_ids)),
  "Metadata", "Sample IDs", length(sample_ids),
    length(intersect(sample_ids, rownames(metadata))),
    length(setdiff(sample_ids, rownames(metadata))),
    length(setdiff(rownames(metadata), sample_ids))
)

stopifnot(
  all(id_gate$Matched == id_gate$Expected),
  all(id_gate$Missing == 0L),
  all(id_gate$Unexpected == 0L)
)
id_gate

taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
sequences <- sequences[feature_ids]

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(taxonomy)),
  phyloseq::sample_data(metadata),
  phyloseq::phy_tree(query_tree),
  sequences
)

component_names <- methods::slotNames(ps)[
  !vapply(
    methods::slotNames(ps),
    function(slot_name) is.null(methods::slot(ps, slot_name)),
    logical(1)
  )
]

stopifnot(
  inherits(ps, "phyloseq"),
  setequal(
    component_names,
    c("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")
  ),
  phyloseq::taxa_are_rows(ps),
  phyloseq::ntaxa(ps) == 366L,
  phyloseq::nsamples(ps) == 4L,
  sum(phyloseq::sample_sums(ps)) == 5213,
  setequal(phyloseq::taxa_names(ps), feature_ids),
  setequal(phyloseq::sample_names(ps), sample_ids),
  setequal(phyloseq::phy_tree(ps)$tip.label, feature_ids),
  setequal(names(phyloseq::refseq(ps)), feature_ids)
)

object_summary <- tibble::tibble(
  Metric = c(
    "Class", "Components", "Features", "Samples", "Reads",
    "Tree tips", "Reference sequences", "Taxa are rows"
  ),
  Value = c(
    class(ps)[[1L]],
    length(component_names),
    phyloseq::ntaxa(ps),
    phyloseq::nsamples(ps),
    sum(phyloseq::sample_sums(ps)),
    ape::Ntip(phyloseq::phy_tree(ps)),
    length(phyloseq::refseq(ps)),
    phyloseq::taxa_are_rows(ps)
  )
)
object_summary

object_taxonomy <- as(
  phyloseq::tax_table(ps),
  "matrix"
)
taxa_reads <- phyloseq::taxa_sums(ps)

taxonomy_coverage <- do.call(
  rbind,
  lapply(colnames(object_taxonomy), function(rank_name) {
    values <- object_taxonomy[, rank_name]
    assigned <- !is.na(values) & nzchar(values) & values != "Unassigned"
    data.frame(
      Rank = rank_name,
      AssignedASVs = sum(assigned),
      TotalASVs = length(assigned),
      AssignedReads = sum(taxa_reads[assigned]),
      TotalReads = sum(taxa_reads),
      stringsAsFactors = FALSE
    )
  })
)

genus_row <- taxonomy_coverage[taxonomy_coverage$Rank == "Genus", ]
species_row <- taxonomy_coverage[taxonomy_coverage$Rank == "Species", ]
stopifnot(
  genus_row$AssignedASVs == 318L,
  genus_row$AssignedReads == 4217,
  species_row$AssignedASVs == 0L,
  species_row$AssignedReads == 0
)

taxonomy_coverage
phyloseq::sample_sums(ps)[sample_ids]
phyloseq::rank_names(ps)
phyloseq::sample_variables(ps)

rds_path <- "results/17-phyloseq-import/atacama-v4-phyloseq.rds"
dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(ps, rds_path, version = 3, compress = "xz")
ps_reloaded <- readRDS(rds_path)

stopifnot(
  inherits(ps_reloaded, "phyloseq"),
  phyloseq::ntaxa(ps_reloaded) == phyloseq::ntaxa(ps),
  phyloseq::nsamples(ps_reloaded) == phyloseq::nsamples(ps),
  identical(
    as.numeric(phyloseq::sample_sums(ps_reloaded)[sample_ids]),
    as.numeric(phyloseq::sample_sums(ps)[sample_ids])
  ),
  setequal(phyloseq::taxa_names(ps_reloaded), phyloseq::taxa_names(ps)),
  setequal(phyloseq::sample_names(ps_reloaded), phyloseq::sample_names(ps))
)

tibble::tibble(
  Metric = c("Class", "Features", "Samples", "Reads", "Taxa IDs", "Sample IDs"),
  Before = c(
    class(ps)[[1L]], phyloseq::ntaxa(ps), phyloseq::nsamples(ps),
    sum(phyloseq::sample_sums(ps)), length(phyloseq::taxa_names(ps)),
    length(phyloseq::sample_names(ps))
  ),
  After = c(
    class(ps_reloaded)[[1L]], phyloseq::ntaxa(ps_reloaded),
    phyloseq::nsamples(ps_reloaded), sum(phyloseq::sample_sums(ps_reloaded)),
    length(phyloseq::taxa_names(ps_reloaded)),
    length(phyloseq::sample_names(ps_reloaded))
  )
)
