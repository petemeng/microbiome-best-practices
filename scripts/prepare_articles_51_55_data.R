#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(phyloseq)
})

root <- normalizePath(".", mustWork = TRUE)
raw_archive <- file.path(root, "data", "raw", "microbiome_1.26.0.tar.gz")
out_dir <- file.path(root, "data", "small", "longitudinal-dietswap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_sha256 <- "67cc6d117ce9dfbbc472fba828e1bb0d48c621fb574f4f8d3424f6cc9a77fdab"
observed_sha256 <- digest::digest(raw_archive, algo = "sha256", file = TRUE)
if (!identical(observed_sha256, expected_sha256)) {
  stop("Unexpected microbiome 1.26.0 source archive checksum: ", observed_sha256)
}

temporary_directory <- tempfile("dietswap-")
dir.create(temporary_directory)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)
utils::untar(
  raw_archive,
  files = "microbiome/data/dietswap.rda",
  exdir = temporary_directory
)
data_environment <- new.env(parent = emptyenv())
load(
  file.path(temporary_directory, "microbiome", "data", "dietswap.rda"),
  envir = data_environment
)
dietswap <- data_environment$dietswap
stopifnot(
  inherits(dietswap, "phyloseq"),
  phyloseq::taxa_are_rows(dietswap),
  phyloseq::ntaxa(dietswap) == 130L,
  phyloseq::nsamples(dietswap) == 222L
)

original_names <- phyloseq::taxa_names(dietswap)
feature_ids <- sprintf("HT_%03d", seq_along(original_names))

otutab <- as(phyloseq::otu_table(dietswap), "matrix")
storage.mode(otutab) <- "numeric"
rownames(otutab) <- feature_ids
stopifnot(all(is.finite(otutab)), all(otutab >= 0))

original_taxonomy <- as.data.frame(
  as(phyloseq::tax_table(dietswap), "matrix"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
taxonomy <- data.frame(
  Kingdom = "Bacteria",
  Phylum = original_taxonomy$Phylum,
  Class = "",
  Order = "",
  Family = original_taxonomy$Family,
  Genus = original_taxonomy$Genus,
  Species = "",
  DisplayName = original_names,
  stringsAsFactors = FALSE,
  row.names = feature_ids
)

original_metadata <- as.data.frame(
  phyloseq::sample_data(dietswap),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
metadata <- data.frame(
  SubjectID = original_metadata$subject,
  Sex = tools::toTitleCase(as.character(original_metadata$sex)),
  Cohort = ifelse(
    as.character(original_metadata$nationality) == "AAM",
    "African American",
    "Rural African"
  ),
  OriginalCohortCode = as.character(original_metadata$nationality),
  PhaseCode = as.character(original_metadata$group),
  Timepoint = as.integer(original_metadata$timepoint),
  TimeWithinPhase = as.integer(original_metadata$timepoint.within.group),
  BMIGroup = tools::toTitleCase(as.character(original_metadata$bmi_group)),
  stringsAsFactors = FALSE,
  row.names = rownames(original_metadata)
)
rownames(metadata) <- original_metadata$sample
stopifnot(
  identical(colnames(otutab), rownames(metadata)),
  length(unique(metadata$SubjectID)) == 38L,
  identical(sort(unique(metadata$Timepoint)), 1:6)
)

write.table(
  cbind(FeatureID = rownames(otutab), as.data.frame(otutab, check.names = FALSE)),
  file.path(out_dir, "otutab.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cbind(FeatureID = rownames(taxonomy), taxonomy),
  file.path(out_dir, "taxonomy.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  cbind(SampleID = rownames(metadata), metadata),
  file.path(out_dir, "metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

summary <- list(
  dataset = "O'Keefe_diet_swap_2015",
  primary_citation = "O'Keefe et al. Nature Communications 2015; doi:10.1038/ncomms7342",
  primary_data_repository = "https://doi.org/10.5061/dryad.1mn1n",
  packaged_source = "Bioconductor microbiome 1.26.0 dietswap data object",
  packaged_source_url = "https://bioconductor.org/packages/3.19/bioc/src/contrib/microbiome_1.26.0.tar.gz",
  packaged_source_license = "BSD-2-Clause",
  packaged_source_sha256 = observed_sha256,
  assay = "HITChip phylogenetic 16S rRNA microarray abundance profiles",
  features = nrow(otutab),
  samples = ncol(otutab),
  subjects = length(unique(metadata$SubjectID)),
  cohort_counts = as.list(table(metadata$Cohort)),
  timepoint_counts = as.list(table(metadata$Timepoint)),
  abundance_total = unname(sum(otutab)),
  contract = "otutab/taxonomy/metadata; features x samples; exact sample order shared across files",
  note = "The assay is a phylogenetic 16S microarray rather than amplicon sequencing; the longitudinal model and plotting workflow transfer, but assay-specific preprocessing does not."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "source-summary.json"),
  useBytes = TRUE
)

message(
  "Prepared DietSwap three-piece bundle: ",
  nrow(otutab), " features x ", ncol(otutab), " samples; ",
  length(unique(metadata$SubjectID)), " subjects."
)
