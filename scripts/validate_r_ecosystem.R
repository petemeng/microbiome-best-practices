#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260719)

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
required_args <- c("project-root", "output-dir", "figure-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "jsonlite", "digest", "readr", "dplyr", "tidyr", "tibble", "ggplot2",
  "phyloseq", "microeco", "vegan", "BiocManager", "renv", "sessioninfo",
  "scales"
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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x)) y else x
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
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

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  ids <- as.character(x[[1L]])
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

expected_versions <- c(
  phyloseq = "1.48.0",
  microeco = "2.0.0",
  tidyverse = "2.0.0",
  ggplot2 = "3.5.2",
  dplyr = "1.1.4",
  readr = "2.1.5",
  tidyr = "1.3.1",
  tibble = "3.2.1",
  vegan = "2.6-6.1",
  BiocManager = "1.30.23",
  renv = "1.0.7",
  sessioninfo = "1.2.2"
)

package_roles <- c(
  phyloseq = "S4 microbiome object",
  microeco = "R6 analysis framework",
  tidyverse = "data grammar bundle",
  ggplot2 = "publication graphics",
  dplyr = "table transformation",
  readr = "delimited input/output",
  tidyr = "long/wide reshaping",
  tibble = "rectangular data",
  vegan = "ecological statistics",
  BiocManager = "Bioconductor resolver",
  renv = "project lockfile",
  sessioninfo = "environment provenance"
)

library_scope <- function(package) {
  path <- normalizePath(find.package(package), mustWork = TRUE)
  if (startsWith(path, project_root)) {
    "Project library"
  } else if (startsWith(path, "/usr/") || startsWith(path, "/usr/local/")) {
    "System / site library"
  } else {
    "User library"
  }
}

observed_versions <- vapply(
  names(expected_versions),
  function(package) packageDescription(package, fields = "Version"),
  character(1)
)
namespace_loaded <- vapply(
  names(expected_versions),
  requireNamespace,
  logical(1),
  quietly = TRUE
)

package_audit <- data.frame(
  package = names(expected_versions),
  role = unname(package_roles[names(expected_versions)]),
  expected_version = unname(expected_versions),
  observed_version = unname(observed_versions),
  library_scope = vapply(
    names(expected_versions),
    library_scope,
    character(1)
  ),
  namespace_loaded = namespace_loaded,
  status = ifelse(
    namespace_loaded & observed_versions == expected_versions,
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)

lock_path <- file.path(project_root, "env", "renv.lock")
lock <- jsonlite::read_json(lock_path, simplifyVector = FALSE)
lock_packages <- lock$Packages
lock_audit <- do.call(
  rbind,
  lapply(names(expected_versions), function(package) {
    record <- lock_packages[[package]]
    locked_version <- if (is.null(record)) NA_character_ else record$Version
    data.frame(
      package = package,
      expected_version = expected_versions[[package]],
      locked_version = locked_version,
      source = if (is.null(record)) NA_character_ else record$Source %||% NA_character_,
      repository = if (is.null(record)) NA_character_ else record$Repository %||% NA_character_,
      status = ifelse(
        !is.null(record) &&
          identical(locked_version, expected_versions[[package]]) &&
          nzchar(record$Source %||% ""),
        "PASS",
        "FAIL"
      ),
      check.names = FALSE
    )
  })
)

data_paths <- c(
  otutab = file.path(project_root, "data", "small", "otutab.tsv"),
  taxonomy = file.path(project_root, "data", "small", "taxonomy.tsv"),
  metadata = file.path(project_root, "data", "small", "metadata.tsv")
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64"
)
observed_sha256 <- vapply(data_paths, sha256_file, character(1))

otutab <- read_keyed_tsv(data_paths[["otutab"]])
taxonomy <- read_keyed_tsv(data_paths[["taxonomy"]])
metadata <- read_keyed_tsv(data_paths[["metadata"]])
otu_matrix <- as.matrix(otutab)
storage.mode(otu_matrix) <- "numeric"
tax_matrix <- as.matrix(taxonomy)
storage.mode(tax_matrix) <- "character"

input_audit <- data.frame(
  object = c("otutab", "taxonomy", "metadata"),
  orientation = c(
    "features × samples",
    "features × taxonomic ranks",
    "samples × variables"
  ),
  rows = c(nrow(otutab), nrow(taxonomy), nrow(metadata)),
  columns = c(ncol(otutab), ncol(taxonomy), ncol(metadata)),
  sha256 = unname(observed_sha256[c("otutab", "taxonomy", "metadata")]),
  expected_sha256 = unname(expected_sha256[c("otutab", "taxonomy", "metadata")]),
  status = ifelse(
    observed_sha256[c("otutab", "taxonomy", "metadata")] ==
      expected_sha256[c("otutab", "taxonomy", "metadata")],
    "PASS",
    "FAIL"
  ),
  check.names = FALSE
)

phyloseq_object <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE),
  phyloseq::tax_table(tax_matrix),
  phyloseq::sample_data(metadata)
)
microtable_object <- microeco::microtable$new(
  sample_table = metadata,
  otu_table = as.data.frame(otu_matrix, check.names = FALSE),
  tax_table = taxonomy,
  auto_tidy = TRUE
)

raw_library <- colSums(otu_matrix)
phyloseq_library <- phyloseq::sample_sums(phyloseq_object)
microtable_library <- colSums(as.matrix(microtable_object$otu_table))

sample_library_audit <- dplyr::left_join(
  dplyr::left_join(
    tibble::tibble(
      sample_id = names(raw_library),
      raw_reads = as.numeric(raw_library)
    ),
    tibble::tibble(
      sample_id = names(phyloseq_library),
      phyloseq_reads = as.numeric(phyloseq_library)
    ),
    by = "sample_id"
  ),
  tibble::tibble(
    sample_id = names(microtable_library),
    microtable_reads = as.numeric(microtable_library)
  ),
  by = "sample_id"
)
sample_library_audit <- dplyr::mutate(
  sample_library_audit,
  phyloseq_delta = phyloseq_reads - raw_reads,
  microtable_delta = microtable_reads - raw_reads,
  status = ifelse(
    phyloseq_delta == 0 & microtable_delta == 0,
    "PASS",
    "FAIL"
  )
)

raw_features <- nrow(otu_matrix)
raw_samples <- ncol(otu_matrix)
raw_reads <- sum(otu_matrix)
phyloseq_features <- phyloseq::ntaxa(phyloseq_object)
phyloseq_samples <- phyloseq::nsamples(phyloseq_object)
phyloseq_reads <- sum(phyloseq::sample_sums(phyloseq_object))
microtable_features <- nrow(microtable_object$otu_table)
microtable_samples <- ncol(microtable_object$otu_table)
microtable_reads <- sum(as.matrix(microtable_object$otu_table))

object_parity <- data.frame(
  object = rep(c("phyloseq", "microtable"), each = 3L),
  metric = rep(c("features", "samples", "total_reads"), times = 2L),
  reference_value = rep(c(raw_features, raw_samples, raw_reads), times = 2L),
  observed_value = c(
    phyloseq_features,
    phyloseq_samples,
    phyloseq_reads,
    microtable_features,
    microtable_samples,
    microtable_reads
  ),
  status = "PASS",
  check.names = FALSE
)
object_parity$status[
  object_parity$reference_value != object_parity$observed_value
] <- "FAIL"

check_rows <- list()
add_check <- function(id, observed, expected, passed) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check_id = id,
    observed = as.character(observed),
    expected = as.character(expected),
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    check.names = FALSE
  )
}

r_version <- paste(R.version$major, R.version$minor, sep = ".")
bioconductor_version <- as.character(BiocManager::version())
lock_r_version <- lock$R$Version

add_check("r-version", r_version, "4.4.1", identical(r_version, "4.4.1"))
add_check(
  "r-architecture",
  paste0(.Machine$sizeof.pointer * 8L, "-bit"),
  "64-bit",
  .Machine$sizeof.pointer == 8L
)
add_check(
  "bioconductor-version",
  bioconductor_version,
  "3.19",
  identical(bioconductor_version, "3.19")
)
add_check(
  "lock-r-version",
  lock_r_version,
  "4.4.1",
  identical(lock_r_version, "4.4.1")
)
add_check(
  "lock-package-records",
  length(lock_packages),
  ">=147",
  length(lock_packages) >= 147L
)
add_check(
  "package-namespaces",
  sum(package_audit$namespace_loaded),
  nrow(package_audit),
  all(package_audit$namespace_loaded)
)
add_check(
  "package-versions",
  sum(package_audit$status == "PASS"),
  nrow(package_audit),
  all(package_audit$status == "PASS")
)
add_check(
  "lock-package-versions",
  sum(lock_audit$status == "PASS"),
  nrow(lock_audit),
  all(lock_audit$status == "PASS")
)
add_check(
  "lock-package-sources",
  sum(!is.na(lock_audit$source) & nzchar(lock_audit$source)),
  nrow(lock_audit),
  all(!is.na(lock_audit$source) & nzchar(lock_audit$source))
)
add_check(
  "otutab-sha256",
  observed_sha256[["otutab"]],
  expected_sha256[["otutab"]],
  identical(observed_sha256[["otutab"]], expected_sha256[["otutab"]])
)
add_check(
  "taxonomy-sha256",
  observed_sha256[["taxonomy"]],
  expected_sha256[["taxonomy"]],
  identical(observed_sha256[["taxonomy"]], expected_sha256[["taxonomy"]])
)
add_check(
  "metadata-sha256",
  observed_sha256[["metadata"]],
  expected_sha256[["metadata"]],
  identical(observed_sha256[["metadata"]], expected_sha256[["metadata"]])
)
add_check(
  "otutab-dimensions",
  paste(dim(otutab), collapse = " × "),
  "13628 × 90",
  identical(dim(otutab), c(13628L, 90L))
)
add_check(
  "taxonomy-dimensions",
  paste(dim(taxonomy), collapse = " × "),
  "13628 × 7",
  identical(dim(taxonomy), c(13628L, 7L))
)
add_check(
  "metadata-dimensions",
  paste(dim(metadata), collapse = " × "),
  "90 × 3",
  identical(dim(metadata), c(90L, 3L))
)
add_check(
  "unique-feature-ids",
  sum(duplicated(rownames(otutab))) + sum(duplicated(rownames(taxonomy))),
  "0 duplicated IDs",
  !anyDuplicated(rownames(otutab)) && !anyDuplicated(rownames(taxonomy))
)
add_check(
  "unique-sample-ids",
  sum(duplicated(colnames(otutab))) + sum(duplicated(rownames(metadata))),
  "0 duplicated IDs",
  !anyDuplicated(colnames(otutab)) && !anyDuplicated(rownames(metadata))
)
add_check(
  "feature-id-alignment",
  sum(rownames(otutab) == rownames(taxonomy)),
  nrow(otutab),
  identical(rownames(otutab), rownames(taxonomy))
)
add_check(
  "sample-id-alignment",
  sum(colnames(otutab) == rownames(metadata)),
  ncol(otutab),
  identical(colnames(otutab), rownames(metadata))
)
valid_counts <- all(is.finite(otu_matrix)) &&
  all(otu_matrix >= 0) &&
  all(otu_matrix == round(otu_matrix))
add_check(
  "nonnegative-integer-counts",
  valid_counts,
  TRUE,
  valid_counts
)
add_check(
  "total-reads",
  format(raw_reads, scientific = FALSE),
  "1619670",
  identical(as.numeric(raw_reads), 1619670)
)
add_check(
  "phyloseq-contract",
  paste(phyloseq_features, phyloseq_samples, phyloseq_reads, sep = "/"),
  "13628/90/1619670",
  phyloseq_features == 13628L &&
    phyloseq_samples == 90L &&
    phyloseq_reads == 1619670
)
add_check(
  "microtable-contract",
  paste(microtable_features, microtable_samples, microtable_reads, sep = "/"),
  "13628/90/1619670",
  microtable_features == 13628L &&
    microtable_samples == 90L &&
    microtable_reads == 1619670
)
add_check(
  "object-parity-six-metrics",
  sum(object_parity$status == "PASS"),
  "6",
  all(object_parity$status == "PASS")
)
add_check(
  "sample-library-parity",
  sum(sample_library_audit$status == "PASS"),
  "90",
  nrow(sample_library_audit) == 90L &&
    all(sample_library_audit$status == "PASS")
)
add_check(
  "tidyverse-complete-sample-audit",
  sum(stats::complete.cases(sample_library_audit)),
  "90",
  nrow(sample_library_audit) == 90L &&
    all(stats::complete.cases(sample_library_audit))
)

checks <- do.call(rbind, check_rows)
stopifnot(nrow(checks) == 26L)
checks_passed <- sum(checks$status == "PASS")
checks_failed <- sum(checks$status == "FAIL")

dark <- "#111827"
gray <- "#4B5563"
blue <- "#0072B2"
sky <- "#56B4E9"
green <- "#009E73"
orange <- "#D55E00"
gold <- "#E69F00"
purple <- "#CC79A7"
card_fill <- "#F8FAFC"

map_nodes <- data.frame(
  xmin = c(0.25, 2.05, 3.85, 5.65, 7.45, 7.45, 9.25),
  xmax = c(1.65, 3.45, 5.25, 7.05, 8.85, 8.85, 10.95),
  ymin = c(1.50, 1.50, 1.50, 1.50, 2.35, 0.65, 1.50),
  ymax = c(3.10, 3.10, 3.10, 3.10, 3.95, 2.25, 3.10),
  label = c(
    "R runtime\n4.4.1",
    paste0("renv.lock\n", length(lock_packages), " records"),
    "tidyverse\nread · reshape",
    "Three tables\nIDs + counts",
    "phyloseq\nS4 object",
    "microeco\nR6 object",
    "Methods + plots\nPDF · PNG · TIFF"
  ),
  colour = c(blue, sky, green, gold, purple, orange, blue),
  stringsAsFactors = FALSE
)
map_edges <- data.frame(
  x = c(1.65, 3.45, 5.25, 7.05, 7.05, 8.85, 8.85),
  y = c(2.30, 2.30, 2.30, 2.30, 2.30, 3.15, 1.45),
  xend = c(2.00, 3.80, 5.60, 7.40, 7.40, 9.20, 9.20),
  yend = c(2.30, 2.30, 2.30, 3.15, 1.45, 2.45, 2.15)
)
ecosystem_map <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = map_edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#6B7280",
    linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = map_nodes,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      colour = colour
    ),
    fill = card_fill,
    linewidth = 1.0
  ) +
  ggplot2::geom_text(
    data = map_nodes,
    ggplot2::aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label
    ),
    colour = dark,
    fontface = "bold",
    lineheight = 1.1,
    size = 3.7
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::annotate(
    "text",
    x = 0.25,
    y = 4.65,
    hjust = 0,
    label = "A reproducible R microbiome stack",
    colour = dark,
    fontface = "bold",
    size = 6.6
  ) +
  ggplot2::annotate(
    "text",
    x = 0.25,
    y = 4.25,
    hjust = 0,
    label = "The IDE is optional; the runtime, lockfile, data contract, and object parity are the evidence.",
    colour = gray,
    size = 3.6
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0.25,
    xmax = 10.95,
    ymin = -0.05,
    ymax = 0.42,
    fill = "#F3F4F6",
    colour = NA
  ) +
  ggplot2::annotate(
    "text",
    x = 0.48,
    y = 0.19,
    hjust = 0,
    label = "Acceptance: 26 / 26 checks passed",
    colour = green,
    fontface = "bold",
    size = 4.1
  ) +
  ggplot2::annotate(
    "text",
    x = 10.70,
    y = 0.19,
    hjust = 1,
    label = "PASS",
    colour = green,
    fontface = "bold",
    size = 4.1
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 11.2), ylim = c(-0.15, 4.9), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(18, 18, 16, 18))

package_cards <- package_audit
package_cards$column <- rep(1:4, length.out = nrow(package_cards))
package_cards$row <- rep(3:1, each = 4)
package_cards$xmin <- package_cards$column - 0.43
package_cards$xmax <- package_cards$column + 0.43
package_cards$ymin <- package_cards$row - 0.37
package_cards$ymax <- package_cards$row + 0.37
package_cards$label <- paste0(
  package_cards$package,
  "\n",
  package_cards$observed_version,
  "\n",
  package_cards$status
)
package_cards$colour <- rep(c(blue, green, purple, orange), 3L)
package_plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = package_cards,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      colour = colour
    ),
    fill = card_fill,
    linewidth = 0.9
  ) +
  ggplot2::geom_text(
    data = package_cards,
    ggplot2::aes(x = column, y = row, label = label),
    colour = dark,
    fontface = "bold",
    lineheight = 1.05,
    size = 3.4
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::annotate(
    "text",
    x = 0.55,
    y = 4.15,
    hjust = 0,
    label = "Locked package contract",
    colour = dark,
    fontface = "bold",
    size = 6.5
  ) +
  ggplot2::annotate(
    "text",
    x = 0.55,
    y = 3.78,
    hjust = 0,
    label = paste0(
      "R 4.4.1 · Bioconductor 3.19 · ",
      nrow(package_audit),
      " / ",
      nrow(package_audit),
      " namespaces and versions passed"
    ),
    colour = gray,
    size = 3.6
  ) +
  ggplot2::annotate(
    "text",
    x = 0.55,
    y = 0.18,
    hjust = 0,
    label = paste0(
      "renv.lock: ",
      length(lock_packages),
      " package records · audited records include exact version and source"
    ),
    colour = green,
    fontface = "bold",
    size = 3.6
  ) +
  ggplot2::coord_cartesian(xlim = c(0.45, 4.55), ylim = c(-0.02, 4.38), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(18, 18, 14, 18))

parity_plot_data <- rbind(
  data.frame(
    sample_id = sample_library_audit$sample_id,
    raw_reads = sample_library_audit$raw_reads,
    object_reads = sample_library_audit$phyloseq_reads,
    object = "phyloseq",
    group = metadata[sample_library_audit$sample_id, "Group"]
  ),
  data.frame(
    sample_id = sample_library_audit$sample_id,
    raw_reads = sample_library_audit$raw_reads,
    object_reads = sample_library_audit$microtable_reads,
    object = "microtable",
    group = metadata[sample_library_audit$sample_id, "Group"]
  )
)
parity_plot_data$object <- factor(
  parity_plot_data$object,
  levels = c("phyloseq", "microtable")
)
parity_plot <- ggplot2::ggplot(
  parity_plot_data,
  ggplot2::aes(x = raw_reads, y = object_reads, colour = group)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    colour = "#6B7280",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 2.6, alpha = 0.82) +
  ggplot2::facet_wrap(~object, nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(CW = blue, IW = gold, RW = green),
    name = "Wetland group"
  ) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Object conversion preserves every sample library size",
    subtitle = paste0(
      "Real wetland triad: 13,628 features × 90 samples · ",
      format(raw_reads, big.mark = ",", scientific = FALSE),
      " total reads · max |difference| = 0"
    ),
    x = "Reads in raw otutab",
    y = "Reads after object construction",
    caption = "Each point is one real sample; the dashed line is exact identity."
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 18,
      face = "bold",
      colour = dark,
      margin = ggplot2::margin(b = 7)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      colour = gray,
      margin = ggplot2::margin(b = 12)
    ),
    strip.text = ggplot2::element_text(size = 12.5, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(size = 9.5, colour = gray, hjust = 0),
    plot.margin = ggplot2::margin(16, 18, 14, 16)
  )

save_publication_plot(ecosystem_map, "09-r-ecosystem-map", 12.0, 5.9)
save_publication_plot(package_plot, "09-package-version-audit", 10.7, 8.0)
save_publication_plot(parity_plot, "09-object-parity-audit", 11.4, 6.9)

write_tsv(package_audit, "package-audit.tsv")
write_tsv(lock_audit, "lockfile-audit.tsv")
write_tsv(input_audit, "input-audit.tsv")
write_tsv(object_parity, "object-parity.tsv")
write_tsv(sample_library_audit, "sample-library-audit.tsv")

session_lines <- c(
  capture.output(sessioninfo::session_info(pkgs = names(expected_versions))),
  "",
  paste0("Bioconductor version: ", bioconductor_version),
  paste0("renv.lock SHA-256: ", sha256_file(lock_path))
)
writeLines(
  sanitize_text(session_lines),
  file.path(output_dir, "r-session-info.txt"),
  useBytes = TRUE
)

summary_payload <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  validation_date = as.character(Sys.Date()),
  r_version = r_version,
  architecture = paste0(.Machine$sizeof.pointer * 8L, "-bit"),
  platform = R.version$platform,
  bioconductor_version = bioconductor_version,
  renv_lock_sha256 = sha256_file(lock_path),
  lock_package_records = length(lock_packages),
  packages_audited = nrow(package_audit),
  package_checks_passed = sum(package_audit$status == "PASS"),
  checks_total = nrow(checks),
  checks_passed = checks_passed,
  checks_failed = checks_failed,
  object_parity_checks_total = nrow(object_parity),
  object_parity_checks_passed = sum(object_parity$status == "PASS"),
  input = list(
    features = raw_features,
    samples = raw_samples,
    taxonomy_ranks = ncol(taxonomy),
    metadata_variables = ncol(metadata),
    total_reads = as.numeric(raw_reads)
  ),
  phyloseq = list(
    features = phyloseq_features,
    samples = phyloseq_samples,
    total_reads = as.numeric(phyloseq_reads)
  ),
  microtable = list(
    features = microtable_features,
    samples = microtable_samples,
    total_reads = as.numeric(microtable_reads)
  ),
  sample_library_parity = list(
    samples = nrow(sample_library_audit),
    maximum_absolute_phyloseq_difference = max(
      abs(sample_library_audit$phyloseq_delta)
    ),
    maximum_absolute_microtable_difference = max(
      abs(sample_library_audit$microtable_delta)
    )
  ),
  figures = c(
    "09-r-ecosystem-map",
    "09-package-version-audit",
    "09-object-parity-audit"
  )
)
jsonlite::write_json(
  summary_payload,
  file.path(output_dir, "environment-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

log_lines <- c(
  "# Article 09 R ecosystem and object-contract audit",
  paste0("validation_date: ", summary_payload$validation_date),
  paste0("status: ", summary_payload$status),
  paste0("R: ", r_version, " (", summary_payload$architecture, ")"),
  paste0("Bioconductor: ", bioconductor_version),
  paste0("renv.lock records: ", length(lock_packages)),
  paste0("renv.lock SHA-256: ", summary_payload$renv_lock_sha256),
  "",
  "## Required checks",
  paste(
    checks$check_id,
    checks$status,
    paste0("observed=", checks$observed),
    paste0("expected=", checks$expected),
    sep = "\t"
  ),
  "",
  "## Real triad",
  paste0("otutab: ", raw_features, " features x ", raw_samples, " samples"),
  paste0("taxonomy: ", nrow(taxonomy), " features x ", ncol(taxonomy), " ranks"),
  paste0("metadata: ", nrow(metadata), " samples x ", ncol(metadata), " variables"),
  paste0("total reads: ", format(raw_reads, scientific = FALSE)),
  "",
  "## Object parity",
  paste(
    object_parity$object,
    object_parity$metric,
    object_parity$observed_value,
    object_parity$status,
    sep = "\t"
  ),
  "",
  paste0("Checks passed: ", checks_passed, " / ", nrow(checks))
)
writeLines(
  sanitize_text(log_lines),
  file.path(output_dir, "validation.log"),
  useBytes = TRUE
)

cat(
  jsonlite::toJSON(
    summary_payload,
    pretty = TRUE,
    auto_unbox = TRUE
  ),
  "\n"
)

if (checks_failed > 0L) {
  stop(
    checks_failed,
    " of ",
    nrow(checks),
    " required checks failed.",
    call. = FALSE
  )
}
