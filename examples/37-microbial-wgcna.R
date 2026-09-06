# WGCNA 用于微生物组：模块、性状关联与候选 hub
# Run sequentially in a new working directory.
# Required packages: WGCNA, ggdendro, ggplot2, ggrepel, knitr, patchwork, ragg, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/environment.tsv",
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(WGCNA)

set.seed(20260737)
font_pub <- "sans"
pal_pub <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
  yellow = "#F0E442", grey = "#6B7280"
)
scale_color_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_color_manual(..., values = values)
}
scale_fill_pub <- function(..., values = unname(pal_pub)) {
  ggplot2::scale_fill_manual(..., values = values)
}
theme_pub <- function(base_size = 10, base_family = font_pub) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "#E6E6E6", linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      axis.ticks = ggplot2::element_line(colour = "#1A1A1A", linewidth = 0.3),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(), legend.position = "top"
    )
}
save_pub <- function(
  plot, file_base, width = 89, height = 70, units = "mm",
  dpi = 600, write_svg = TRUE, write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot, width = width, height = height,
    units = units, device = grDevices::cairo_pdf,
    family = base_family, bg = "white"
  )
  if (write_svg) ggplot2::ggsave(
    paste0(file_base, ".svg"), plot, width = width, height = height,
    units = units, device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot, width = width, height = height,
    units = units, dpi = dpi, device = ragg::agg_png, bg = "white"
  )
  if (write_tiff) ggplot2::ggsave(
    paste0(file_base, ".tiff"), plot, width = width, height = height,
    units = units, dpi = dpi, device = ragg::agg_tiff,
    compression = "lzw", bg = "white"
  )
  invisible(plot)
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path, show_col_types = FALSE, progress = FALSE,
    name_repair = "minimal", na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}
otutab_raw <- read_keyed_tsv("data/small/otutab.tsv")
taxonomy <- read_keyed_tsv("data/small/taxonomy.tsv")
metadata <- read_keyed_tsv("data/small/metadata.tsv")
environment <- read_keyed_tsv("data/small/environment.tsv")
otutab <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
metadata$Group <- factor(metadata$Group, levels = c("CW", "IW", "TW"))
environment <- environment[rownames(metadata), , drop = FALSE]

rank_names <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
clean_taxon <- function(x) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  x[
    is.na(x) | x == "" |
      grepl(
        "unclassified|unknown|unassigned|uncultured|metagenome",
        x, ignore.case = TRUE
      )
  ] <- NA_character_
  x
}
taxonomy_clean <- as.data.frame(
  lapply(taxonomy[rownames(otutab), rank_names], clean_taxon),
  check.names = FALSE
)
known_genus <- !is.na(taxonomy_clean$Genus)
lineage_key <- apply(
  taxonomy_clean[known_genus, rank_names[1:6], drop = FALSE], 1L,
  function(x) paste(ifelse(is.na(x), "?", x), collapse = ";")
)
genus_counts <- rowsum(
  otutab[known_genus, , drop = FALSE], lineage_key, reorder = FALSE
)
lineage_matrix <- do.call(
  rbind, strsplit(rownames(genus_counts), ";", fixed = TRUE)
)
colnames(lineage_matrix) <- rank_names[1:6]
feature_ids <- make.unique(paste0("Genus:", lineage_matrix[, "Genus"]))
rownames(genus_counts) <- feature_ids

eligible <- rowMeans(genus_counts > 0) >= 0.20 & rowSums(genus_counts) >= 100
counts <- genus_counts[eligible, rownames(metadata), drop = FALSE]
make_clr <- function(pseudocount) {
  log_values <- log2(t(counts) + pseudocount)
  log_values - rowMeans(log_values)
}
dat_expr <- make_clr(0.5)
quality <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)

stopifnot(
  nrow(otutab) == 13628L,
  ncol(otutab) == 90L,
  sum(otutab) == 1619670L,
  nrow(dat_expr) == 90L,
  ncol(dat_expr) == 422L,
  quality$allOK,
  identical(rownames(dat_expr), rownames(metadata)),
  all(is.finite(dat_expr))
)

data.frame(
  Samples = nrow(dat_expr),
  NamedGenera = ncol(dat_expr),
  RetainedAllReadFraction = sum(counts) / sum(otutab),
  FailedSamples = sum(!quality$goodSamples),
  FailedGenera = sum(!quality$goodGenes)
)

powers <- c(1:12, 14, 16, 18, 20)
soft_threshold <- WGCNA::pickSoftThreshold(
  dat_expr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.10),
  verbose = 0
)
power_table <- soft_threshold$fitIndices
eligible_power <- with(
  power_table,
  Power[SFT.R.sq >= 0.80 & mean.k. >= 5]
)
if (length(eligible_power)) {
  selected_power <- min(eligible_power)
  power_gate <- "PASS"
} else {
  connectivity_ok <- power_table$mean.k. >= 5
  stopifnot(any(connectivity_ok))
  selected_power <- power_table$Power[
    which.max(ifelse(connectivity_ok, power_table$SFT.R.sq, -Inf))
  ]
  power_gate <- "FALLBACK"
}
selected_power_row <- power_table[power_table$Power == selected_power, ]

data.frame(
  SelectedPower = selected_power,
  SignedR2 = selected_power_row$SFT.R.sq,
  MeanConnectivity = selected_power_row$mean.k.,
  Gate = power_gate
)

set.seed(20260737)
wgcna_fit <- WGCNA::blockwiseModules(
  dat_expr,
  power = selected_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.10,
  deepSplit = 2,
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  randomSeed = 20260737L,
  nThreads = 4,
  verbose = 0
)
module_colours <- WGCNA::labels2colors(wgcna_fit$colors)
module_eigengenes <- WGCNA::orderMEs(WGCNA::moduleEigengenes(
  dat_expr,
  colors = module_colours,
  excludeGrey = TRUE
)$eigengenes)
module_sizes <- as.data.frame(table(module_colours), stringsAsFactors = FALSE)
names(module_sizes) <- c("Module", "Genera")

stopifnot(
  length(module_colours) == ncol(dat_expr),
  ncol(module_eigengenes) >= 2L,
  all(rownames(module_eigengenes) == rownames(dat_expr)),
  sum(module_sizes$Genera) == ncol(dat_expr)
)
module_sizes

trait_data <- data.frame(
  CW = as.integer(metadata$Group == "CW"),
  IW = as.integer(metadata$Group == "IW"),
  TW = as.integer(metadata$Group == "TW"),
  pH = as.numeric(environment$pH),
  TOC = as.numeric(environment$TOC),
  Conductivity = as.numeric(environment$Conductivity),
  Temperature = as.numeric(environment$Temperature),
  Precipitation = as.numeric(environment$Precipitation),
  row.names = rownames(metadata),
  check.names = FALSE
)
module_trait_r <- WGCNA::bicor(
  module_eigengenes,
  trait_data,
  use = "pairwise.complete.obs",
  maxPOutliers = 0.10,
  robustY = FALSE
)
module_trait_p <- WGCNA::corPvalueStudent(
  module_trait_r, nSamples = nrow(dat_expr)
)
module_trait_q <- matrix(
  stats::p.adjust(module_trait_p, method = "BH"),
  nrow = nrow(module_trait_p),
  dimnames = dimnames(module_trait_p)
)
module_trait_ledger <- expand.grid(
  Module = rownames(module_trait_r),
  Trait = colnames(module_trait_r),
  stringsAsFactors = FALSE
)
module_trait_ledger$Correlation <- as.vector(module_trait_r)
module_trait_ledger$PValue <- as.vector(module_trait_p)
module_trait_ledger$QValue <- as.vector(module_trait_q)
module_trait_ledger$Significant <- module_trait_ledger$QValue < 0.05

data.frame(
  Modules = nrow(module_trait_r),
  Traits = ncol(module_trait_r),
  Tests = length(module_trait_p),
  BHSignificant = sum(module_trait_q < 0.05)
)

adjacency_matrix <- WGCNA::adjacency(
  dat_expr,
  power = selected_power,
  type = "signed",
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.10)
)
module_membership <- WGCNA::signedKME(
  dat_expr,
  module_eigengenes,
  corFnc = "bicor",
  corOptions = 'use="p", maxPOutliers=0.1'
)
hub_ledger <- data.frame(
  FeatureID = colnames(dat_expr),
  DisplayTaxon = sub("Genus:", "", colnames(dat_expr), fixed = TRUE),
  Module = module_colours,
  stringsAsFactors = FALSE
)
hub_ledger$WithinConnectivity <- vapply(
  seq_len(ncol(dat_expr)),
  function(index) {
    sum(adjacency_matrix[index, module_colours == module_colours[[index]]]) - 1
  },
  numeric(1)
)
membership_columns <- match(
  paste0("kME", hub_ledger$Module), colnames(module_membership)
)
hub_ledger$KME <- as.numeric(module_membership[cbind(
  seq_len(nrow(hub_ledger)), membership_columns
)])
hub_ledger$CandidateHub <- with(
  hub_ledger,
  Module != "grey" & abs(KME) >= 0.80 &
    ave(
      WithinConnectivity, Module,
      FUN = function(x) x >= stats::quantile(x, 0.90)
    )
)

set.seed(20260737)
wgcna_pseudocount_one <- WGCNA::blockwiseModules(
  make_clr(1.0),
  power = selected_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.10,
  deepSplit = 2,
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  randomSeed = 20260737L,
  nThreads = 4,
  verbose = 0
)
module_colours_pc_one <- WGCNA::labels2colors(
  wgcna_pseudocount_one$colors
)
adjusted_rand <- function(a, b) {
  contingency <- table(a, b)
  choose_two <- function(x) x * (x - 1) / 2
  n <- sum(contingency)
  index <- sum(choose_two(contingency))
  row_index <- sum(choose_two(rowSums(contingency)))
  column_index <- sum(choose_two(colSums(contingency)))
  expected <- row_index * column_index / choose_two(n)
  (index - expected) / (0.5 * (row_index + column_index) - expected)
}
pseudocount_ari <- adjusted_rand(
  module_colours, module_colours_pc_one
)
pair_coclustering_agreement <- mean(
  outer(module_colours, module_colours, "==") ==
    outer(module_colours_pc_one, module_colours_pc_one, "==")
)

data.frame(
  CandidateHubs = sum(hub_ledger$CandidateHub),
  Pseudocount05Modules = length(setdiff(unique(module_colours), "grey")),
  Pseudocount10Modules = length(setdiff(unique(module_colours_pc_one), "grey")),
  AdjustedRand = pseudocount_ari,
  PairCoclusteringAgreement = pair_coclustering_agreement
)

wgcna_audit <- data.frame(
  Dimension = c(
    "Composition", "Taxonomic coverage", "Soft power", "Grey module",
    "Trait family", "Confounding", "Hub", "Sensitivity"
  ),
  FragilePractice = c(
    "WGCNA on percentages",
    "Pretend genus table covers all reads",
    "Pick highest scale-free R2",
    "Force every taxon into a module",
    "Report nominal p-values",
    "Read module-trait r as mechanism",
    "Rank degree only",
    "Show one zero replacement"
  ),
  UpgradedContract = c(
    "Pre-filtered CLR with declared replacement",
    "Report all-read denominator coverage",
    "R2 plus connectivity gate",
    "Retain and report grey fraction",
    "One BH family across module × trait",
    "Design and collinearity boundary",
    "|kME| plus within-module connectivity",
    "Refit modules and report ARI"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(wgcna_audit, caption = "Microbiome WGCNA audit contract")

r2_panel <- ggplot(
  power_table,
  aes(x = Power, y = SFT.R.sq)
) +
  geom_hline(yintercept = 0.80, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_line(colour = pal_pub[["blue"]], linewidth = 0.7) +
  geom_point(
    aes(fill = Power == selected_power),
    shape = 21, size = 2.8, colour = "white", stroke = 0.35
  ) +
  scale_fill_manual(values = c("FALSE" = pal_pub[["grey"]], "TRUE" = pal_pub[["vermillion"]]), guide = "none") +
  labs(
    title = "A. Scale-free topology fit",
    x = "Soft-threshold power", y = "Signed scale-free R²"
  ) +
  theme_pub(base_size = 8.2)

connectivity_panel <- ggplot(
  power_table,
  aes(x = Power, y = mean.k.)
) +
  geom_hline(yintercept = 5, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_line(colour = pal_pub[["orange"]], linewidth = 0.7) +
  geom_point(
    aes(fill = Power == selected_power),
    shape = 21, size = 2.8, colour = "white", stroke = 0.35
  ) +
  scale_fill_manual(values = c("FALSE" = pal_pub[["grey"]], "TRUE" = pal_pub[["vermillion"]]), guide = "none") +
  scale_y_log10() +
  labs(
    title = "B. Mean connectivity",
    x = "Soft-threshold power", y = "Mean connectivity (log scale)"
  ) +
  theme_pub(base_size = 8.2)

soft_threshold_plot <- patchwork::wrap_plots(
  r2_panel, connectivity_panel, nrow = 1L
) +
  patchwork::plot_annotation(
    title = sprintf("Power %d is the first value passing both gates", selected_power),
    subtitle = "Selection requires signed R² >= 0.80 and mean connectivity >= 5",
    caption = "A high R² alone is not sufficient when the resulting network is nearly disconnected."
  )
save_pub(
  soft_threshold_plot, "figures/37-soft-threshold",
  width = 180, height = 122
)
soft_threshold_plot

gene_tree <- wgcna_fit$dendrograms[[1L]]
block_genes <- wgcna_fit$blockGenes[[1L]]
dendrogram_data <- ggdendro::dendro_data(
  stats::as.dendrogram(gene_tree), type = "rectangle"
)
ordered_gene_index <- block_genes[gene_tree$order]
ordered_modules <- module_colours[ordered_gene_index]
module_strip <- data.frame(
  X = seq_along(ordered_modules),
  Y = -0.035 * max(dendrogram_data$segments$y),
  Module = ordered_modules,
  stringsAsFactors = FALSE
)
dendrogram_plot <- ggplot() +
  geom_segment(
    data = dendrogram_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#4D4D4D", linewidth = 0.25
  ) +
  geom_tile(
    data = module_strip,
    aes(x = X, y = Y, fill = Module),
    width = 1.05, height = 0.045 * max(dendrogram_data$segments$y)
  ) +
  scale_fill_identity() +
  coord_cartesian(clip = "off") +
  labs(
    title = "Dynamic tree cut separates four non-grey modules",
    subtitle = sprintf(
      "%d genera remain grey; merge cut height = 0.25",
      sum(module_colours == "grey")
    ),
    x = "Genera ordered by signed TOM dissimilarity",
    y = "Clustering height",
    caption = "Module colours are arbitrary identifiers, not biological annotations."
  ) +
  theme_pub(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 6, 12, 6)
  )
save_pub(
  dendrogram_plot, "figures/37-module-dendrogram",
  width = 180, height = 132
)
dendrogram_plot

module_trait_ledger$ModuleLabel <- sub("^ME", "", module_trait_ledger$Module)
module_trait_ledger$Label <- sprintf(
  "%.2f\nq=%s",
  module_trait_ledger$Correlation,
  ifelse(
    module_trait_ledger$QValue < 0.001,
    "<.001", sprintf("%.3f", module_trait_ledger$QValue)
  )
)
module_trait_ledger$TextColour <- ifelse(
  abs(module_trait_ledger$Correlation) >= 0.55, "white", "black"
)
module_trait_plot <- ggplot(
  module_trait_ledger,
  aes(x = Trait, y = ModuleLabel, fill = Correlation)
) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(
    aes(label = Label, colour = TextColour),
    family = font_pub, size = 2.55, lineheight = 0.88
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = pal_pub[["blue"]], mid = "white", high = pal_pub[["vermillion"]],
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(
    title = "Module eigengenes mostly track the wetland gradient",
    subtitle = sprintf(
      "%d of %d module-trait tests have BH q < 0.05",
      sum(module_trait_q < 0.05), length(module_trait_q)
    ),
    x = NULL, y = "Module eigengene", fill = "Bicor r",
    caption = "CW/IW/TW are one-vs-rest indicators; environmental correlations remain confounded with group and geography."
  ) +
  theme_pub(base_size = 8.2) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )
save_pub(
  module_trait_plot, "figures/37-module-trait",
  width = 180, height = 138
)
module_trait_plot

plot_hubs <- hub_ledger[hub_ledger$Module != "grey", ]
label_hubs <- head(
  plot_hubs[
    order(!plot_hubs$CandidateHub, -abs(plot_hubs$KME), -plot_hubs$WithinConnectivity),
  ],
  14L
)
hub_plot <- ggplot(
  plot_hubs,
  aes(
    x = abs(KME), y = WithinConnectivity,
    fill = Module, shape = CandidateHub
  )
) +
  geom_vline(xintercept = 0.80, linetype = 2, colour = "#555555") +
  geom_point(
    size = 3, alpha = 0.82, colour = "white", stroke = 0.4
  ) +
  ggrepel::geom_text_repel(
    data = label_hubs,
    aes(label = DisplayTaxon),
    family = font_pub, size = 2.45, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 24)) +
  labs(
    title = "Hub status requires membership, connectivity and sensitivity",
    subtitle = sprintf(
      "%d candidates; pseudocount 0.5 vs 1.0 module ARI = %.3f",
      sum(hub_ledger$CandidateHub), pseudocount_ari
    ),
    x = "Absolute module membership (|kME|)",
    y = "Within-module connectivity",
    shape = "Candidate hub",
    caption = "Triangles pass |kME| >= 0.80 and the within-module top-decile rule; labels are candidates, not validated keystones."
  ) +
  theme_pub(base_size = 8.3) +
  theme(legend.position = "bottom")
save_pub(
  hub_plot, "figures/37-hub-sensitivity",
  width = 180, height = 138
)
hub_plot

result_dir <- "results/37-microbial-wgcna"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  power_table,
  file.path(result_dir, "soft-threshold.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  module_sizes,
  file.path(result_dir, "module-sizes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
module_trait_ledger_export <- module_trait_ledger
module_trait_ledger_export$Label <- gsub(
  "\n", "; ", module_trait_ledger_export$Label, fixed = TRUE
)
utils::write.table(
  module_trait_ledger_export,
  file.path(result_dir, "module-trait-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  hub_ledger,
  file.path(result_dir, "hub-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expected_bases <- c(
  "figures/37-soft-threshold",
  "figures/37-module-dendrogram",
  "figures/37-module-trait",
  "figures/37-hub-sensitivity"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  as.character(utils::packageVersion("WGCNA")) == "1.74",
  selected_power == 7L,
  power_gate == "PASS",
  length(setdiff(unique(module_colours), "grey")) == 4L,
  nrow(module_trait_ledger) == 32L,
  pseudocount_ari >= 0 && pseudocount_ari <= 1
)
