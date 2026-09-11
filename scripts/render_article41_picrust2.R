#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpicrust2)
  library(patchwork)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_root)
set.seed(20260741)

font_pub <- "DejaVu Sans"
pal_pub <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
  yellow = "#F0E442", grey = "#6B7280"
)
theme_pub <- function(base_size = 10, base_family = font_pub) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.25),
      axis.text = element_text(colour = "#1A1A1A"),
      axis.title = element_text(colour = "#1A1A1A"),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(colour = "#4D4D4D"),
      plot.caption = element_text(colour = "#666666", hjust = 0),
      strip.background = element_rect(fill = "#F2F2F2", colour = "#B3B3B3"),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(), legend.position = "top"
    )
}
save_pub <- function(plot, file_base, width = 89, height = 70, units = "mm", dpi = 600) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = units, device = grDevices::cairo_pdf, family = font_pub, bg = "white")
  ggsave(paste0(file_base, ".svg"), plot, width = width, height = height,
         units = units, device = svglite::svglite, bg = "white")
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = units, dpi = dpi, device = ragg::agg_png, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = units, dpi = dpi, device = ragg::agg_tiff,
         compression = "lzw", bg = "white")
  invisible(plot)
}
read_keyed <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, progress = FALSE,
                name_repair = "minimal", na = character())
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_dir <- "data/small/picrust2-chemerin"
output_dir <- file.path(input_dir, "picrust2-out-v2.6.3")
result_dir <- "results/41-picrust2"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_keyed(file.path(input_dir, "metadata.tsv"))
metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "KO"))
metadata$Facility <- factor(metadata$Facility, levels = c("Dalhousie", "PaloAlto"))
counts <- as.matrix(data.frame(
  lapply(read_keyed(file.path(input_dir, "otutab.tsv")), as.integer),
  check.names = FALSE
))
rownames(counts) <- rownames(read_keyed(file.path(input_dir, "otutab.tsv")))
counts <- counts[, rownames(metadata), drop = FALSE]

nsti <- read_keyed(file.path(output_dir, "combined_marker_predicted_and_nsti.tsv.gz"))
nsti$metadata_NSTI <- as.numeric(nsti$metadata_NSTI)
weighted_nsti <- read_keyed(file.path(output_dir, "EC_metagenome_out", "weighted_nsti.tsv.gz"))
weighted_nsti$weighted_NSTI <- as.numeric(weighted_nsti$weighted_NSTI)
weighted_nsti <- weighted_nsti[rownames(metadata), , drop = FALSE]
nsti_ledger <- data.frame(
  SampleID = rownames(metadata), Genotype = metadata$Genotype,
  Facility = metadata$Facility, WeightedNSTI = weighted_nsti$weighted_NSTI,
  stringsAsFactors = FALSE
)

pathway_abundance <- as.matrix(data.frame(
  lapply(read_keyed(file.path(output_dir, "pathways_out", "path_abun_unstrat.tsv.gz")), as.numeric),
  check.names = FALSE
))
rownames(pathway_abundance) <- rownames(
  read_keyed(file.path(output_dir, "pathways_out", "path_abun_unstrat.tsv.gz"))
)
pathway_abundance <- pathway_abundance[, rownames(metadata), drop = FALSE]
pathway_coverage <- as.matrix(data.frame(
  lapply(read_keyed(file.path(output_dir, "pathways_out", "path_cov_unstrat.tsv.gz")), as.numeric),
  check.names = FALSE
))
rownames(pathway_coverage) <- rownames(
  read_keyed(file.path(output_dir, "pathways_out", "path_cov_unstrat.tsv.gz"))
)
pathway_coverage <- pathway_coverage[, rownames(metadata), drop = FALSE]

annotation <- ggpicrust2::metacyc_reference
description <- setNames(as.character(annotation$description), annotation$id)
pathway_ids <- rownames(pathway_abundance)
pathway_description <- unname(description[pathway_ids])
pathway_description[is.na(pathway_description) | pathway_description == ""] <-
  pathway_ids[is.na(pathway_description) | pathway_description == ""]
names(pathway_description) <- pathway_ids

prevalence <- rowMeans(pathway_abundance > 0)
eligible <- prevalence >= 0.20 & rowSums(pathway_abundance) > 0
analysis_abundance <- pathway_abundance[eligible, , drop = FALSE]
clr <- log(analysis_abundance + 0.5)
clr <- sweep(clr, 2L, colMeans(clr), "-")

test_pathway <- function(pathway) {
  model_data <- data.frame(
    Value = clr[pathway, ], Genotype = metadata$Genotype,
    Facility = metadata$Facility
  )
  fit <- lm(Value ~ Genotype + Facility, data = model_data)
  coefficient <- summary(fit)$coefficients["GenotypeKO", ]
  coverage_values <- if (pathway %in% rownames(pathway_coverage)) {
    pathway_coverage[pathway, ]
  } else {
    rep(NA_real_, nrow(metadata))
  }
  data.frame(
    Pathway = pathway,
    Description = unname(pathway_description[pathway]),
    KOminusWT = unname(coefficient["Estimate"]),
    SE = unname(coefficient["Std. Error"]),
    PValue = unname(coefficient["Pr(>|t|)"]),
    Prevalence = prevalence[[pathway]],
    MeanRelativeAbundance = mean(
      sweep(pathway_abundance, 2L, colSums(pathway_abundance), "/")[pathway, ]
    ),
    MedianCoverage = stats::median(coverage_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
pathway_tests <- do.call(rbind, lapply(rownames(analysis_abundance), test_pathway))
pathway_tests$AdjustedP <- p.adjust(pathway_tests$PValue, method = "BH")
pathway_tests <- pathway_tests[order(pathway_tests$AdjustedP, -abs(pathway_tests$KOminusWT)), ]

coverage_summary <- data.frame(
  Pathway = rownames(pathway_coverage),
  Description = unname(description[rownames(pathway_coverage)]),
  MedianCoverage = apply(pathway_coverage, 1L, median, na.rm = TRUE),
  MinimumCoverage = apply(pathway_coverage, 1L, min, na.rm = TRUE),
  MaximumCoverage = apply(pathway_coverage, 1L, max, na.rm = TRUE),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(metadata) == 24L, nrow(counts) == 37L, sum(counts) == 108718L,
  nrow(nsti) == 37L, all(nsti$best_domain == "bac"),
  max(nsti$metadata_NSTI) < 2,
  ncol(pathway_abundance) == 24L,
  nrow(pathway_tests) > 100L,
  all(is.finite(pathway_tests$KOminusWT))
)

# Figure 1: analysis workflow.
workflow <- data.frame(
  x = 1:6, y = 1,
  Label = c("ASV counts", "Representative\nsequences", "EPA-ng\nplacement",
            "Hidden-state\nprediction", "EC / KO\nprofiles", "MetaCyc\npathways"),
  Detail = c("37 × 24", "400 bp", "Bacteria + Archaea", "NSTI + 16S copy",
             "Copy-normalized", "Abundance + coverage"),
  stringsAsFactors = FALSE
)
workflow_edges <- data.frame(x = 1:5 + 0.37, xend = 2:6 - 0.37, y = 1, yend = 1)
p_workflow <- ggplot() +
  geom_segment(data = workflow_edges, aes(x, y, xend = xend, yend = yend),
               linewidth = 0.65, colour = pal_pub[["grey"]],
               arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  geom_label(data = workflow, aes(x, y, label = Label),
             fill = c("#EAF4FB", "#EAF4FB", "#FFF4DC", "#FFF4DC", "#E9F7F1", "#E9F7F1"),
             colour = "#1A1A1A", label.size = 0.3, size = 3.1, lineheight = 0.95) +
  geom_text(data = workflow, aes(x, y - 0.25, label = Detail), size = 2.5, colour = "#4D4D4D") +
  annotate("text", x = 3.5, y = 1.43, label = "PICRUSt2 2.6.3", fontface = "bold", size = 4) +
  coord_cartesian(xlim = c(0.55, 6.45), ylim = c(0.55, 1.55), clip = "off") +
  labs(caption = "Predictions are functional potential inferred from phylogenetic placement, not measured genes or activity.") +
  theme_void(base_family = font_pub) +
  theme(plot.caption = element_text(colour = "#666666", hjust = 0, margin = margin(t = 8)))
save_pub(p_workflow, "figures/41-picrust2-workflow", 178, 68)

# Figure 2: weighted and feature-level NSTI audit.
p_weighted <- ggplot(nsti_ledger, aes(Genotype, WeightedNSTI, colour = Facility)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "#4D4D4D") +
  geom_jitter(width = 0.10, height = 0, size = 1.7, alpha = 0.75) +
  scale_colour_manual(values = c(Dalhousie = pal_pub[["blue"]], PaloAlto = pal_pub[["orange"]])) +
  labs(title = "Sample-weighted NSTI", x = NULL, y = "Weighted NSTI", colour = "Facility") +
  theme_pub()
nsti_feature <- data.frame(ASV = rownames(nsti), NSTI = nsti$metadata_NSTI)
p_feature <- ggplot(nsti_feature, aes(NSTI)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = pal_pub[["green"]], colour = "white") +
  geom_rug(colour = pal_pub[["green"]], alpha = 0.65) +
  labs(title = "ASV-level NSTI", x = "ASV NSTI", y = "ASVs",
       caption = sprintf("37/37 ASVs retained; max %.3f (cutoff 2.0).", max(nsti$metadata_NSTI))) +
  theme_pub()
p_nsti <- p_weighted + p_feature + patchwork::plot_layout(widths = c(1, 1))
save_pub(p_nsti, "figures/41-nsti-audit", 178, 88)

# Figure 3: facility-aware pathway heatmap.
selected <- head(pathway_tests$Pathway, 12L)
group_id <- interaction(metadata$Genotype, metadata$Facility, sep = " · ")
heat <- do.call(rbind, lapply(selected, function(pathway) {
  means <- tapply(clr[pathway, ], group_id, mean)
  data.frame(
    Pathway = pathway,
    Label = paste0(substr(pathway_description[pathway], 1, 46), " [", pathway, "]"),
    Group = names(means), Value = as.numeric(scale(means)), stringsAsFactors = FALSE
  )
}))
heat$Label <- factor(heat$Label, levels = rev(unique(heat$Label)))
p_heat <- ggplot(heat, aes(Group, Label, fill = Value)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = pal_pub[["blue"]], mid = "white",
                       high = pal_pub[["vermillion"]], midpoint = 0) +
  labs(title = "Predicted MetaCyc pathway profiles",
       subtitle = "Rows selected by facility-adjusted genotype test",
       x = NULL, y = NULL, fill = "Row z-score",
       caption = "MetaCyc labels: ggpicrust2 2.5.17.") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_pub(p_heat, "figures/41-pathway-heatmap", 178, 112)

# Figure 4: effect and multiplicity audit.
label_effect <- head(pathway_tests[order(pathway_tests$AdjustedP, -abs(pathway_tests$KOminusWT)), ], 8L)
p_effect <- ggplot(pathway_tests, aes(KOminusWT, -log10(pmax(AdjustedP, 1e-12)))) +
  geom_vline(xintercept = 0, linetype = 3, colour = "#B3B3B3") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "#B3B3B3") +
  geom_point(aes(colour = MedianCoverage, size = MeanRelativeAbundance), alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = label_effect, aes(label = Pathway), size = 2.6,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  scale_colour_gradient(low = pal_pub[["orange"]], high = pal_pub[["blue"]], na.value = pal_pub[["grey"]]) +
  scale_size_continuous(labels = scales::label_percent(accuracy = 0.01)) +
  labs(title = "Facility-adjusted genotype effects",
       x = "KO − WT CLR coefficient", y = "−log10 BH-adjusted p",
       colour = "Median coverage", size = "Mean abundance",
       caption = sprintf("BH correction covers all %d prevalence-filtered pathways.", nrow(pathway_tests))) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical")
save_pub(p_effect, "figures/41-pathway-effect", 178, 105)

write_tsv(nsti_ledger, file.path(result_dir, "sample-weighted-nsti.tsv"))
write_tsv(data.frame(
  ASV = rownames(nsti), NSTI = nsti$metadata_NSTI,
  Domain = nsti$best_domain,
  ClosestReferenceGenome = nsti$closest_reference_genome
), file.path(result_dir, "asv-nsti.tsv"))
write_tsv(pathway_tests, file.path(result_dir, "pathway-tests.tsv"))
write_tsv(coverage_summary, file.path(result_dir, "pathway-coverage-summary.tsv"))
write_tsv(data.frame(
  PICRUSt2 = "2.6.3", ggpicrust2 = as.character(packageVersion("ggpicrust2")),
  Samples = nrow(metadata), ASVs = nrow(nsti), BacterialASVs = sum(nsti$best_domain == "bac"),
  ArchaealASVs = sum(nsti$best_domain == "arc"), MaxNSTI = max(nsti$metadata_NSTI),
  MedianWeightedNSTI = median(nsti_ledger$WeightedNSTI),
  PredictedPathways = nrow(pathway_abundance), TestedPathways = nrow(pathway_tests)
), file.path(result_dir, "run-summary.tsv"))

message("Article 41 outputs and four publication figures completed.")
