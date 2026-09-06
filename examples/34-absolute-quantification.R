# 绝对定量陷阱：相对丰度 vs QMP/spike-in
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, knitr, ragg, scales, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/absolute-quantification/cell-load.tsv",
  "data/small/absolute-quantification/metadata.tsv",
  "data/small/absolute-quantification/otutab.tsv",
  "data/small/absolute-quantification/source-summary.json",
  "data/small/absolute-quantification/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)

set.seed(20260734)
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
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.15)
      ),
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
data_dir <- "data/small/absolute-quantification"
otutab_raw <- read_keyed_tsv(file.path(data_dir, "otutab.tsv"))
taxonomy <- read_keyed_tsv(file.path(data_dir, "taxonomy.tsv"))
metadata_all <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
cell_load_all <- read_keyed_tsv(file.path(data_dir, "cell-load.tsv"))
source_summary <- jsonlite::read_json(file.path(data_dir, "source-summary.json"))

otutab_all <- as.matrix(data.frame(
  lapply(otutab_raw, as.integer), check.names = FALSE
))
rownames(otutab_all) <- rownames(otutab_raw)
cell_load_all$AbsoluteCellLoad <- as.numeric(cell_load_all$AbsoluteCellLoad)

stopifnot(
  identical(rownames(otutab_all), rownames(taxonomy)),
  identical(colnames(otutab_all), rownames(metadata_all)),
  identical(colnames(otutab_all), rownames(cell_load_all)),
  nrow(otutab_all) == 234L,
  ncol(otutab_all) == 135L,
  sum(otutab_all) == 4080996L,
  all(otutab_all >= 0L),
  all(otutab_all == round(otutab_all)),
  all(cell_load_all$AbsoluteCellLoad > 0),
  identical(unique(cell_load_all$Unit), "cells per gram frozen feces"),
  identical(source_summary$package$version, "1.16.0")
)

head(otutab_all[, 1:5])
head(taxonomy)
head(metadata_all)
head(cell_load_all)

sample_ids <- rownames(metadata_all)[metadata_all$Cohort == "Disease"]
metadata <- metadata_all[sample_ids, , drop = FALSE]
metadata$HealthStatus <- relevel(
  factor(metadata$HealthStatus), ref = "Healthy"
)
cell_load <- cell_load_all[sample_ids, , drop = FALSE]
counts <- otutab_all[, sample_ids, drop = FALSE]

stopifnot(
  length(sample_ids) == 95L,
  sum(metadata$HealthStatus == "CD") == 29L,
  sum(metadata$HealthStatus == "Healthy") == 66L,
  identical(colnames(counts), rownames(metadata)),
  identical(colnames(counts), rownames(cell_load))
)

table(metadata$Cohort, metadata$HealthStatus)
summary(cell_load$AbsoluteCellLoad)

library_size <- colSums(counts)
relative_abundance <- sweep(counts, 2L, library_size, "/")
quantitative_abundance <- sweep(
  relative_abundance,
  2L,
  cell_load$AbsoluteCellLoad,
  "*"
)

stopifnot(
  max(abs(colSums(relative_abundance) - 1)) < 1e-12,
  max(abs(colSums(quantitative_abundance) - cell_load$AbsoluteCellLoad)) < 1e-3
)

load_cd <- cell_load$AbsoluteCellLoad[metadata$HealthStatus == "CD"]
load_healthy <- cell_load$AbsoluteCellLoad[metadata$HealthStatus == "Healthy"]
load_wilcox <- stats::wilcox.test(load_cd, load_healthy, exact = FALSE)
load_log2_median_ratio <- log2(median(load_cd) / median(load_healthy))

set.seed(20260734)
bootstrap_replicates <- 2000L
load_bootstrap <- replicate(bootstrap_replicates, {
  log2(
    median(sample(load_cd, replace = TRUE)) /
      median(sample(load_healthy, replace = TRUE))
  )
})
load_bootstrap_ci <- stats::quantile(
  load_bootstrap, probs = c(0.025, 0.975), names = FALSE
)
load_cliffs_delta <- mean(outer(load_cd, load_healthy, ">")) -
  mean(outer(load_cd, load_healthy, "<"))

load_summary <- data.frame(
  Contrast = "CD versus Healthy (Disease cohort)",
  CDMedianCellsPerGram = median(load_cd),
  HealthyMedianCellsPerGram = median(load_healthy),
  Log2MedianRatio = load_log2_median_ratio,
  BootstrapCI95Lower = load_bootstrap_ci[[1L]],
  BootstrapCI95Upper = load_bootstrap_ci[[2L]],
  CliffsDelta = load_cliffs_delta,
  WilcoxonP = load_wilcox$p.value,
  stringsAsFactors = FALSE
)
load_summary

prevalence <- rowMeans(counts > 0)
eligible <- prevalence >= 0.10 & rowSums(counts) >= 20
feature_ids <- rownames(counts)[eligible]
group <- metadata$HealthStatus

safe_log2_median_effect <- function(x, group) {
  positive <- x[x > 0 & is.finite(x)]
  pseudocount <- if (length(positive)) min(positive) / 2 else 0.5
  median(log2(x[group == "CD"] + pseudocount)) -
    median(log2(x[group == "Healthy"] + pseudocount))
}
safe_wilcox <- function(x, group) {
  if (length(unique(x)) < 2L) return(1)
  stats::wilcox.test(
    x[group == "CD"], x[group == "Healthy"], exact = FALSE
  )$p.value
}

relative_effect <- vapply(feature_ids, function(id) {
  safe_log2_median_effect(relative_abundance[id, ], group)
}, numeric(1))
quantitative_effect <- vapply(feature_ids, function(id) {
  safe_log2_median_effect(quantitative_abundance[id, ], group)
}, numeric(1))
relative_p <- vapply(feature_ids, function(id) {
  safe_wilcox(relative_abundance[id, ], group)
}, numeric(1))
quantitative_p <- vapply(feature_ids, function(id) {
  safe_wilcox(quantitative_abundance[id, ], group)
}, numeric(1))

scale_comparison <- data.frame(
  FeatureID = feature_ids,
  Genus = taxonomy[feature_ids, "Genus"],
  Prevalence = prevalence[feature_ids],
  RelativeEffect = relative_effect,
  QuantitativeEffect = quantitative_effect,
  RelativeP = relative_p,
  QuantitativeP = quantitative_p,
  stringsAsFactors = FALSE
)
scale_comparison$ReportableGenus <- !grepl(
  "unclassified|unknown|unassigned|uncultured|metagenome|^other$",
  scale_comparison$Genus,
  ignore.case = TRUE
) & nzchar(trimws(scale_comparison$Genus))
scale_comparison$RelativeQ <- stats::p.adjust(scale_comparison$RelativeP, method = "BH")
scale_comparison$QuantitativeQ <- stats::p.adjust(
  scale_comparison$QuantitativeP, method = "BH"
)
scale_comparison$DirectionConcordant <- sign(scale_comparison$RelativeEffect) ==
  sign(scale_comparison$QuantitativeEffect)
scale_comparison$EvidenceClass <- with(
  scale_comparison,
  ifelse(
    RelativeQ < 0.05 & QuantitativeQ < 0.05, "Both scales",
    ifelse(
      RelativeQ < 0.05, "Relative only",
      ifelse(QuantitativeQ < 0.05, "Quantitative only", "Neither")
    )
  )
)
scale_comparison$EffectShift <- scale_comparison$QuantitativeEffect -
  scale_comparison$RelativeEffect
scale_comparison <- scale_comparison[order(
  -abs(scale_comparison$EffectShift),
  pmin(scale_comparison$RelativeQ, scale_comparison$QuantitativeQ)
), ]

stopifnot(
  nrow(scale_comparison) >= 50L,
  all(is.finite(scale_comparison$RelativeEffect)),
  all(is.finite(scale_comparison$QuantitativeEffect)),
  all(scale_comparison$RelativeQ >= 0 & scale_comparison$RelativeQ <= 1),
  all(scale_comparison$QuantitativeQ >= 0 & scale_comparison$QuantitativeQ <= 1),
  sum(scale_comparison$ReportableGenus) >= 20L
)

table(scale_comparison$EvidenceClass)
table(scale_comparison$DirectionConcordant)
head(scale_comparison, 15)

exemplar_pool <- scale_comparison[
  scale_comparison$ReportableGenus & (
    scale_comparison$RelativeQ < 0.05 |
      scale_comparison$QuantitativeQ < 0.05
  ),
]
stopifnot(nrow(exemplar_pool) > 0L)
exemplar_pool$SelectionPriority <- ifelse(
  !exemplar_pool$DirectionConcordant,
  2L,
  ifelse(exemplar_pool$EvidenceClass != "Both scales", 1L, 0L)
)
exemplar_pool <- exemplar_pool[order(
  -exemplar_pool$SelectionPriority,
  -abs(exemplar_pool$EffectShift),
  exemplar_pool$Genus
), ]
exemplar <- exemplar_pool[1L, ]
exemplar_id <- exemplar$FeatureID[[1L]]
exemplar_genus <- exemplar$Genus[[1L]]

exemplar_long <- rbind(
  data.frame(
    SampleID = sample_ids,
    HealthStatus = metadata$HealthStatus,
    Scale = "Relative abundance (%)",
    Value = 100 * relative_abundance[exemplar_id, ],
    stringsAsFactors = FALSE
  ),
  data.frame(
    SampleID = sample_ids,
    HealthStatus = metadata$HealthStatus,
    Scale = "Quantitative abundance (cells/g)",
    Value = quantitative_abundance[exemplar_id, ],
    stringsAsFactors = FALSE
  )
)
exemplar

absolute_audit <- data.frame(
  Dimension = c(
    "Sample denominator", "Load measurement", "QMP formula",
    "Cohort", "Taxon bias", "Claim"
  ),
  FragilePractice = c(
    "Use reads as cells",
    "Omit units and calibration",
    "Multiply after unrelated normalization",
    "Mix controls from different cohorts",
    "Assume QMP removes every bias",
    "Call association an absolute causal change"
  ),
  UpgradedContract = c(
    "Use sample mass/volume and an independent scale",
    "Report platform, gating/curve, QC and units",
    "Multiply within-sample proportions by matched load",
    "Primary analysis within a shared collection cohort",
    "Retain extraction/PCR/copy-number limitations",
    "Report calibrated association plus orthogonal validation"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(absolute_audit, caption = "Quantitative microbiome profiling audit")

load_plot_data <- data.frame(
  SampleID = sample_ids,
  HealthStatus = metadata$HealthStatus,
  CellLoad = cell_load$AbsoluteCellLoad,
  stringsAsFactors = FALSE
)
load_plot <- ggplot(
  load_plot_data,
  aes(x = HealthStatus, y = CellLoad, fill = HealthStatus)
) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.68) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 20260734),
    shape = 21, colour = "white", stroke = 0.25,
    size = 1.65, alpha = 0.78
  ) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_fill_manual(values = c(
    Healthy = pal_pub[["blue"]], CD = pal_pub[["vermillion"]]
  )) +
  labs(
    title = "Crohn's disease samples have a different microbial-load scale",
    subtitle = "Primary comparison restricted to the shared Disease cohort",
    x = NULL, y = "Flow-cytometry cell load (cells/g frozen feces)",
    fill = "Health status",
    caption = sprintf(
      "Median log2 ratio (CD / Healthy) = %.2f; bootstrap 95%% CI %.2f to %.2f; Wilcoxon p = %.2g.",
      load_log2_median_ratio, load_bootstrap_ci[[1L]],
      load_bootstrap_ci[[2L]], load_wilcox$p.value
    )
  ) +
  theme_pub(base_size = 8.5) +
  theme(legend.position = "none")

save_pub(load_plot, "figures/34-microbial-load", width = 152, height = 122)
load_plot

effect_plot_data <- scale_comparison
label_data <- head(
  effect_plot_data[effect_plot_data$ReportableGenus, ],
  12L
)
effect_plot <- ggplot(
  effect_plot_data,
  aes(
    x = RelativeEffect, y = QuantitativeEffect,
    colour = EvidenceClass, shape = DirectionConcordant
  )
) +
  geom_hline(yintercept = 0, colour = "#888888", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#777777") +
  geom_point(size = 2.1, alpha = 0.80, stroke = 0.65) +
  ggrepel::geom_text_repel(
    data = label_data, aes(label = Genus),
    family = font_pub, size = 2.4, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Both scales" = pal_pub[["purple"]],
    "Relative only" = pal_pub[["orange"]],
    "Quantitative only" = pal_pub[["green"]],
    "Neither" = "#B8B8B8"
  )) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "Same direction", `FALSE` = "Opposite direction")
  ) +
  labs(
    title = "Microbial load changes taxon-level effect estimates",
    subtitle = "Feature-wise median log2 effects: CD versus Healthy",
    x = "Relative-abundance effect",
    y = "Quantitative-abundance effect",
    colour = "BH evidence", shape = "Direction audit",
    caption = "Dashed identity line marks no effect shift after restoring sample load."
  ) +
  theme_pub(base_size = 8.3) +
  theme(legend.position = "bottom", legend.box = "vertical")

save_pub(
  effect_plot, "figures/34-relative-quantitative-effects",
  width = 170, height = 132
)
effect_plot

exemplar_plot <- ggplot(
  exemplar_long,
  aes(x = HealthStatus, y = Value, fill = HealthStatus)
) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.68) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0, seed = 20260734),
    shape = 21, colour = "white", stroke = 0.25,
    size = 1.45, alpha = 0.75
  ) +
  facet_wrap(~Scale, scales = "free_y", nrow = 1L) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_fill_manual(values = c(
    Healthy = pal_pub[["blue"]], CD = pal_pub[["vermillion"]]
  )) +
  labs(
    title = paste0(exemplar_genus, ": the denominator changes the story"),
    subtitle = sprintf(
      "Relative effect %.2f (q %.2g); quantitative effect %.2f (q %.2g)",
      exemplar$RelativeEffect, exemplar$RelativeQ,
      exemplar$QuantitativeEffect, exemplar$QuantitativeQ
    ),
    x = NULL, y = "Observed value", fill = "Health status",
    caption = "Quantitative values equal within-sample sequencing proportions × matched flow-cytometry load."
  ) +
  theme_pub(base_size = 8.2) +
  theme(legend.position = "none")

save_pub(exemplar_plot, "figures/34-qmp-exemplar", width = 180, height = 112)
exemplar_plot

expected_bases <- c(
  "figures/34-microbial-load",
  "figures/34-relative-quantitative-effects",
  "figures/34-qmp-exemplar"
)
expected_files <- as.vector(outer(
  expected_bases, c(".pdf", ".svg", ".png", ".tiff"), paste0
))
stopifnot(
  all(file.exists(expected_files)),
  bootstrap_replicates == 2000L,
  load_log2_median_ratio < 0,
  nrow(exemplar) == 1L,
  exemplar_id %in% rownames(counts),
  any(!scale_comparison$DirectionConcordant)
)
