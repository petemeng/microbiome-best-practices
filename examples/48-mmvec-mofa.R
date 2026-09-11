# MMvec/MOFA（菌×代谢物）
# Run sequentially in a new working directory.
# Required packages: ggplot2, jsonlite, patchwork, ragg, readr, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/run_article48_native_mmvec.py",
  "scripts/run_article48_mmvec_mofa.py",
  "env/mmvec-native.yml",
  "env/multiomics.yml",
  "data/small/paired-ibd-multiomics/metabolite-annotation.tsv",
  "data/small/paired-ibd-multiomics/metabolites.tsv",
  "data/small/paired-ibd-multiomics/metadata.tsv",
  "data/small/paired-ibd-multiomics/otutab.tsv",
  "data/small/paired-ibd-multiomics/source-summary.json",
  "data/small/paired-ibd-multiomics/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(patchwork)

set.seed(20260748)
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
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.25),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#B3B3B3"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(), legend.position = "top"
    )
}
save_pub <- function(plot, file_base, width = 89, height = 70, units = "mm", dpi = 600) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
                  units = units, device = grDevices::cairo_pdf, family = font_pub, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".svg"), plot, width = width, height = height,
                  units = units, device = svglite::svglite, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
                  units = units, dpi = dpi, device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
                  units = units, dpi = dpi, device = ragg::agg_tiff,
                  compression = "lzw", bg = "white")
  invisible(plot)
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                       name_repair = "minimal", na = character())
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}
data_dir <- "data/small/paired-ibd-multiomics"
otutab_raw <- read_keyed_tsv(file.path(data_dir, "otutab.tsv"))
taxonomy <- read_keyed_tsv(file.path(data_dir, "taxonomy.tsv"))
metadata <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
metabolites_raw <- read_keyed_tsv(file.path(data_dir, "metabolites.tsv"))
metabolite_annotation <- read_keyed_tsv(
  file.path(data_dir, "metabolite-annotation.tsv")
)
otutab <- as.matrix(data.frame(lapply(otutab_raw, as.numeric), check.names = FALSE))
metabolites <- as.matrix(data.frame(
  lapply(metabolites_raw, as.numeric), check.names = FALSE
))
rownames(otutab) <- rownames(otutab_raw)
rownames(metabolites) <- rownames(metabolites_raw)
stopifnot(
  nrow(otutab) == 250L, ncol(otutab) == 220L,
  nrow(metabolites) == 277L, ncol(metabolites) == 220L,
  identical(colnames(otutab), rownames(metadata)),
  identical(colnames(metabolites), rownames(metadata)),
  all(otutab >= 0), all(metabolites >= 0)
)
head(otutab[, 1:5])
head(taxonomy)
head(metadata)
head(metabolites[, 1:5])

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(Sys.getenv("MMVEC_PYTHON", "mmvec-env/bin/python"), c("scripts/run_article48_native_mmvec.py"))
run_analysis(Sys.getenv("MULTIOMICS_PYTHON", "multiomics-env/bin/python"), c("scripts/run_article48_mmvec_mofa.py"))

microbe_features <- names(sort(rowSums(otutab), decreasing = TRUE))[1:100]
metabolite_variance <- apply(log1p(t(metabolites)), 2, var)
metabolite_features <- names(sort(metabolite_variance, decreasing = TRUE))[1:120]
analysis_microbes <- otutab[microbe_features, , drop = FALSE]
analysis_metabolites <- metabolites[metabolite_features, , drop = FALSE]
stopifnot(
  dim(analysis_microbes)[1] == 100L,
  dim(analysis_metabolites)[1] == 120L,
  identical(colnames(analysis_microbes), colnames(analysis_metabolites))
)

result_dir <- "results/48-mmvec-mofa"
cca_scores <- readr::read_tsv(
  file.path(result_dir, "cca-scores.tsv"), show_col_types = FALSE
)
mmvec_loss <- readr::read_tsv(
  file.path(result_dir, "mmvec-loss.tsv"), show_col_types = FALSE
)
mmvec_pairs <- readr::read_tsv(
  file.path(result_dir, "mmvec-top-pairs.tsv"), show_col_types = FALSE
)
mofa_variance <- readr::read_tsv(
  file.path(result_dir, "mofa-variance-explained.tsv"), show_col_types = FALSE
)
mofa_scores <- readr::read_tsv(
  file.path(result_dir, "mofa-factor-scores.tsv"), show_col_types = FALSE
)
mofa_tests <- readr::read_tsv(
  file.path(result_dir, "mofa-factor-group-tests.tsv"), show_col_types = FALSE
)
integration_summary <- jsonlite::read_json(
  file.path(result_dir, "summary.json"), simplifyVector = TRUE
)
stopifnot(
  integration_summary$spearman_tested_pairs == 12000L,
  integration_summary$spearman_fdr_hits == 4166L,
  integration_summary$mmvec_train_samples == 153L,
  integration_summary$mmvec_holdout_samples == 67L,
  integration_summary$mmvec_epochs_run == 100L,
  integration_summary$mmvec_best_epoch_one_based == 100L,
  integration_summary$mofa_factors_retained == 4L,
  integration_summary$spearman_mmvec_top50_pair_overlap == 0L
)
integration_summary[c(
  "cca_canonical_correlations", "mmvec_best_validation_mae",
  "mofa_total_variance_explained_by_view", "mofa_group_associated_factors_fdr"
)]

cca_long <- do.call(rbind, lapply(1:3, function(component) {
  data.frame(
    StudyGroup = cca_scores$StudyGroup,
    Component = paste0("CC", component),
    Microbiome = cca_scores[[paste0("MicrobiomeCC", component)]],
    Metabolome = cca_scores[[paste0("MetabolomeCC", component)]]
  )
}))
cca_cor <- aggregate(cbind(Microbiome, Metabolome) ~ Component,
                     data = cca_long, FUN = mean)
cca_labels <- setNames(
  sprintf("r = %.4f", integration_summary$cca_canonical_correlations),
  paste0("CC", 1:3)
)
group_colors <- c(
  Control = unname(pal_pub["blue"]),
  CD = unname(pal_pub["vermillion"]),
  UC = unname(pal_pub["orange"])
)
p48_1 <- ggplot(cca_long, aes(Microbiome, Metabolome, colour = StudyGroup)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777", linewidth = 0.4) +
  geom_point(size = 1.2, alpha = 0.65) +
  facet_wrap(~ Component, scales = "free") +
  geom_text(
    data = data.frame(Component = paste0("CC", 1:3), x = -Inf, y = Inf,
                      label = unname(cca_labels)),
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.3, family = font_pub, size = 3
  ) +
  scale_color_manual(values = group_colors) +
  labs(x = "Microbiome canonical score", y = "Metabolome canonical score",
       colour = "Study group", title = "Naive CCA: an overfitting warning",
       subtitle = "In-sample correlations require held-out validation") + theme_pub()
save_pub(p48_1, "figures/48-mmvec-mofa/48-1-cca-overfit-audit", 180, 78)
p48_1

loss_long <- rbind(
  data.frame(Epoch = mmvec_loss$Epoch,
             Metric = "Training objective", Value = mmvec_loss$TrainLoss),
  data.frame(Epoch = mmvec_loss$Epoch,
             Metric = "Holdout MAE", Value = mmvec_loss$ValidationMAE)
)
loss_long$Metric <- factor(
  loss_long$Metric, levels = c("Training objective", "Holdout MAE")
)
best_epoch <- integration_summary$mmvec_best_epoch_one_based
best_labels <- data.frame(
  Epoch = best_epoch,
  Value = Inf,
  Metric = factor("Holdout MAE",
                  levels = levels(loss_long$Metric)),
  Label = paste0("Best: ", best_epoch)
)
p48_2 <- ggplot(loss_long, aes(Epoch, Value)) +
  geom_line(linewidth = 0.7, colour = pal_pub["blue"]) +
  geom_vline(xintercept = best_epoch, linetype = 2, colour = "#555555") +
  geom_text(data = best_labels, aes(Epoch, Value, label = Label),
            inherit.aes = FALSE, hjust = 1.04, vjust = 1.25,
            family = font_pub, size = 2.7) +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
  labs(x = "Epoch", y = "Metric value",
       title = "Native MMvec optimization",
       subtitle = "Best holdout value occurs at the epoch limit") +
  theme_pub(9)
save_pub(p48_2, "figures/48-mmvec-mofa/48-2-mmvec-loss", 89, 88)
p48_2

top_mmvec <- head(mmvec_pairs, 15)
top_mmvec$Pair <- paste0(top_mmvec$MicrobeName, "  |  ", top_mmvec$MetaboliteName)
top_mmvec$Pair <- factor(top_mmvec$Pair, levels = rev(top_mmvec$Pair))
p48_3 <- ggplot(top_mmvec, aes(score, Pair)) +
  geom_segment(aes(x = 0, xend = score, yend = Pair), colour = pal_pub["sky"], linewidth = 0.7) +
  geom_point(colour = pal_pub["blue"], size = 2.8) +
  labs(x = "MMvec conditional co-occurrence score", y = NULL,
       title = "Top MMvec pairs",
       subtitle = "Exact overlap with Spearman top 50: 0; estimands differ") +
  theme_pub(9)
save_pub(p48_3, "figures/48-mmvec-mofa/48-3-mmvec-top-pairs", 180, 105)
p48_3

mofa_variance$Percent <- 100 * mofa_variance$VarianceExplained
p_variance <- ggplot(mofa_variance, aes(Factor, View, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), family = font_pub, size = 3) +
  scale_fill_gradient(low = "#F7FBFF", high = pal_pub["blue"]) +
  labs(x = NULL, y = NULL, fill = "Variance (%)", title = "Variance explained") +
  theme_pub(9) + theme(legend.position = "right")

mofa_long <- reshape(
  as.data.frame(mofa_scores),
  varying = paste0("Factor", 1:4), v.names = "Score",
  timevar = "Factor", times = paste0("Factor", 1:4), direction = "long"
)
mofa_long$Factor <- factor(mofa_long$Factor, levels = paste0("Factor", 1:4))
q_labels <- setNames(sprintf("BH q = %.2g", mofa_tests$QValue), mofa_tests$Factor)
p_scores <- ggplot(mofa_long, aes(StudyGroup, Score, fill = StudyGroup)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.80) +
  geom_jitter(width = 0.14, size = 0.55, alpha = 0.35) +
  facet_wrap(~ Factor, scales = "free_y", labeller = as_labeller(function(x) {
    paste0(x, "\n", q_labels[x])
  })) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Factor score", fill = "Study group",
       title = "External association audit") + theme_pub(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p48_4 <- p_variance / p_scores + plot_layout(heights = c(0.65, 1.35))
save_pub(p48_4, "figures/48-mmvec-mofa/48-4-mofa-summary", 180, 135)
p48_4
