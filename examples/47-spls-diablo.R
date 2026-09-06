# sPLS/DIABLO（mixOmics）
# Run sequentially in a new working directory.
# Required packages: ggplot2, jsonlite, mixOmics, patchwork, ragg, readr, scales, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "scripts/run_article47_diablo.R",
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

set.seed(20260747)
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
metadata$StudyGroup <- factor(metadata$StudyGroup, levels = c("Control", "CD", "UC"))
stopifnot(
  nrow(otutab) == 250L, nrow(metabolites) == 277L,
  identical(colnames(otutab), rownames(metadata)),
  identical(colnames(metabolites), rownames(metadata))
)

set.seed(20260726)
train_ids <- unlist(lapply(levels(metadata$StudyGroup), function(group_name) {
  ids <- rownames(metadata)[metadata$StudyGroup == group_name]
  sample(ids, floor(0.70 * length(ids)), replace = FALSE)
}), use.names = FALSE)
train_ids <- rownames(metadata)[rownames(metadata) %in% train_ids]
test_ids <- setdiff(rownames(metadata), train_ids)
split_ledger <- rbind(
  Train = table(metadata[train_ids, "StudyGroup"]),
  Test = table(metadata[test_ids, "StudyGroup"])
)
stopifnot(nrow(metadata) == 220L, length(train_ids) == 153L, length(test_ids) == 67L)
split_ledger

run_analysis <- function(command, args) {
  output <- system2(command, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(paste(tail(output, 20), collapse = "\n"))
  invisible(output)
}
run_analysis(file.path(R.home("bin"), "Rscript"), c("--vanilla", "scripts/run_article47_diablo.R"))

min_train_samples <- ceiling(0.10 * length(train_ids))
prevalent_microbes <- rownames(otutab)[
  rowSums(otutab[, train_ids, drop = FALSE] > 0) >= min_train_samples
]
train_log <- log(t(otutab[prevalent_microbes, train_ids, drop = FALSE]) + 0.5)
train_clr_all <- train_log - rowMeans(train_log)
microbe_features <- names(sort(apply(train_clr_all, 2, var), decreasing = TRUE))[1:60]
train_metabolite_all <- log1p(t(metabolites[, train_ids, drop = FALSE]))
metabolite_features <- names(sort(
  apply(train_metabolite_all, 2, var), decreasing = TRUE
))[1:80]

clr_transform <- function(x) {
  logged <- log(t(x) + 0.5)
  logged - rowMeans(logged)
}
scale_from_train <- function(train, test) {
  center <- colMeans(train)
  spread <- apply(train, 2, sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  list(
    train = sweep(sweep(train, 2, center, "-"), 2, spread, "/"),
    test = sweep(sweep(test, 2, center, "-"), 2, spread, "/")
  )
}
microbe_scaled <- scale_from_train(
  clr_transform(otutab[microbe_features, train_ids, drop = FALSE]),
  clr_transform(otutab[microbe_features, test_ids, drop = FALSE])
)
metabolite_scaled <- scale_from_train(
  log1p(t(metabolites[metabolite_features, train_ids, drop = FALSE])),
  log1p(t(metabolites[metabolite_features, test_ids, drop = FALSE]))
)
train_blocks <- list(
  Microbiome = microbe_scaled$train,
  Metabolome = metabolite_scaled$train
)
test_blocks <- list(
  Microbiome = microbe_scaled$test,
  Metabolome = metabolite_scaled$test
)

result_dir <- "results/47-spls-diablo"
spls_correlations <- readr::read_tsv(
  file.path(result_dir, "spls-correlations.tsv"), show_col_types = FALSE
)
diablo_scores <- readr::read_tsv(
  file.path(result_dir, "diablo-train-scores.tsv"), show_col_types = FALSE
)
test_predictions <- readr::read_tsv(
  file.path(result_dir, "test-predictions.tsv"), show_col_types = FALSE
)
test_auc <- readr::read_tsv(
  file.path(result_dir, "test-auc.tsv"), show_col_types = FALSE
)
stability <- readr::read_tsv(
  file.path(result_dir, "diablo-selection-stability.tsv"), show_col_types = FALSE
)
diablo_summary <- jsonlite::read_json(
  file.path(result_dir, "summary.json"), simplifyVector = TRUE
)
stopifnot(
  nrow(test_predictions) == 67L,
  abs(diablo_summary$test_accuracy - 0.6567164) < 1e-4,
  abs(diablo_summary$test_macro_auc - 0.7984386) < 1e-4,
  diablo_summary$stable_selected_features_rate_ge_0_8 == 26L
)
diablo_summary[c(
  "repeated_cv_weighted_vote_ber_component2", "test_accuracy",
  "test_balanced_accuracy", "test_macro_auc"
)]

spls_long <- rbind(
  data.frame(Component = spls_correlations$Component, Split = "Training",
             Correlation = spls_correlations$TrainCorrelation),
  data.frame(Component = spls_correlations$Component, Split = "Locked test",
             Correlation = spls_correlations$TestCorrelation)
)
p47_1 <- ggplot(spls_long, aes(Component, Correlation, fill = Split)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = sprintf("%.2f", Correlation)),
            position = position_dodge(width = 0.72), vjust = -0.35,
            size = 3, family = font_pub) +
  scale_fill_manual(values = c(
    Training = unname(pal_pub["grey"]),
    `Locked test` = unname(pal_pub["blue"])
  )) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Cross-block correlation", fill = "Split",
       title = "sPLS generalization audit") + theme_pub()
save_pub(p47_1, "figures/47-spls-diablo/47-1-spls-generalization", 89, 70)
p47_1

group_colors <- c(
  Control = unname(pal_pub["blue"]),
  CD = unname(pal_pub["vermillion"]),
  UC = unname(pal_pub["orange"])
)
p47_2 <- ggplot(diablo_scores, aes(Component1, Component2, colour = StudyGroup)) +
  geom_hline(yintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
  geom_point(size = 1.5, alpha = 0.76) +
  stat_ellipse(linewidth = 0.55, level = 0.80, show.legend = FALSE) +
  facet_wrap(~ Block, scales = "free") +
  scale_color_manual(values = group_colors) +
  labs(x = "Component 1", y = "Component 2", colour = "Study group",
       title = "Supervised training representation",
       subtitle = "Interpretation plot; performance is evaluated on the locked test set") +
  theme_pub()
save_pub(p47_2, "figures/47-spls-diablo/47-2-diablo-scores", 180, 82)
p47_2

confusion <- as.data.frame(table(
  Truth = factor(test_predictions$Truth, levels = c("Control", "CD", "UC")),
  Predicted = factor(test_predictions$Predicted, levels = c("Control", "CD", "UC"))
))
confusion$WithinTruth <- ave(confusion$Freq, confusion$Truth, FUN = function(x) x / sum(x))
p_confusion <- ggplot(confusion, aes(Predicted, Truth, fill = WithinTruth)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%d\n%.0f%%", Freq, 100 * WithinTruth)),
            family = font_pub, size = 3.2) +
  scale_fill_gradient(low = "#F7FBFF", high = pal_pub["blue"], limits = c(0, 1)) +
  labs(x = "Predicted class", y = "True class", fill = "Row fraction",
       title = "Confusion matrix") + theme_pub(9) + theme(legend.position = "none")

test_auc$Class <- factor(test_auc$Class, levels = rev(c("Control", "CD", "UC")))
p_auc <- ggplot(test_auc, aes(AUC, Class)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.16,
                 colour = pal_pub["grey"], linewidth = 0.55) +
  geom_point(size = 2.8, colour = pal_pub["vermillion"]) +
  scale_x_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "One-vs-rest AUC (95% CI)", y = NULL,
       title = "Class discrimination") + theme_pub(9) + theme(legend.position = "none")
p47_3 <- p_confusion + p_auc + plot_annotation(
  title = "Locked-test performance",
  subtitle = "Accuracy 0.657; balanced accuracy 0.676; macro AUC 0.798"
)
save_pub(p47_3, "figures/47-spls-diablo/47-3-test-performance", 180, 82)
p47_3

stable_plot <- stability[order(-stability$FoldSelectionRate), ]
stable_plot <- stable_plot[!duplicated(paste(stable_plot$Block, stable_plot$FeatureID)), ]
stable_plot <- head(stable_plot, 20)
stable_plot$Label <- ifelse(
  nchar(stable_plot$DisplayName) > 28,
  paste0(substr(stable_plot$DisplayName, 1, 25), "..."),
  stable_plot$DisplayName
)
stable_plot$Label <- factor(stable_plot$Label, levels = rev(stable_plot$Label))
p47_4 <- ggplot(stable_plot, aes(FoldSelectionRate, Label, colour = Block)) +
  geom_vline(xintercept = 0.80, linetype = 2, colour = "#777777") +
  geom_segment(aes(x = 0, xend = FoldSelectionRate, yend = Label), linewidth = 0.6) +
  geom_point(size = 2.7) +
  scale_color_manual(values = c(
    Microbiome = unname(pal_pub["blue"]),
    Metabolome = unname(pal_pub["orange"])
  )) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Mean fold-selection rate", y = NULL, colour = "Block",
       title = "Feature stability", subtitle = "Dashed line: pre-specified 80% threshold") +
  theme_pub(9) + theme(legend.position = "right")
save_pub(p47_4, "figures/47-spls-diablo/47-4-feature-stability", 180, 105)
p47_4
