# 跨队列 meta 分析与外部验证
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, metafor, pROC, ragg, ranger, readr, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/cross-cohort-crc/crc_xiang/metadata.tsv",
  "data/small/cross-cohort-crc/crc_xiang/otutab.tsv",
  "data/small/cross-cohort-crc/crc_xiang/taxonomy.tsv",
  "data/small/cross-cohort-crc/crc_zackular/metadata.tsv",
  "data/small/cross-cohort-crc/crc_zackular/otutab.tsv",
  "data/small/cross-cohort-crc/crc_zackular/taxonomy.tsv",
  "data/small/cross-cohort-crc/crc_zhao/metadata.tsv",
  "data/small/cross-cohort-crc/crc_zhao/otutab.tsv",
  "data/small/cross-cohort-crc/crc_zhao/taxonomy.tsv",
  "data/small/cross-cohort-crc/source-summary.json"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(metafor)

set.seed(20260744)
font_pub <- "sans"
pal_pub <- c(
  blue = "#0072B2", orange = "#E69F00", green = "#009E73",
  vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
  yellow = "#F0E442", grey = "#6B7280"
)
scale_color_pub <- function(..., values = unname(pal_pub)) ggplot2::scale_color_manual(..., values = values)
scale_fill_pub <- function(..., values = unname(pal_pub)) ggplot2::scale_fill_manual(..., values = values)
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
clean_genus <- function(x) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  x[is.na(x) | x == "" | grepl("unknown|unclassified|uncultured|metagenome", x,
                                  ignore.case = TRUE)] <- NA_character_
  x
}
read_study <- function(study) {
  directory <- file.path("data/small/cross-cohort-crc", study)
  raw <- read_keyed_tsv(file.path(directory, "otutab.tsv"))
  taxonomy <- read_keyed_tsv(file.path(directory, "taxonomy.tsv"))
  metadata <- read_keyed_tsv(file.path(directory, "metadata.tsv"))
  counts <- as.matrix(data.frame(lapply(raw, as.integer), check.names = FALSE))
  rownames(counts) <- rownames(raw)
  counts <- counts[, rownames(metadata), drop = FALSE]
  genus <- clean_genus(taxonomy[rownames(counts), "Genus"])
  keep <- !is.na(genus) & rowSums(counts) > 0
  genus_counts <- rowsum(counts[keep, , drop = FALSE], group = genus[keep], reorder = TRUE)
  metadata$Group <- factor(metadata$Group, levels = c("Control", "CRC"))
  list(counts = genus_counts, metadata = metadata, taxonomy = taxonomy)
}
study_names <- c("crc_xiang", "crc_zhao", "crc_zackular")
studies <- setNames(lapply(study_names, read_study), study_names)
study_label <- c(crc_xiang = "Xiang", crc_zhao = "Zhao", crc_zackular = "Zackular")

dimension_ledger <- do.call(rbind, lapply(study_names, function(study) {
  data.frame(
    Study = study_label[[study]], Genera = nrow(studies[[study]]$counts),
    Samples = ncol(studies[[study]]$counts), ReadsAssignedToGenus = sum(studies[[study]]$counts),
    Controls = sum(studies[[study]]$metadata$Group == "Control"),
    CRC = sum(studies[[study]]$metadata$Group == "CRC")
  )
}))
stopifnot(
  identical(dimension_ledger$Samples, c(43L, 102L, 60L)),
  identical(dimension_ledger$Controls, c(22L, 56L, 30L)),
  identical(dimension_ledger$CRC, c(21L, 46L, 30L))
)
dimension_ledger

eligible_by_study <- lapply(studies, function(study) {
  rownames(study$counts)[rowMeans(study$counts > 0) >= 0.10 & rowSums(study$counts) >= 20]
})
eligibility_count <- table(unlist(eligible_by_study))
meta_genera <- names(eligibility_count)[eligibility_count >= 2L]
stopifnot(length(meta_genera) >= 20L)

align_counts <- function(counts, features) {
  aligned <- matrix(0L, nrow = length(features), ncol = ncol(counts),
                    dimnames = list(features, colnames(counts)))
  common <- intersect(features, rownames(counts))
  aligned[common, ] <- counts[common, , drop = FALSE]
  aligned
}
clr_transform <- function(counts, pseudocount = 0.5) {
  logged <- log(counts + pseudocount)
  sweep(logged, 2L, colMeans(logged), "-")
}
study_clr <- lapply(studies, function(study) clr_transform(align_counts(study$counts, meta_genera)))

study_effects <- do.call(rbind, lapply(study_names, function(study) {
  group <- studies[[study]]$metadata$Group
  do.call(rbind, lapply(meta_genera, function(genus) {
    fit <- stats::lm(study_clr[[study]][genus, ] ~ group)
    coefficient <- summary(fit)$coefficients["groupCRC", ]
    data.frame(
      Study = study_label[[study]], Genus = genus,
      Effect = unname(coefficient["Estimate"]), SE = unname(coefficient["Std. Error"]),
      PValue = unname(coefficient["Pr(>|t|)"]), stringsAsFactors = FALSE
    )
  }))
}))

meta_one_genus <- function(genus) {
  x <- study_effects[study_effects$Genus == genus & is.finite(study_effects$SE) & study_effects$SE > 0, ]
  fit <- metafor::rma.uni(yi = x$Effect, sei = x$SE, method = "REML")
  data.frame(
    Genus = genus, Studies = fit$k, PooledEffect = as.numeric(fit$b),
    SE = fit$se, CILower = fit$ci.lb, CIUpper = fit$ci.ub, PValue = fit$pval,
    Tau2 = fit$tau2, I2 = fit$I2, DirectionAgreement = abs(sum(sign(x$Effect))) / nrow(x),
    stringsAsFactors = FALSE
  )
}
meta_results <- do.call(rbind, lapply(meta_genera, meta_one_genus))
meta_results$AdjustedP <- stats::p.adjust(meta_results$PValue, method = "BH")
meta_results <- meta_results[order(meta_results$AdjustedP, -abs(meta_results$PooledEffect)), ]
stopifnot(nrow(meta_results) == length(meta_genera), all(meta_results$I2 >= 0 & meta_results$I2 <= 100))
head(meta_results, 12)

auc_score <- function(truth, probability) {
  as.numeric(pROC::auc(truth, probability, levels = c("Control", "CRC"),
                       direction = "<", quiet = TRUE))
}
select_training_features <- function(training_studies, max_features = 80L) {
  eligible <- lapply(training_studies, function(study) {
    rownames(study$counts)[rowMeans(study$counts > 0) >= 0.10 & rowSums(study$counts) >= 20]
  })
  shared <- Reduce(intersect, eligible)
  combined <- do.call(cbind, lapply(training_studies, function(study) align_counts(study$counts, shared)))
  transformed <- t(clr_transform(combined))
  variance <- apply(transformed, 2L, stats::var)
  names(head(sort(variance, decreasing = TRUE), max_features))
}
prepare_model_matrix <- function(study, features) {
  t(clr_transform(align_counts(study$counts, features)))
}
grid <- expand.grid(MtryFraction = c(0.10, 0.30, 0.60), MinNodeSize = c(3L, 10L))
external_predictions <- list()
tuning_results <- list()

for (held_out in study_names) {
  training_names <- setdiff(study_names, held_out)
  inner_score <- numeric(nrow(grid))
  for (grid_id in seq_len(nrow(grid))) {
    direction_auc <- numeric(length(training_names))
    for (inner_id in seq_along(training_names)) {
      inner_valid <- training_names[inner_id]
      inner_train <- setdiff(training_names, inner_valid)
      features <- select_training_features(studies[inner_train])
      x_train <- do.call(rbind, lapply(studies[inner_train], prepare_model_matrix, features = features))
      y_train <- factor(unlist(lapply(studies[inner_train], function(x) as.character(x$metadata$Group))),
                        levels = c("Control", "CRC"))
      x_valid <- prepare_model_matrix(studies[[inner_valid]], features)
      y_valid <- studies[[inner_valid]]$metadata$Group
      fit <- ranger::ranger(
        x = x_train, y = y_train, probability = TRUE, num.trees = 500L,
        mtry = max(1L, round(grid$MtryFraction[grid_id] * ncol(x_train))),
        min.node.size = grid$MinNodeSize[grid_id], importance = "none",
        seed = 20260744 + 1000 * match(held_out, study_names) + 100 * grid_id + inner_id,
        num.threads = 1L
      )
      direction_auc[inner_id] <- auc_score(y_valid, predict(fit, data = x_valid)$predictions[, "CRC"])
    }
    inner_score[grid_id] <- mean(direction_auc)
    tuning_results[[length(tuning_results) + 1L]] <- data.frame(
      HeldOut = held_out, MtryFraction = grid$MtryFraction[grid_id],
      MinNodeSize = grid$MinNodeSize[grid_id], MeanTrainingStudyAUC = inner_score[grid_id]
    )
  }
  best_id <- which.max(inner_score)
  features <- select_training_features(studies[training_names])
  x_train <- do.call(rbind, lapply(studies[training_names], prepare_model_matrix, features = features))
  y_train <- factor(unlist(lapply(studies[training_names], function(x) as.character(x$metadata$Group))),
                    levels = c("Control", "CRC"))
  x_test <- prepare_model_matrix(studies[[held_out]], features)
  y_test <- studies[[held_out]]$metadata$Group
  fit <- ranger::ranger(
    x = x_train, y = y_train, probability = TRUE, num.trees = 1500L,
    mtry = max(1L, round(grid$MtryFraction[best_id] * ncol(x_train))),
    min.node.size = grid$MinNodeSize[best_id], importance = "permutation",
    seed = 20260744 + 100000 * match(held_out, study_names), num.threads = 1L
  )
  external_predictions[[held_out]] <- data.frame(
    SampleID = rownames(x_test), HeldOutStudy = study_label[[held_out]],
    Truth = as.character(y_test), ProbabilityCRC = predict(fit, data = x_test)$predictions[, "CRC"],
    TrainingStudies = paste(study_label[training_names], collapse = " + "),
    Features = length(features), stringsAsFactors = FALSE
  )
}
external_predictions <- do.call(rbind, external_predictions)
tuning_results <- do.call(rbind, tuning_results)
external_performance <- do.call(rbind, lapply(unique(external_predictions$HeldOutStudy), function(study) {
  x <- external_predictions[external_predictions$HeldOutStudy == study, ]
  data.frame(HeldOutStudy = study, Samples = nrow(x), AUC = auc_score(
    factor(x$Truth, levels = c("Control", "CRC")), x$ProbabilityCRC
  ), Brier = mean((as.numeric(x$Truth == "CRC") - x$ProbabilityCRC)^2))
}))
external_performance$MacroMeanAUC <- mean(external_performance$AUC)
stopifnot(nrow(external_predictions) == 205L, nrow(external_performance) == 3L)
external_performance

ordination_genera <- meta_genera
combined_clr <- do.call(cbind, lapply(study_names, function(study) {
  clr_transform(align_counts(studies[[study]]$counts, ordination_genera))
}))
pca <- stats::prcomp(t(combined_clr), center = FALSE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
ordination <- data.frame(
  SampleID = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  Study = rep(study_label[study_names], times = vapply(studies, function(x) ncol(x$counts), integer(1))),
  Group = unlist(lapply(studies, function(x) as.character(x$metadata$Group))),
  stringsAsFactors = FALSE
)
p_pcoa <- ggplot(ordination, aes(PC1, PC2, colour = Study, shape = Group)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_pub(values = c(Xiang = pal_pub[["blue"]], Zhao = pal_pub[["orange"]], Zackular = pal_pub[["green"]])) +
  labs(title = "Study structure precedes modeling",
       x = sprintf("Aitchison PC1 (%.1f%%)", variance[1]),
       y = sprintf("Aitchison PC2 (%.1f%%)", variance[2]),
       caption = "The same outcome-blind genus universe is used across studies.") + theme_pub()
save_pub(p_pcoa, "figures/44-cohort-pcoa", 115, 95)
p_pcoa

forest_genus <- meta_results$Genus[1]
forest_data <- study_effects[study_effects$Genus == forest_genus, ]
forest_data$CILower <- forest_data$Effect - 1.96 * forest_data$SE
forest_data$CIUpper <- forest_data$Effect + 1.96 * forest_data$SE
pooled_row <- meta_results[meta_results$Genus == forest_genus, ]
forest_plot_data <- rbind(
  data.frame(Label = forest_data$Study, Effect = forest_data$Effect,
             CILower = forest_data$CILower, CIUpper = forest_data$CIUpper, Type = "Study"),
  data.frame(Label = "REML pooled", Effect = pooled_row$PooledEffect,
             CILower = pooled_row$CILower, CIUpper = pooled_row$CIUpper, Type = "Pooled")
)
forest_plot_data$Label <- factor(forest_plot_data$Label, levels = rev(forest_plot_data$Label))
p_forest <- ggplot(forest_plot_data, aes(Effect, Label, colour = Type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_pub(values = c(Study = pal_pub[["blue"]], Pooled = pal_pub[["vermillion"]])) +
  labs(title = "Random-effects meta-analysis",
       subtitle = sprintf("%s · I² = %.1f%% · BH p = %.3g", forest_genus, pooled_row$I2, pooled_row$AdjustedP),
       x = "CRC − control CLR mean difference", y = NULL,
       caption = "Study intervals are estimate ± 1.96 SE; pooled interval is REML.") +
  theme_pub() + theme(legend.position = "none")
save_pub(p_forest, "figures/44-meta-forest", 115, 82)
p_forest

label_meta <- unique(rbind(
  head(meta_results[order(meta_results$AdjustedP, -abs(meta_results$PooledEffect)), ], 3),
  head(meta_results[order(-meta_results$I2, na.last = TRUE), ], 3)
))
p_heterogeneity <- ggplot(meta_results, aes(PooledEffect, I2)) +
  geom_hline(yintercept = 50, linetype = 3, colour = "#B3B3B3") +
  geom_vline(xintercept = 0, linetype = 3, colour = "#B3B3B3") +
  geom_point(aes(colour = -log10(pmax(AdjustedP, 1e-12)), size = abs(PooledEffect)), alpha = 0.75) +
  ggrepel::geom_text_repel(data = label_meta, aes(label = Genus), size = 2.6,
                           max.overlaps = Inf, min.segment.length = 0,
                           box.padding = 0.4) +
  scale_colour_gradient(low = pal_pub[["sky"]], high = pal_pub[["vermillion"]]) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = "Effect–heterogeneity audit", x = "Random-effects pooled CLR difference",
       y = "I² (%)", colour = "−log10 BH p", size = "|Effect|",
       caption = "With three studies, I² is descriptive.") +
  theme_pub()
save_pub(p_heterogeneity, "figures/44-heterogeneity", 135, 100)
p_heterogeneity

external_roc <- do.call(rbind, lapply(unique(external_predictions$HeldOutStudy), function(study) {
  x <- external_predictions[external_predictions$HeldOutStudy == study, ]
  roc <- pROC::roc(factor(x$Truth, levels = c("Control", "CRC")), x$ProbabilityCRC,
                   levels = c("Control", "CRC"), direction = "<", quiet = TRUE)
  data.frame(HeldOutStudy = study, FalsePositiveRate = 1 - roc$specificities,
             TruePositiveRate = roc$sensitivities)
}))
auc_labels <- merge(external_performance, data.frame(
  HeldOutStudy = unique(external_predictions$HeldOutStudy), x = 0.58, y = c(0.18, 0.30, 0.42)
), by = "HeldOutStudy")
p_external <- ggplot(external_roc, aes(FalsePositiveRate, TruePositiveRate, colour = HeldOutStudy)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_path(linewidth = 1) +
  geom_text(data = auc_labels, aes(x, y, label = sprintf("%s AUC %.2f", HeldOutStudy, AUC)),
            inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_color_pub(values = c(Xiang = pal_pub[["blue"]], Zhao = pal_pub[["orange"]], Zackular = pal_pub[["green"]])) +
  coord_equal() +
  labs(title = "Leave-one-study-out external validation",
       subtitle = sprintf("Macro-mean AUC = %.2f", mean(external_performance$AUC)),
       x = "False-positive rate", y = "True-positive rate",
       caption = "Each held-out study was unseen during model development.") +
  theme_pub() + theme(legend.position = "none")
save_pub(p_external, "figures/44-external-validation", 110, 100)
p_external

dir.create("results/44-cross-cohort-validation", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(dimension_ledger, "results/44-cross-cohort-validation/cohort-dimensions.tsv")
readr::write_tsv(study_effects, "results/44-cross-cohort-validation/study-genus-effects.tsv")
readr::write_tsv(meta_results, "results/44-cross-cohort-validation/random-effects-meta.tsv")
readr::write_tsv(external_predictions, "results/44-cross-cohort-validation/external-predictions.tsv")
readr::write_tsv(external_performance, "results/44-cross-cohort-validation/external-performance.tsv")
readr::write_tsv(tuning_results, "results/44-cross-cohort-validation/training-study-tuning.tsv")
