# 微生物组与时间结局：Cox、Kaplan–Meier 与时间依赖 ROC
# Run sequentially in a new working directory.
# Required packages: ggplot2, glmnet, ragg, readr, scales, survival, svglite, timeROC.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/survival-t1d/metadata.tsv",
  "data/small/survival-t1d/otutab.tsv",
  "data/small/survival-t1d/source-summary.json",
  "data/small/survival-t1d/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(ggplot2)
library(survival)

set.seed(20260745)
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
data_dir <- "data/small/survival-t1d"
otutab_raw <- read_keyed_tsv(file.path(data_dir, "otutab.tsv"))
taxonomy <- read_keyed_tsv(file.path(data_dir, "taxonomy.tsv"))
metadata <- read_keyed_tsv(file.path(data_dir, "metadata.tsv"))
otutab <- as.matrix(data.frame(lapply(otutab_raw, as.integer), check.names = FALSE))
rownames(otutab) <- rownames(otutab_raw)
otutab <- otutab[, rownames(metadata), drop = FALSE]
metadata$TimeWeeks <- as.numeric(metadata$TimeWeeks)
metadata$SamplingWeek <- as.numeric(metadata$SamplingWeek)
metadata$FollowUpWeeks <- as.numeric(metadata$FollowUpWeeks)
metadata$Event <- as.integer(metadata$Event)
metadata$Sex <- relevel(factor(metadata$Sex), ref = "F")
metadata$Treatment <- relevel(factor(metadata$Treatment), ref = "Control")

stopifnot(
  nrow(otutab) == 348L, ncol(otutab) == 173L, sum(otutab) == 3073108L,
  sum(metadata$Event) == 118L, sum(metadata$Event == 0L) == 55L,
  all(metadata$FollowUpWeeks > 0),
  identical(sort(unique(metadata$SamplingWeek)), c(6, 7, 8)),
  all(metadata$FollowUpWeeks == metadata$TimeWeeks - metadata$SamplingWeek),
  identical(colnames(otutab), rownames(metadata))
)
head(otutab[, 1:5])
head(taxonomy)
head(metadata)

clean_genus <- function(x) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  x[is.na(x) | x == "" | grepl("unknown|unclassified|uncultured|metagenome", x,
                                  ignore.case = TRUE)] <- NA_character_
  x
}
genus <- clean_genus(taxonomy[rownames(otutab), "Genus"])
keep_named <- !is.na(genus)
genus_counts <- rowsum(otutab[keep_named, , drop = FALSE], group = genus[keep_named], reorder = TRUE)
eligible <- rowMeans(genus_counts > 0) >= 0.10 & rowSums(genus_counts) >= 20
analysis_counts <- genus_counts[eligible, , drop = FALSE]
clr <- log(analysis_counts + 0.5)
clr <- sweep(clr, 2L, colMeans(clr), "-")

stopifnot(nrow(analysis_counts) >= 20L, all(colnames(analysis_counts) == rownames(metadata)))
data.frame(NamedGenera = nrow(genus_counts), EligibleGenera = nrow(analysis_counts),
           ReadsAssignedToNamedGenus = sum(genus_counts))

survival_object <- survival::Surv(metadata$FollowUpWeeks, metadata$Event)
km_treatment <- survival::survfit(survival_object ~ Treatment, data = metadata, conf.type = "log-log")
logrank <- survival::survdiff(survival_object ~ Treatment, data = metadata)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)

km_summary <- summary(km_treatment, censored = TRUE)
km_data <- data.frame(
  Time = km_summary$time, Survival = km_summary$surv,
  Lower = km_summary$lower, Upper = km_summary$upper,
  NRisk = km_summary$n.risk, NEvent = km_summary$n.event,
  Treatment = sub("Treatment=", "", km_summary$strata, fixed = TRUE),
  stringsAsFactors = FALSE
)
km_data <- rbind(
  data.frame(Time = 0, Survival = 1, Lower = 1, Upper = 1,
             NRisk = as.numeric(table(metadata$Treatment)), NEvent = 0,
             Treatment = names(table(metadata$Treatment))),
  km_data
)
data.frame(LogRankChiSquare = logrank$chisq, LogRankP = logrank_p)

fit_one_genus <- function(genus_name) {
  genus_z <- as.numeric(scale(clr[genus_name, ]))
  model_data <- data.frame(
    FollowUpWeeks = metadata$FollowUpWeeks, Event = metadata$Event,
    GenusZ = genus_z, Sex = metadata$Sex, Treatment = metadata$Treatment
  )
  fit <- survival::coxph(
    survival::Surv(FollowUpWeeks, Event) ~ GenusZ + Sex + Treatment,
    data = model_data, ties = "efron", x = TRUE
  )
  coefficient <- summary(fit)$coefficients["GenusZ", ]
  interval <- summary(fit)$conf.int["GenusZ", ]
  ph <- survival::cox.zph(fit, transform = "km")
  data.frame(
    Genus = genus_name, HR = interval["exp(coef)"],
    CILower = interval["lower .95"], CIUpper = interval["upper .95"],
    Coefficient = coefficient["coef"], SE = coefficient["se(coef)"],
    PValue = coefficient["Pr(>|z|)"], PHPValue = ph$table["GenusZ", "p"],
    stringsAsFactors = FALSE
  )
}
cox_results <- do.call(rbind, lapply(rownames(clr), fit_one_genus))
cox_results$AdjustedP <- stats::p.adjust(cox_results$PValue, method = "BH")
cox_results$PHAdjustedP <- stats::p.adjust(cox_results$PHPValue, method = "BH")
cox_results <- cox_results[order(cox_results$AdjustedP, -abs(cox_results$Coefficient)), ]
stopifnot(all(cox_results$HR > 0), all(cox_results$CILower > 0))
head(cox_results, 12)

make_stratified_folds <- function(strata, k, seed) {
  set.seed(seed)
  folds <- integer(length(strata))
  for (level in unique(as.character(strata))) {
    index <- which(as.character(strata) == level)
    index <- sample(index, length(index))
    folds[index] <- rep(seq_len(k), length.out = length(index))
  }
  folds
}
set.seed(20260745)
train_index <- unlist(lapply(split(seq_len(nrow(metadata)), metadata$Event), function(index) {
  sample(index, floor(0.70 * length(index)), replace = FALSE)
}))
train_index <- sort(train_index)
test_index <- setdiff(seq_len(nrow(metadata)), train_index)
train_ids <- rownames(metadata)[train_index]
test_ids <- rownames(metadata)[test_index]

training_keep <- rowMeans(genus_counts[, train_ids, drop = FALSE] > 0) >= 0.10 &
  rowSums(genus_counts[, train_ids, drop = FALSE]) >= 20
training_features <- rownames(genus_counts)[training_keep]
training_logged <- t(log(genus_counts[training_features, train_ids, drop = FALSE] + 0.5))
training_clr <- training_logged - rowMeans(training_logged)
training_variance <- apply(training_clr, 2L, stats::var)
training_features <- names(head(sort(training_variance, decreasing = TRUE), 50L))

make_matrix <- function(sample_ids, features) {
  x <- t(log(genus_counts[features, sample_ids, drop = FALSE] + 0.5))
  x - rowMeans(x)
}
x_train <- make_matrix(train_ids, training_features)
x_test <- make_matrix(test_ids, training_features)
y_train <- survival::Surv(metadata$FollowUpWeeks[train_index], metadata$Event[train_index])
fold_id <- make_stratified_folds(metadata$Event[train_index], 5L, 20260745)

set.seed(20260745)
ridge_cv <- glmnet::cv.glmnet(
  x_train, y_train, family = "cox", alpha = 0,
  foldid = fold_id, type.measure = "deviance", standardize = TRUE
)
lambda_locked <- ridge_cv$lambda.1se
train_risk <- as.numeric(predict(ridge_cv, newx = x_train, s = lambda_locked, type = "link"))
test_risk <- as.numeric(predict(ridge_cv, newx = x_test, s = lambda_locked, type = "link"))
training_cutoff <- median(train_risk)
test_risk_group <- factor(ifelse(test_risk > training_cutoff, "High", "Low"), levels = c("Low", "High"))

test_concordance <- survival::concordance(
  survival::Surv(metadata$FollowUpWeeks[test_index], metadata$Event[test_index]) ~ test_risk,
  reverse = TRUE
)$concordance
horizons <- c(10, 15, 20)
test_time_roc <- timeROC::timeROC(
  T = metadata$FollowUpWeeks[test_index], delta = metadata$Event[test_index],
  marker = test_risk, cause = 1, weighting = "marginal", times = horizons, iid = TRUE
)

split_ledger <- data.frame(
  Set = c("Training", "Test"),
  Samples = c(length(train_index), length(test_index)),
  Events = c(sum(metadata$Event[train_index]), sum(metadata$Event[test_index])),
  Features = c(length(training_features), length(training_features))
)
stopifnot(length(intersect(train_ids, test_ids)) == 0L, all(is.finite(test_time_roc$AUC)))
data.frame(Lambda1SE = lambda_locked, TestCIndex = test_concordance,
           HorizonWeeks = horizons, TestAUC = test_time_roc$AUC)

p_km <- ggplot(km_data, aes(Time, Survival, colour = Treatment, fill = Treatment)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.13, colour = NA) +
  geom_step(linewidth = 0.95) +
  scale_color_pub(values = c(Control = pal_pub[["blue"]], Tylosin = pal_pub[["orange"]])) +
  scale_fill_pub(values = c(Control = pal_pub[["blue"]], Tylosin = pal_pub[["orange"]])) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Time to T1D onset", subtitle = sprintf("Log-rank p = %.3g", logrank_p),
       x = "Weeks since baseline microbiome sampling", y = "T1D-free probability",
       caption = "Pointwise 95% CI; baseline samples collected at weeks 6–8.") +
  theme_pub()
save_pub(p_km, "figures/45-treatment-km", 110, 90)
p_km

forest_genus <- head(cox_results, 12)
forest_genus$Genus <- factor(forest_genus$Genus, levels = rev(forest_genus$Genus))
p_cox <- ggplot(forest_genus, aes(HR, Genus, colour = AdjustedP < 0.10)) +
  geom_vline(xintercept = 1, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_errorbarh(aes(xmin = CILower, xmax = CIUpper), height = 0.18, linewidth = 0.65) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  scale_color_pub(values = c(`TRUE` = pal_pub[["vermillion"]], `FALSE` = pal_pub[["blue"]])) +
  labs(title = "Genus-level adjusted Cox models", x = "Hazard ratio per 1-SD CLR", y = NULL,
       colour = "BH p < 0.10",
       caption = "Adjusted for Sex and Treatment; 95% Wald CIs.") +
  theme_pub()
save_pub(p_cox, "figures/45-taxon-cox", 120, 100)
p_cox

diagnostic_genus <- as.character(cox_results$Genus[1])
diagnostic_data <- data.frame(
  FollowUpWeeks = metadata$FollowUpWeeks, Event = metadata$Event,
  GenusZ = as.numeric(scale(clr[diagnostic_genus, ])),
  Sex = metadata$Sex, Treatment = metadata$Treatment
)
diagnostic_fit <- survival::coxph(
  survival::Surv(FollowUpWeeks, Event) ~ GenusZ + Sex + Treatment,
  data = diagnostic_data, ties = "efron", x = TRUE
)
diagnostic_zph <- survival::cox.zph(diagnostic_fit, transform = "km")
ph_plot_data <- data.frame(
  TransformedTime = diagnostic_zph$x,
  ScaledSchoenfeldResidual = diagnostic_zph$y[, "GenusZ"]
)
p_ph <- ggplot(ph_plot_data, aes(TransformedTime, ScaledSchoenfeldResidual)) +
  geom_hline(yintercept = coef(diagnostic_fit)["GenusZ"], linetype = 2, colour = pal_pub[["grey"]]) +
  geom_point(colour = pal_pub[["blue"]], alpha = 0.55, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, colour = pal_pub[["vermillion"]],
              fill = pal_pub[["vermillion"]], alpha = 0.13) +
  labs(title = paste("Proportional-hazards audit:", diagnostic_genus),
       subtitle = sprintf("Schoenfeld test p = %.3g", diagnostic_zph$table["GenusZ", "p"]),
       x = "KM-transformed event time", y = "Scaled Schoenfeld residual",
       caption = "A systematic time trend suggests a time-varying coefficient.") + theme_pub()
save_pub(p_ph, "figures/45-ph-diagnostics", 110, 90)
p_ph

time_roc_data <- do.call(rbind, lapply(seq_along(horizons), function(index) {
  data.frame(FalsePositiveRate = test_time_roc$FP[, index],
             TruePositiveRate = test_time_roc$TP[, index],
             Horizon = paste0(horizons[index], " weeks"))
}))
auc_text <- data.frame(
  Horizon = paste0(horizons, " weeks"),
  Label = sprintf("%d weeks: AUC %.2f", horizons, test_time_roc$AUC),
  x = 0.55, y = c(0.18, 0.30, 0.42)
)
p_time_roc <- ggplot(time_roc_data, aes(FalsePositiveRate, TruePositiveRate, colour = Horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = pal_pub[["grey"]]) +
  geom_path(linewidth = 1) +
  geom_text(data = auc_text, aes(x, y, label = Label, colour = Horizon),
            inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_color_pub(values = c(`10 weeks` = pal_pub[["blue"]],
                             `15 weeks` = pal_pub[["orange"]],
                             `20 weeks` = pal_pub[["green"]])) +
  coord_equal() +
  labs(title = "Locked-test time-dependent ROC",
       subtitle = sprintf("Test C-index = %.2f; n = %d", test_concordance, length(test_index)),
       x = "False-positive rate", y = "True-positive rate",
       caption = "Training-only preprocessing and tuning; test evaluated once.") +
  theme_pub() + theme(legend.position = "none")
save_pub(p_time_roc, "figures/45-time-dependent-roc", 110, 100)
p_time_roc

dir.create("results/45-survival-analysis", recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(km_data, "results/45-survival-analysis/treatment-km.tsv")
readr::write_tsv(cox_results, "results/45-survival-analysis/genus-adjusted-cox.tsv")
readr::write_tsv(split_ledger, "results/45-survival-analysis/prediction-split.tsv")
readr::write_tsv(data.frame(
  SampleID = test_ids, FollowUpWeeks = metadata$FollowUpWeeks[test_index],
  Event = metadata$Event[test_index], RiskScore = test_risk,
  LockedRiskGroup = test_risk_group
), "results/45-survival-analysis/locked-test-predictions.tsv")
readr::write_tsv(data.frame(
  Lambda1SE = lambda_locked, TrainingMedianCutoff = training_cutoff,
  TestCIndex = test_concordance, HorizonWeeks = horizons, TestAUC = test_time_roc$AUC
), "results/45-survival-analysis/prediction-performance.tsv")
