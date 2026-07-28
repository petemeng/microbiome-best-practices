#!/usr/bin/env Rscript

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(emmeans)
  library(ggplot2)
  library(jsonlite)
  library(nlme)
  library(patchwork)
  library(svglite)
  library(vegan)
})
source("R/theme_pub.R")

set.seed(20260751)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "longitudinal-dietswap")
out_dir <- file.path(root, "results", "51-longitudinal-analysis")
figure_dir <- file.path(root, "figures", "51-longitudinal-analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

otutab <- as.matrix(read.delim(
  file.path(data_dir, "otutab.tsv"),
  row.names = 1,
  check.names = FALSE
))
taxonomy <- read.delim(
  file.path(data_dir, "taxonomy.tsv"),
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  file.path(data_dir, "metadata.tsv"),
  row.names = 1,
  check.names = FALSE
)
storage.mode(otutab) <- "numeric"
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(colSums(otutab) > 0),
  all(otutab >= 0)
)

relative <- sweep(otutab, 2L, colSums(otutab), "/")
shannon <- -colSums(ifelse(relative > 0, relative * log(relative), 0))
metadata$Shannon <- shannon[rownames(metadata)]
metadata$Cohort <- factor(
  metadata$Cohort,
  levels = c("African American", "Rural African")
)
metadata$TimeFactor <- factor(metadata$Timepoint, levels = 1:6)
metadata$TimeLabel <- factor(
  paste0("T", metadata$Timepoint),
  levels = paste0("T", 1:6)
)

fit_ri_ml <- nlme::lme(
  Shannon ~ Cohort * TimeFactor,
  random = ~ 1 | SubjectID,
  data = metadata,
  method = "ML",
  na.action = na.omit,
  control = nlme::lmeControl(opt = "optim")
)
fit_ar1_ml <- update(
  fit_ri_ml,
  correlation = nlme::corAR1(form = ~ Timepoint | SubjectID)
)
fit_ri <- update(fit_ri_ml, method = "REML")
fit_ar1 <- update(fit_ar1_ml, method = "REML")
comparison <- as.data.frame(anova(fit_ri_ml, fit_ar1_ml))
comparison$Model <- c("Random intercept", "Random intercept + AR(1)")
comparison$AR1Phi <- c(NA_real_, unname(coef(
  fit_ar1$modelStruct$corStruct,
  unconstrained = FALSE
)))
comparison$LikelihoodRatioP <- c(NA_real_, comparison$`p-value`[2])
comparison <- comparison[, c(
  "Model", "df", "AIC", "BIC", "logLik", "AR1Phi",
  "LikelihoodRatioP"
)]

extract_fixed <- function(model, model_name) {
  estimates <- summary(model)$tTable
  critical <- stats::qt(0.975, df = estimates[, "DF"])
  data.frame(
    Model = model_name,
    Term = rownames(estimates),
    Estimate = estimates[, "Value"],
    SE = estimates[, "Std.Error"],
    DF = estimates[, "DF"],
    Lower = estimates[, "Value"] - critical * estimates[, "Std.Error"],
    Upper = estimates[, "Value"] + critical * estimates[, "Std.Error"],
    PValue = estimates[, "p-value"],
    row.names = NULL,
    check.names = FALSE
  )
}
fixed_effects <- rbind(
  extract_fixed(fit_ri, "Random intercept"),
  extract_fixed(fit_ar1, "Random intercept + AR(1)")
)

model_means <- as.data.frame(emmeans::emmeans(
  fit_ri,
  ~ Cohort * TimeFactor
))
model_means$Timepoint <- as.integer(as.character(model_means$TimeFactor))

bray <- vegan::vegdist(t(relative), method = "bray")
pcoa <- stats::cmdscale(bray, k = 5, eig = TRUE, add = TRUE)
pcoa_scores <- data.frame(
  SampleID = rownames(pcoa$points),
  Axis1 = pcoa$points[, 1],
  Axis2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)
pcoa_scores <- cbind(
  pcoa_scores,
  metadata[pcoa_scores$SampleID, c(
    "SubjectID", "Cohort", "Timepoint", "TimeLabel"
  ), drop = FALSE]
)
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
pcoa_variance <- data.frame(
  Axis = paste0("Axis", 1:5),
  Percent = 100 * pcoa$eig[1:5] / sum(positive_eigenvalues)
)
centroids <- aggregate(
  cbind(Axis1, Axis2) ~ Cohort + Timepoint,
  data = pcoa_scores,
  FUN = mean
)

bray_matrix <- as.matrix(bray)
volatility_rows <- lapply(split(metadata, metadata$SubjectID), function(subject) {
  subject <- subject[order(subject$Timepoint), , drop = FALSE]
  if (nrow(subject) < 2L) {
    return(NULL)
  }
  pairs <- lapply(seq_len(nrow(subject) - 1L), function(index) {
    if (subject$Timepoint[index + 1L] - subject$Timepoint[index] != 1L) {
      return(NULL)
    }
    from_id <- rownames(subject)[index]
    to_id <- rownames(subject)[index + 1L]
    data.frame(
      SubjectID = subject$SubjectID[index],
      Cohort = subject$Cohort[index],
      From = subject$Timepoint[index],
      To = subject$Timepoint[index + 1L],
      Window = paste0("T", subject$Timepoint[index], "-T", subject$Timepoint[index + 1L]),
      Volatility = unname(bray_matrix[from_id, to_id]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pairs)
})
volatility <- do.call(rbind, volatility_rows)
volatility$Cohort <- factor(
  volatility$Cohort,
  levels = levels(metadata$Cohort)
)
volatility$Window <- factor(
  volatility$Window,
  levels = paste0("T", 1:5, "-T", 2:6)
)
fit_volatility <- nlme::lme(
  Volatility ~ Cohort * Window,
  random = ~ 1 | SubjectID,
  data = volatility,
  method = "REML",
  control = nlme::lmeControl(opt = "optim")
)
volatility_effects <- extract_fixed(fit_volatility, "Volatility LMM")

paired_change <- do.call(rbind, lapply(
  split(metadata, metadata$SubjectID),
  function(subject) {
    if (!all(c(1L, 4L) %in% subject$Timepoint)) {
      return(NULL)
    }
    data.frame(
      SubjectID = subject$SubjectID[1],
      Cohort = subject$Cohort[1],
      DeltaT4T1 = subject$Shannon[subject$Timepoint == 4L][1] -
        subject$Shannon[subject$Timepoint == 1L][1]
    )
  }
))

set.seed(20260751)
bootstrap_draws <- do.call(rbind, lapply(
  levels(metadata$Cohort),
  function(cohort_name) {
    values <- paired_change$DeltaT4T1[paired_change$Cohort == cohort_name]
    data.frame(
      Cohort = cohort_name,
      Draw = seq_len(5000),
      MeanChange = replicate(
        5000,
        mean(sample(values, length(values), replace = TRUE))
      )
    )
  }
))
bootstrap_summary <- do.call(rbind, lapply(
  split(bootstrap_draws$MeanChange, bootstrap_draws$Cohort),
  function(values) data.frame(
    Mean = mean(values),
    Lower = unname(stats::quantile(values, 0.025)),
    Upper = unname(stats::quantile(values, 0.975))
  )
))
bootstrap_summary$Cohort <- rownames(bootstrap_summary)
rownames(bootstrap_summary) <- NULL

write.table(
  cbind(SampleID = rownames(metadata), metadata),
  file.path(out_dir, "sample-metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (object in c(
  "comparison", "fixed_effects", "model_means", "pcoa_scores",
  "pcoa_variance", "volatility", "volatility_effects",
  "paired_change", "bootstrap_summary"
)) {
  write.table(
    get(object),
    file.path(out_dir, paste0(gsub("_", "-", object), ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

cohort_colours <- c(
  "African American" = pal_pub[["blue"]],
  "Rural African" = pal_pub[["orange"]]
)

figure_1 <- ggplot(metadata, aes(Timepoint, Shannon, colour = Cohort)) +
  geom_line(aes(group = SubjectID), linewidth = 0.35, alpha = 0.17) +
  geom_point(size = 0.9, alpha = 0.25) +
  geom_ribbon(
    data = model_means,
    aes(x = Timepoint, ymin = lower.CL, ymax = upper.CL, fill = Cohort),
    inherit.aes = FALSE,
    alpha = 0.16,
    colour = NA
  ) +
  geom_line(
    data = model_means,
    aes(Timepoint, emmean, colour = Cohort, group = Cohort),
    linewidth = 1.05
  ) +
  geom_point(
    data = model_means,
    aes(Timepoint, emmean, colour = Cohort),
    size = 2.1
  ) +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Within-subject diversity trajectories",
    subtitle = "Thin lines: participants; thick lines and bands: mixed-model means and 95% CI",
    x = "Study timepoint",
    y = "Shannon diversity",
    colour = "Cohort",
    fill = "Cohort"
  ) +
  theme_pub(10)

figure_2 <- ggplot(
  pcoa_scores,
  aes(Axis1, Axis2, colour = Cohort)
) +
  geom_path(
    aes(group = SubjectID),
    linewidth = 0.3,
    alpha = 0.13,
    arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  geom_point(alpha = 0.18, size = 0.8) +
  geom_path(
    data = centroids,
    aes(group = Cohort),
    linewidth = 1.2,
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
  ) +
  geom_point(data = centroids, aes(fill = Cohort), shape = 21, size = 3, colour = "white") +
  geom_text(
    data = centroids[centroids$Timepoint %in% c(1L, 4L, 6L), ],
    aes(label = paste0("T", Timepoint)),
    colour = "#202020",
    nudge_y = 0.008,
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Bray-Curtis trajectories retain subject history",
    subtitle = "Faint arrows show participants; labelled paths show cohort centroids",
    x = sprintf("PCoA 1 (%.1f%%)", pcoa_variance$Percent[1]),
    y = sprintf("PCoA 2 (%.1f%%)", pcoa_variance$Percent[2]),
    colour = "Cohort",
    fill = "Cohort"
  ) +
  coord_equal() +
  theme_pub(10)

figure_3 <- ggplot(volatility, aes(Window, Volatility, fill = Cohort)) +
  geom_boxplot(
    position = position_dodge(width = 0.72),
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.72
  ) +
  geom_point(
    aes(colour = Cohort),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72),
    size = 0.8,
    alpha = 0.45
  ) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Volatility is a within-person distance",
    subtitle = "Bray-Curtis distance between consecutive observed timepoints",
    x = "Consecutive time window",
    y = "Bray-Curtis volatility",
    colour = "Cohort",
    fill = "Cohort"
  ) +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

interaction_effects <- fixed_effects[
  grepl("CohortRural African:TimeFactor", fixed_effects$Term),
]
interaction_effects$Timepoint <- sub(
  "CohortRural African:TimeFactor",
  "T",
  interaction_effects$Term,
  fixed = TRUE
)
audit_a <- ggplot(
  interaction_effects,
  aes(Estimate, Timepoint, colour = Model)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper),
    position = position_dodge(width = 0.45),
    height = 0
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8) +
  scale_colour_manual(values = c(
    "Random intercept" = pal_pub[["blue"]],
    "Random intercept + AR(1)" = pal_pub[["vermillion"]]
  )) +
  labs(
    title = "A. Correlation-structure audit",
    x = "Cohort-by-time interaction estimate (95% CI)",
    y = "Timepoint vs T1",
    colour = "Model"
  ) +
  theme_pub(9)

audit_b <- ggplot(bootstrap_draws, aes(MeanChange, fill = Cohort)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#555555") +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "B. Participant-level bootstrap",
    subtitle = "Mean Shannon change, T4 minus T1",
    x = "Bootstrapped mean change",
    y = "Density",
    fill = "Cohort"
  ) +
  theme_pub(9)

figure_4 <- audit_a + audit_b + plot_layout(widths = c(1.25, 1))

save_pub(figure_1, file.path(figure_dir, "51-1-alpha-trajectories"), 180, 125)
save_pub(figure_2, file.path(figure_dir, "51-2-pcoa-trajectories"), 180, 135)
save_pub(figure_3, file.path(figure_dir, "51-3-volatility"), 180, 125)
save_pub(figure_4, file.path(figure_dir, "51-4-model-audit"), 180, 115)

summary <- list(
  samples = nrow(metadata),
  subjects = length(unique(metadata$SubjectID)),
  complete_consecutive_pairs = nrow(volatility),
  random_intercept_aic = unname(AIC(fit_ri_ml)),
  ar1_aic = unname(AIC(fit_ar1_ml)),
  ar1_phi = unname(comparison$AR1Phi[2]),
  ar1_likelihood_ratio_p = unname(comparison$LikelihoodRatioP[2]),
  cohort_time_t4_interaction = unname(
    summary(fit_ri)$tTable["CohortRural African:TimeFactor4", "Value"]
  ),
  cohort_time_t4_p = unname(
    summary(fit_ri)$tTable["CohortRural African:TimeFactor4", "p-value"]
  ),
  note = "Association estimates describe the observed six-timepoint cohort and do not by themselves establish a diet-mediated causal effect."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "summary.json")
)

message("Article 51 complete: ", nrow(metadata), " samples, ", nrow(volatility), " consecutive pairs.")
