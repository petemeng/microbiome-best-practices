#!/usr/bin/env Rscript

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(mediation)
  library(patchwork)
  library(svglite)
})
source("R/theme_pub.R")

set.seed(20260753)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "paired-ibd-multiomics")
out_dir <- file.path(root, "results", "53-mediation-analysis")
figure_dir <- file.path(root, "figures", "53-mediation-analysis")
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
  taxonomy["MB_G005", "Genus"] == "Faecalibacterium"
)

complete <- metadata$StudyGroup %in% c("CD", "UC") &
  !is.na(metadata$FecalCalprotectin) &
  metadata$steroids %in% c("No", "Yes") &
  !is.na(metadata$Age) &
  metadata$antibiotic %in% c("No", "Yes") &
  metadata$immunosuppressant %in% c("No", "Yes") &
  metadata$mesalamine %in% c("No", "Yes")
sample_ids <- rownames(metadata)[complete]
stopifnot(length(sample_ids) == 108L)

z_score <- function(value) as.numeric(scale(value))
make_mediator <- function(pseudocount) {
  logged <- log(otutab[, sample_ids, drop = FALSE] + pseudocount)
  clr <- sweep(logged, 2L, colMeans(logged), "-")
  z_score(clr["MB_G005", ])
}

analysis_data <- data.frame(
  SampleID = sample_ids,
  SteroidExposure = factor(
    metadata[sample_ids, "steroids"],
    levels = c("No", "Yes")
  ),
  StudyGroup = factor(
    metadata[sample_ids, "StudyGroup"],
    levels = c("CD", "UC")
  ),
  AgeZ = z_score(metadata[sample_ids, "Age"]),
  Antibiotic = factor(metadata[sample_ids, "antibiotic"], levels = c("No", "Yes")),
  Immunosuppressant = factor(
    metadata[sample_ids, "immunosuppressant"],
    levels = c("No", "Yes")
  ),
  Mesalamine = factor(metadata[sample_ids, "mesalamine"], levels = c("No", "Yes")),
  FaecalibacteriumCLR = make_mediator(0.5),
  LogCalprotectin = log1p(metadata[sample_ids, "FecalCalprotectin"]),
  stringsAsFactors = FALSE
)

fit_models <- function(data, simulations, seed) {
  mediator_model <- stats::lm(
    FaecalibacteriumCLR ~ SteroidExposure + StudyGroup + AgeZ +
      Antibiotic + Immunosuppressant + Mesalamine,
    data = data
  )
  outcome_model <- stats::lm(
    LogCalprotectin ~ SteroidExposure + FaecalibacteriumCLR +
      StudyGroup + AgeZ + Antibiotic + Immunosuppressant + Mesalamine,
    data = data
  )
  set.seed(seed)
  mediation_fit <- mediation::mediate(
    mediator_model,
    outcome_model,
    treat = "SteroidExposure",
    mediator = "FaecalibacteriumCLR",
    control.value = "No",
    treat.value = "Yes",
    sims = simulations,
    boot = FALSE,
    robustSE = TRUE
  )
  list(
    mediator = mediator_model,
    outcome = outcome_model,
    mediation = mediation_fit
  )
}

main_fit <- fit_models(analysis_data, simulations = 2000, seed = 20260753)
main_summary <- summary(main_fit$mediation)
effects <- data.frame(
  Effect = c("ACME (indirect)", "ADE (direct)", "Total effect", "Proportion mediated"),
  Estimate = c(
    main_summary$d.avg,
    main_summary$z.avg,
    main_summary$tau.coef,
    main_summary$n.avg
  ),
  Lower = c(
    main_summary$d.avg.ci[1],
    main_summary$z.avg.ci[1],
    main_summary$tau.ci[1],
    main_summary$n.avg.ci[1]
  ),
  Upper = c(
    main_summary$d.avg.ci[2],
    main_summary$z.avg.ci[2],
    main_summary$tau.ci[2],
    main_summary$n.avg.ci[2]
  ),
  PValue = c(
    main_summary$d.avg.p,
    main_summary$z.avg.p,
    main_summary$tau.p,
    main_summary$n.avg.p
  )
)

set.seed(20260753)
sensitivity_fit <- mediation::medsens(
  main_fit$mediation,
  rho.by = 0.05,
  effect.type = "indirect"
)
sensitivity <- data.frame(
  Rho = sensitivity_fit$rho,
  ACME = (sensitivity_fit$d0 + sensitivity_fit$d1) / 2,
  Lower = (sensitivity_fit$lower.d0 + sensitivity_fit$lower.d1) / 2,
  Upper = (sensitivity_fit$upper.d0 + sensitivity_fit$upper.d1) / 2
)

pseudocounts <- c(0.25, 0.5, 1.0)
transformation_sensitivity <- do.call(rbind, lapply(
  seq_along(pseudocounts),
  function(index) {
    value <- pseudocounts[index]
    data_variant <- analysis_data
    data_variant$FaecalibacteriumCLR <- make_mediator(value)
    fit <- fit_models(
      data_variant,
      simulations = 1000,
      seed = 20260753 + index
    )
    fit_summary <- summary(fit$mediation)
    data.frame(
      Pseudocount = value,
      ACME = fit_summary$d.avg,
      Lower = fit_summary$d.avg.ci[1],
      Upper = fit_summary$d.avg.ci[2],
      PValue = fit_summary$d.avg.p
    )
  }
))

extract_model <- function(model, model_name) {
  value <- summary(model)$coefficients
  critical <- stats::qt(0.975, df = stats::df.residual(model))
  data.frame(
    Model = model_name,
    Term = rownames(value),
    Estimate = value[, "Estimate"],
    SE = value[, "Std. Error"],
    Lower = value[, "Estimate"] - critical * value[, "Std. Error"],
    Upper = value[, "Estimate"] + critical * value[, "Std. Error"],
    PValue = value[, "Pr(>|t|)"],
    row.names = NULL
  )
}
regression_coefficients <- rbind(
  extract_model(main_fit$mediator, "Mediator model"),
  extract_model(main_fit$outcome, "Outcome model")
)

write.table(
  analysis_data,
  file.path(out_dir, "analysis-data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (object in c(
  "effects", "sensitivity", "transformation_sensitivity",
  "regression_coefficients"
)) {
  write.table(
    get(object),
    file.path(out_dir, paste0(gsub("_", "-", object), ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

path_a <- summary(main_fit$mediator)$coefficients[
  "SteroidExposureYes", "Estimate"
]
path_a_p <- summary(main_fit$mediator)$coefficients[
  "SteroidExposureYes", "Pr(>|t|)"
]
path_b <- summary(main_fit$outcome)$coefficients[
  "FaecalibacteriumCLR", "Estimate"
]
path_b_p <- summary(main_fit$outcome)$coefficients[
  "FaecalibacteriumCLR", "Pr(>|t|)"
]
path_direct <- summary(main_fit$outcome)$coefficients[
  "SteroidExposureYes", "Estimate"
]
path_direct_p <- summary(main_fit$outcome)$coefficients[
  "SteroidExposureYes", "Pr(>|t|)"
]

nodes <- data.frame(
  X = c(0, 1.6, 3.2),
  Y = c(1, 1.55, 1),
  Label = c(
    "Exposure\nCurrent steroid use",
    "Mediator\nFaecalibacterium CLR",
    "Outcome\nlog Calprotectin"
  )
)
edges <- data.frame(
  X = c(0.28, 1.87, 0.32),
  Y = c(1.12, 1.42, 0.92),
  XEnd = c(1.32, 2.92, 2.88),
  YEnd = c(1.43, 1.12, 0.92),
  LabelX = c(0.8, 2.4, 1.6),
  LabelY = c(1.48, 1.48, 0.78),
  Label = c(
    sprintf("a = %.3f\nP = %.3f", path_a, path_a_p),
    sprintf("b = %.3f\nP = %.3f", path_b, path_b_p),
    sprintf("c' = %.3f\nP = %.3f", path_direct, path_direct_p)
  ),
  Role = c("Indirect path", "Indirect path", "Direct path")
)
figure_1 <- ggplot() +
  geom_segment(
    data = edges,
    aes(X, Y, xend = XEnd, yend = YEnd, colour = Role),
    linewidth = 1.05,
    arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")
  ) +
  geom_label(
    data = nodes,
    aes(X, Y, label = Label),
    size = 3.5,
    lineheight = 1.05,
    fill = "white",
    label.size = 0.35
  ) +
  geom_label(
    data = edges,
    aes(LabelX, LabelY, label = Label),
    size = 3,
    lineheight = 1,
    fill = "#F7F7F7",
    label.size = 0.25,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 1.6,
    y = 0.48,
    label = sprintf(
      "ACME = %.3f (95%% CI %.3f to %.3f)",
      main_summary$d.avg,
      main_summary$d.avg.ci[1],
      main_summary$d.avg.ci[2]
    ),
    size = 3.3,
    colour = "#444444"
  ) +
  scale_colour_manual(values = c(
    "Indirect path" = pal_pub[["blue"]],
    "Direct path" = pal_pub[["vermillion"]]
  )) +
  coord_cartesian(xlim = c(-0.45, 3.65), ylim = c(0.3, 1.95), clip = "off") +
  labs(
    title = "Cross-sectional mediation decomposes a model, not time",
    subtitle = "Adjusted for disease group, age and three concomitant medication indicators",
    colour = "Path"
  ) +
  theme_void(base_family = font_pub) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    plot.margin = margin(10, 15, 10, 15)
  )

plot_effects <- effects[effects$Effect != "Proportion mediated", , drop = FALSE]
plot_effects$Effect <- factor(
  plot_effects$Effect,
  levels = rev(c("ACME (indirect)", "ADE (direct)", "Total effect", "Proportion mediated"))
)
figure_2 <- ggplot(plot_effects, aes(Estimate, Effect)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = Effect == "ACME (indirect)"), size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = pal_pub[["grey"]], `TRUE` = pal_pub[["blue"]])) +
  labs(
    title = "None of the decomposed effects excludes the null",
    subtitle = "Quasi-Bayesian intervals from 2,000 simulation draws",
    x = "Effect on log calprotectin scale (95% interval)",
    y = NULL
  ) +
  theme_pub(10) +
  theme(legend.position = "none")

figure_3 <- ggplot(sensitivity, aes(Rho, ACME)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = pal_pub[["sky"]], alpha = 0.25) +
  geom_line(colour = pal_pub[["blue"]], linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#555555") +
  geom_vline(xintercept = 0, linetype = 3, colour = "#777777") +
  labs(
    title = "Residual mediator-outcome confounding can change the ACME",
    subtitle = "rho is the assumed correlation between the two regression error terms",
    x = "Assumed residual correlation (rho)",
    y = "Sensitivity-adjusted ACME"
  ) +
  theme_pub(10)

steroid_colours <- c(No = pal_pub[["blue"]], Yes = pal_pub[["vermillion"]])
observed_a <- ggplot(
  analysis_data,
  aes(SteroidExposure, FaecalibacteriumCLR, fill = SteroidExposure)
) +
  geom_violin(trim = FALSE, alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_point(
    position = position_jitter(width = 0.09),
    size = 0.8,
    alpha = 0.45
  ) +
  scale_fill_manual(values = steroid_colours) +
  labs(
    title = "A. Exposure-to-mediator contrast",
    x = "Current steroid use",
    y = "Faecalibacterium CLR (z-score)",
    fill = "Steroid use"
  ) +
  theme_pub(9)
observed_b <- ggplot(
  analysis_data,
  aes(FaecalibacteriumCLR, LogCalprotectin, colour = SteroidExposure, shape = StudyGroup)
) +
  geom_point(size = 1.6, alpha = 0.72) +
  geom_smooth(aes(group = SteroidExposure), method = "lm", se = TRUE, linewidth = 0.75) +
  scale_colour_manual(values = steroid_colours) +
  labs(
    title = "B. Mediator-to-outcome relation",
    x = "Faecalibacterium CLR (z-score)",
    y = "log(1 + fecal calprotectin)",
    colour = "Steroid use",
    shape = "Disease group"
  ) +
  theme_pub(9)
figure_4 <- observed_a + observed_b + plot_layout(widths = c(0.8, 1.2), guides = "collect") &
  theme(legend.position = "top")

save_pub(figure_1, file.path(figure_dir, "53-1-mediation-dag"), 180, 105)
save_pub(figure_2, file.path(figure_dir, "53-2-effect-decomposition"), 180, 110)
save_pub(figure_3, file.path(figure_dir, "53-3-unmeasured-confounding"), 180, 110)
save_pub(figure_4, file.path(figure_dir, "53-4-observed-data"), 180, 115)

closest_zero <- sensitivity[which.min(abs(sensitivity$ACME)), ]
summary <- list(
  samples = nrow(analysis_data),
  steroid_exposed = sum(analysis_data$SteroidExposure == "Yes"),
  steroid_unexposed = sum(analysis_data$SteroidExposure == "No"),
  mediator = "Faecalibacterium CLR, pseudocount 0.5, z-scored",
  outcome = "log(1 + fecal calprotectin)",
  acme = unname(main_summary$d.avg),
  acme_lower = unname(main_summary$d.avg.ci[1]),
  acme_upper = unname(main_summary$d.avg.ci[2]),
  acme_p = unname(main_summary$d.avg.p),
  ade = unname(main_summary$z.avg),
  total_effect = unname(main_summary$tau.coef),
  sensitivity_rho_nearest_zero = unname(closest_zero$Rho),
  interpretation = "The cross-sectional decomposition does not support a non-zero indirect effect and cannot verify exposure-mediator-outcome time order."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "summary.json")
)

message(
  "Article 53 complete: n = ", nrow(analysis_data),
  "; ACME = ", format(summary$acme, digits = 3),
  ", P = ", format(summary$acme_p, digits = 3), "."
)
