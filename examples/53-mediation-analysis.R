# 中介分析：暴露 → 微生物 → 结局
# Run sequentially in a new working directory.
# Required packages: ggplot2, mediation, patchwork, ragg, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
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

.libPaths(c(".r-lib", .libPaths()))
library(mediation)
library(ggplot2)
library(patchwork)

set.seed(20260753)
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
      legend.position = "top", legend.key = ggplot2::element_blank()
    )
}
save_pub <- function(plot, file_base, width = 180, height = 110,
                     units = "mm", dpi = 600) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width,
                  height = height, units = units,
                  device = grDevices::cairo_pdf, family = font_pub, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".svg"), plot, width = width,
                  height = height, units = units, device = svglite::svglite,
                  bg = "white")
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width,
                  height = height, units = units, dpi = dpi,
                  device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".tiff"), plot, width = width,
                  height = height, units = units, dpi = dpi,
                  device = ragg::agg_tiff, compression = "lzw", bg = "white")
  invisible(plot)
}

base <- "data/small/paired-ibd-multiomics"
otutab <- as.matrix(read.delim(
  file.path(base, "otutab.tsv"), row.names = 1, check.names = FALSE
))
taxonomy <- read.delim(
  file.path(base, "taxonomy.tsv"), row.names = 1, check.names = FALSE
)
metadata <- read.delim(
  file.path(base, "metadata.tsv"), row.names = 1, check.names = FALSE
)
storage.mode(otutab) <- "numeric"
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  taxonomy["MB_G005", "Genus"] == "Faecalibacterium"
)

head(otutab[, 1:5])
head(taxonomy)
head(metadata)

complete <- metadata$StudyGroup %in% c("CD", "UC") &
  !is.na(metadata$FecalCalprotectin) &
  metadata$steroids %in% c("No", "Yes") &
  !is.na(metadata$Age) &
  metadata$antibiotic %in% c("No", "Yes") &
  metadata$immunosuppressant %in% c("No", "Yes") &
  metadata$mesalamine %in% c("No", "Yes")
sample_ids <- rownames(metadata)[complete]
stopifnot(length(sample_ids) == 108)
table(metadata[sample_ids, "StudyGroup"])
table(metadata[sample_ids, "steroids"])

z <- function(x) as.numeric(scale(x))
make_mediator <- function(pseudocount) {
  logged <- log(otutab[, sample_ids, drop = FALSE] + pseudocount)
  clr <- sweep(logged, 2, colMeans(logged), "-")
  z(clr["MB_G005", ])
}

analysis_data <- data.frame(
  SampleID = sample_ids,
  SteroidExposure = factor(metadata[sample_ids, "steroids"], levels = c("No", "Yes")),
  StudyGroup = factor(metadata[sample_ids, "StudyGroup"], levels = c("CD", "UC")),
  AgeZ = z(metadata[sample_ids, "Age"]),
  Antibiotic = factor(metadata[sample_ids, "antibiotic"], levels = c("No", "Yes")),
  Immunosuppressant = factor(
    metadata[sample_ids, "immunosuppressant"], levels = c("No", "Yes")
  ),
  Mesalamine = factor(metadata[sample_ids, "mesalamine"], levels = c("No", "Yes")),
  FaecalibacteriumCLR = make_mediator(0.5),
  LogCalprotectin = log1p(metadata[sample_ids, "FecalCalprotectin"])
)
head(analysis_data)

fit_mediation <- function(data, simulations, seed) {
  mediator_model <- lm(
    FaecalibacteriumCLR ~ SteroidExposure + StudyGroup + AgeZ +
      Antibiotic + Immunosuppressant + Mesalamine,
    data = data
  )
  outcome_model <- lm(
    LogCalprotectin ~ SteroidExposure + FaecalibacteriumCLR +
      StudyGroup + AgeZ + Antibiotic + Immunosuppressant + Mesalamine,
    data = data
  )
  set.seed(seed)
  mediation_fit <- mediation::mediate(
    mediator_model, outcome_model,
    treat = "SteroidExposure", mediator = "FaecalibacteriumCLR",
    control.value = "No", treat.value = "Yes",
    sims = simulations, boot = FALSE, robustSE = TRUE
  )
  list(mediator = mediator_model, outcome = outcome_model, mediation = mediation_fit)
}

main_fit <- fit_mediation(analysis_data, simulations = 2000, seed = 20260753)
main_summary <- summary(main_fit$mediation)
main_summary

path_a <- coef(summary(main_fit$mediator))["SteroidExposureYes", "Estimate"]
path_a_p <- coef(summary(main_fit$mediator))["SteroidExposureYes", "Pr(>|t|)"]
path_b <- coef(summary(main_fit$outcome))["FaecalibacteriumCLR", "Estimate"]
path_b_p <- coef(summary(main_fit$outcome))["FaecalibacteriumCLR", "Pr(>|t|)"]
path_direct <- coef(summary(main_fit$outcome))["SteroidExposureYes", "Estimate"]
path_direct_p <- coef(summary(main_fit$outcome))["SteroidExposureYes", "Pr(>|t|)"]

nodes <- data.frame(
  X = c(0, 1.6, 3.2), Y = c(1, 1.55, 1),
  Label = c(
    "Exposure\nCurrent steroid use",
    "Mediator\nFaecalibacterium CLR",
    "Outcome\nlog Calprotectin"
  )
)
edges <- data.frame(
  X = c(0.28, 1.87, 0.32), Y = c(1.12, 1.42, 0.92),
  XEnd = c(1.32, 2.92, 2.88), YEnd = c(1.43, 1.12, 0.92),
  LabelX = c(0.8, 2.4, 1.6), LabelY = c(1.48, 1.48, 0.78),
  Label = c(
    sprintf("a = %.3f\nP = %.3f", path_a, path_a_p),
    sprintf("b = %.3f\nP = %.3f", path_b, path_b_p),
    sprintf("c' = %.3f\nP = %.3f", path_direct, path_direct_p)
  ),
  Role = c("Indirect path", "Indirect path", "Direct path")
)
p_dag <- ggplot() +
  geom_segment(
    data = edges, aes(X, Y, xend = XEnd, yend = YEnd, colour = Role),
    linewidth = 1.05,
    arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")
  ) +
  geom_label(
    data = nodes, aes(X, Y, label = Label), size = 3.5,
    lineheight = 1.05, fill = "white", label.size = 0.35
  ) +
  geom_label(
    data = edges, aes(LabelX, LabelY, label = Label), size = 3,
    lineheight = 1, fill = "#F7F7F7", label.size = 0.25,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = 1.6, y = 0.48,
    label = sprintf("ACME = %.3f (95%% CI %.3f to %.3f)",
                    main_summary$d.avg, main_summary$d.avg.ci[1],
                    main_summary$d.avg.ci[2]),
    size = 3.3, colour = "#444444"
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
    legend.position = "top", plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    plot.margin = margin(10, 15, 10, 15)
  )
p_dag
save_pub(p_dag, "figures/53-mediation-analysis/53-1-mediation-dag", height = 105)

effects <- data.frame(
  Effect = c("ACME (indirect)", "ADE (direct)", "Total effect"),
  Estimate = c(main_summary$d.avg, main_summary$z.avg, main_summary$tau.coef),
  Lower = c(main_summary$d.avg.ci[1], main_summary$z.avg.ci[1], main_summary$tau.ci[1]),
  Upper = c(main_summary$d.avg.ci[2], main_summary$z.avg.ci[2], main_summary$tau.ci[2]),
  PValue = c(main_summary$d.avg.p, main_summary$z.avg.p, main_summary$tau.p)
)
effects$Effect <- factor(
  effects$Effect, levels = rev(c("ACME (indirect)", "ADE (direct)", "Total effect"))
)
p_effects <- ggplot(effects, aes(Estimate, Effect)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0,
                 colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = Effect == "ACME (indirect)"), size = 2.5) +
  scale_colour_manual(values = c(
    `FALSE` = pal_pub[["grey"]], `TRUE` = pal_pub[["blue"]]
  )) +
  labs(
    title = "None of the decomposed effects excludes the null",
    subtitle = "Quasi-Bayesian intervals from 2,000 simulation draws",
    x = "Effect on log calprotectin scale (95% interval)", y = NULL
  ) + theme_pub(10) + theme(legend.position = "none")
p_effects
save_pub(p_effects, "figures/53-mediation-analysis/53-2-effect-decomposition")

set.seed(20260753)
sensitivity_fit <- mediation::medsens(
  main_fit$mediation, rho.by = 0.05, effect.type = "indirect"
)
sensitivity <- data.frame(
  Rho = sensitivity_fit$rho,
  ACME = (sensitivity_fit$d0 + sensitivity_fit$d1) / 2,
  Lower = (sensitivity_fit$lower.d0 + sensitivity_fit$lower.d1) / 2,
  Upper = (sensitivity_fit$upper.d0 + sensitivity_fit$upper.d1) / 2
)
p_sensitivity <- ggplot(sensitivity, aes(Rho, ACME)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper),
              fill = pal_pub[["sky"]], alpha = 0.25) +
  geom_line(colour = pal_pub[["blue"]], linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#555555") +
  geom_vline(xintercept = 0, linetype = 3, colour = "#777777") +
  labs(
    title = "Residual mediator-outcome confounding can change the ACME",
    subtitle = "rho is the assumed correlation between the two regression error terms",
    x = "Assumed residual correlation (rho)", y = "Sensitivity-adjusted ACME"
  ) + theme_pub(10)
p_sensitivity
save_pub(p_sensitivity, "figures/53-mediation-analysis/53-3-unmeasured-confounding")

pseudocounts <- c(0.25, 0.5, 1)
transformation_sensitivity <- do.call(rbind, lapply(seq_along(pseudocounts), function(i) {
  d <- analysis_data
  d$FaecalibacteriumCLR <- make_mediator(pseudocounts[i])
  s <- summary(fit_mediation(d, 1000, 20260753 + i)$mediation)
  data.frame(
    Pseudocount = pseudocounts[i], ACME = s$d.avg,
    Lower = s$d.avg.ci[1], Upper = s$d.avg.ci[2], PValue = s$d.avg.p
  )
}))
transformation_sensitivity

steroid_colours <- c(No = pal_pub[["blue"]], Yes = pal_pub[["vermillion"]])
p_obs_a <- ggplot(
  analysis_data,
  aes(SteroidExposure, FaecalibacteriumCLR, fill = SteroidExposure)
) +
  geom_violin(trim = FALSE, alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_point(position = position_jitter(width = 0.09), size = 0.8, alpha = 0.45) +
  scale_fill_manual(values = steroid_colours) +
  labs(
    title = "A. Exposure-to-mediator contrast", x = "Current steroid use",
    y = "Faecalibacterium CLR (z-score)", fill = "Steroid use"
  ) + theme_pub(9)
p_obs_b <- ggplot(
  analysis_data,
  aes(FaecalibacteriumCLR, LogCalprotectin,
      colour = SteroidExposure, shape = StudyGroup)
) +
  geom_point(size = 1.6, alpha = 0.72) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = FaecalibacteriumCLR, y = LogCalprotectin,
        colour = SteroidExposure, group = SteroidExposure),
    method = "lm", se = TRUE, linewidth = 0.75
  ) +
  scale_colour_manual(values = steroid_colours) +
  labs(
    title = "B. Mediator-to-outcome relation",
    x = "Faecalibacterium CLR (z-score)", y = "log(1 + fecal calprotectin)",
    colour = "Steroid use", shape = "Disease group"
  ) + theme_pub(9)
p_observed <- p_obs_a + p_obs_b +
  patchwork::plot_layout(widths = c(0.8, 1.2), guides = "collect") &
  theme(legend.position = "top")
p_observed
save_pub(p_observed, "figures/53-mediation-analysis/53-4-observed-data",
         height = 115)
