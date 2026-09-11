# 孟德尔随机化：工具变量、水平多效性与敏感性分析
# Run sequentially in a new working directory.
# Required packages: MRPRESSO, TwoSampleMR, ggplot2, ggrepel, patchwork, ragg, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/mr-mibiogen/exposure.tsv",
  "data/small/mr-mibiogen/outcome.tsv",
  "data/small/mr-mibiogen/selection-ledger.tsv",
  "data/small/mr-mibiogen/source-summary.json"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

.libPaths(c(".r-lib", .libPaths()))
library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(20260754)
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
save_pub <- function(plot, file_base, width = 180, height = 115,
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

base <- "data/small/mr-mibiogen"
exposure <- read.delim(
  file.path(base, "exposure.tsv"), check.names = FALSE, stringsAsFactors = FALSE
)
outcome <- read.delim(
  file.path(base, "outcome.tsv"), check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(
  nrow(exposure) == 6,
  identical(exposure$SNP, outcome$SNP),
  length(unique(exposure$chr.exposure)) == nrow(exposure),
  all(exposure$FStatistic > 10)
)
exposure
outcome

harmonised <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
harmonisation_ledger <- harmonised[, c(
  "SNP", "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome",
  "beta.exposure", "beta.outcome", "palindromic", "ambiguous",
  "mr_keep", "remove", "action"
)]
harmonisation_ledger
analysis_data <- harmonised[harmonised$mr_keep, , drop = FALSE]
stopifnot(nrow(analysis_data) == 5)

methods <- c(
  "mr_ivw", "mr_egger_regression", "mr_weighted_median",
  "mr_simple_mode", "mr_weighted_mode"
)
set.seed(20260754)
mr_results <- TwoSampleMR::mr(analysis_data, method_list = methods)
mr_results$Lower <- mr_results$b - qnorm(0.975) * mr_results$se
mr_results$Upper <- mr_results$b + qnorm(0.975) * mr_results$se
mr_results$OR <- exp(mr_results$b)
mr_results$ORLower <- exp(mr_results$Lower)
mr_results$ORUpper <- exp(mr_results$Upper)
mr_results[, c("method", "nsnp", "b", "se", "pval", "Lower", "Upper", "OR")]

ivw <- mr_results[mr_results$method == "Inverse variance weighted", ]
egger <- mr_results[mr_results$method == "MR Egger", ]
median <- mr_results[mr_results$method == "Weighted median", ]
egger_intercept <- TwoSampleMR::mr_pleiotropy_test(analysis_data)
line_data <- data.frame(
  Method = c("IVW", "MR Egger", "Weighted median"),
  Intercept = c(0, egger_intercept$egger_intercept, 0),
  Slope = c(ivw$b, egger$b, median$b)
)
method_colours <- c(
  IVW = pal_pub[["blue"]],
  `MR Egger` = pal_pub[["vermillion"]],
  `Weighted median` = pal_pub[["green"]]
)
p_scatter <- ggplot(
  analysis_data, aes(beta.exposure, beta.outcome, label = SNP)
) +
  geom_errorbar(
    aes(ymin = beta.outcome - 1.96 * se.outcome,
        ymax = beta.outcome + 1.96 * se.outcome),
    width = 0, colour = "#888888"
  ) +
  geom_errorbarh(
    aes(xmin = beta.exposure - 1.96 * se.exposure,
        xmax = beta.exposure + 1.96 * se.exposure),
    height = 0, colour = "#888888"
  ) +
  geom_point(size = 2.3, colour = "#202020") +
  ggrepel::geom_text_repel(
    seed = 20260754, size = 2.8, max.overlaps = Inf,
    min.segment.length = 0, box.padding = 0.35, point.padding = 0.25
  ) +
  geom_abline(
    data = line_data,
    aes(intercept = Intercept, slope = Slope, colour = Method), linewidth = 0.9
  ) +
  scale_colour_manual(values = method_colours) +
  coord_cartesian(xlim = range(analysis_data$beta.exposure) * 1.15) +
  labs(
    title = "Harmonised SNP effects do not share a stable causal slope",
    subtitle = "Five non-ambiguous instruments; error bars are 95% confidence intervals",
    x = "SNP effect on Bifidobacterium abundance",
    y = "SNP effect on log-odds of IBD", colour = "MR method"
  ) + theme_pub(10)
p_scatter
save_pub(p_scatter, "figures/54-mendelian-randomization/54-1-harmonised-scatter",
         height = 120)

mr_results$method <- factor(
  mr_results$method,
  levels = rev(c(
    "Inverse variance weighted", "MR Egger", "Weighted median",
    "Simple mode", "Weighted mode"
  ))
)
p_methods <- ggplot(mr_results, aes(b, method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0,
                 colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = method == "Inverse variance weighted"), size = 2.5) +
  scale_colour_manual(values = c(
    `FALSE` = pal_pub[["grey"]], `TRUE` = pal_pub[["blue"]]
  )) +
  labs(
    title = "Sensitivity estimators all include the null",
    subtitle = "Effects are log-odds of IBD per unit higher genetically predicted abundance",
    x = "MR estimate (95% confidence interval)", y = NULL
  ) + theme_pub(10) + theme(legend.position = "none")
p_methods
save_pub(p_methods, "figures/54-mendelian-randomization/54-2-method-forest",
         height = 110)

heterogeneity <- TwoSampleMR::mr_heterogeneity(analysis_data)
egger_intercept <- TwoSampleMR::mr_pleiotropy_test(analysis_data)
heterogeneity
egger_intercept

set.seed(20260754)
presso <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = analysis_data, NbDistribution = 2000, SignifThreshold = 0.05
)
presso[["Main MR results"]]
presso[["MR-PRESSO results"]][["Global Test"]]

leave_one_out <- TwoSampleMR::mr_leaveoneout(analysis_data)
leave_one_out$Lower <- leave_one_out$b - qnorm(0.975) * leave_one_out$se
leave_one_out$Upper <- leave_one_out$b + qnorm(0.975) * leave_one_out$se
leave_one_out$Label <- ifelse(
  leave_one_out$SNP == "All", "All instruments", paste("Remove", leave_one_out$SNP)
)
leave_one_out$Label <- factor(leave_one_out$Label, levels = rev(leave_one_out$Label))
p_loo <- ggplot(leave_one_out, aes(b, Label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0,
                 colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = SNP == "All"), size = 2.3) +
  scale_colour_manual(values = c(
    `FALSE` = pal_pub[["orange"]], `TRUE` = pal_pub[["blue"]]
  )) +
  labs(
    title = "Leave-one-out estimates remain imprecise",
    subtitle = "A stable sign alone would not overcome weak instrument design or pleiotropy",
    x = "IVW estimate after exclusion (95% confidence interval)", y = NULL
  ) + theme_pub(10) + theme(legend.position = "none")
p_loo
save_pub(p_loo, "figures/54-mendelian-randomization/54-3-leave-one-out",
         height = 110)

presso_global_p <- presso[["MR-PRESSO results"]][["Global Test"]]$Pvalue
assumption_audit <- data.frame(
  Assumption = c(
    "Instrument strength", "LD independence", "Allele harmonisation",
    "Directional pleiotropy", "Outlier pleiotropy", "Steiger directionality",
    "Population overlap", "Locus colocalisation"
  ),
  Status = c(
    "Supported", "Partial", "Partial", "Low power",
    "No signal", "Unavailable", "Unresolved", "Unavailable"
  ),
  Evidence = c(
    sprintf("Minimum F = %.1f", min(exposure$FStatistic)),
    "One SNP per chromosome;\nno ancestry-matched clumping",
    sprintf("%d of %d SNPs retained", nrow(analysis_data), nrow(harmonised)),
    sprintf("Egger intercept P = %.3f", egger_intercept$pval),
    sprintf("MR-PRESSO global P = %.3f", presso_global_p),
    "Effect-allele frequencies absent",
    "Cohort overlap unresolved\nfrom compact files",
    "Top-hit file lacks full\nlocus statistics"
  )
)

strength <- merge(
  exposure[, c("SNP", "FStatistic")],
  harmonisation_ledger[, c("SNP", "mr_keep")],
  by = "SNP", all.x = TRUE, sort = FALSE
)
strength$Status <- ifelse(strength$mr_keep, "Retained", "Dropped: palindromic")
p_strength <- ggplot(strength, aes(FStatistic, reorder(SNP, FStatistic), colour = Status)) +
  geom_vline(xintercept = 10, linetype = 2, colour = pal_pub[["vermillion"]]) +
  geom_segment(aes(x = 0, xend = FStatistic, yend = reorder(SNP, FStatistic)),
               linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(
    Retained = pal_pub[["blue"]],
    `Dropped: palindromic` = pal_pub[["orange"]]
  )) +
  labs(
    title = "A. Candidate instruments", subtitle = "The dashed line marks F = 10",
    x = "Approximate F statistic", y = NULL, colour = "Harmonisation"
  ) + theme_pub(9)

assumption_audit$Assumption <- factor(
  assumption_audit$Assumption, levels = rev(assumption_audit$Assumption)
)
status_colours <- c(
  Supported = pal_pub[["green"]], Partial = pal_pub[["orange"]],
  `No signal` = pal_pub[["sky"]], `Low power` = pal_pub[["purple"]],
  Unavailable = pal_pub[["vermillion"]], Unresolved = pal_pub[["grey"]]
)
p_assumptions <- ggplot(assumption_audit, aes(1, Assumption, fill = Status)) +
  geom_tile(width = 0.28, height = 0.7, colour = "white") +
  geom_text(aes(x = 1.2, label = Evidence), hjust = 0, size = 2.35,
            lineheight = 0.9, colour = "#252525") +
  scale_fill_manual(values = status_colours) +
  coord_cartesian(xlim = c(0.82, 2.55), clip = "off") +
  labs(title = "B. Assumption audit", x = NULL, y = NULL, fill = "Status") +
  theme_pub(8) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(), plot.margin = margin(5, 85, 5, 5)
  )
p_audit <- p_strength + p_assumptions +
  patchwork::plot_layout(widths = c(0.85, 1.5), guides = "collect") &
  theme(legend.position = "bottom")
p_audit
save_pub(p_audit, "figures/54-mendelian-randomization/54-4-assumption-audit",
         height = 130)
