#!/usr/bin/env Rscript

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(MRPRESSO)
  library(patchwork)
  library(svglite)
  library(TwoSampleMR)
})
source("R/theme_pub.R")

set.seed(20260754)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "mr-mibiogen")
out_dir <- file.path(root, "results", "54-mendelian-randomization")
figure_dir <- file.path(root, "figures", "54-mendelian-randomization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

exposure <- read.delim(
  file.path(data_dir, "exposure.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
outcome <- read.delim(
  file.path(data_dir, "outcome.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selection_ledger <- read.delim(
  file.path(data_dir, "selection-ledger.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(exposure) == 6L,
  identical(exposure$SNP, outcome$SNP),
  length(unique(exposure$chr.exposure)) == nrow(exposure),
  all(exposure$FStatistic > 10)
)

harmonised <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
harmonisation_ledger <- harmonised[, c(
  "SNP", "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome",
  "beta.exposure", "beta.outcome", "palindromic", "ambiguous",
  "mr_keep", "remove", "action"
)]
analysis_data <- harmonised[harmonised$mr_keep, , drop = FALSE]
stopifnot(nrow(analysis_data) == 5L)

method_list <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_simple_mode",
  "mr_weighted_mode"
)
set.seed(20260754)
mr_results <- TwoSampleMR::mr(analysis_data, method_list = method_list)
mr_results$Lower <- mr_results$b - stats::qnorm(0.975) * mr_results$se
mr_results$Upper <- mr_results$b + stats::qnorm(0.975) * mr_results$se
mr_results$OR <- exp(mr_results$b)
mr_results$ORLower <- exp(mr_results$Lower)
mr_results$ORUpper <- exp(mr_results$Upper)
heterogeneity <- TwoSampleMR::mr_heterogeneity(analysis_data)
egger_intercept <- TwoSampleMR::mr_pleiotropy_test(analysis_data)
leave_one_out <- TwoSampleMR::mr_leaveoneout(analysis_data)
leave_one_out$Lower <- leave_one_out$b - stats::qnorm(0.975) * leave_one_out$se
leave_one_out$Upper <- leave_one_out$b + stats::qnorm(0.975) * leave_one_out$se

set.seed(20260754)
presso <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = analysis_data,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)
presso_main <- as.data.frame(presso[["Main MR results"]])
presso_main$Analysis <- rownames(presso_main)
rownames(presso_main) <- NULL
presso_global <- data.frame(
  RSSObserved = unname(presso[["MR-PRESSO results"]][["Global Test"]]$RSSobs),
  PValue = unname(presso[["MR-PRESSO results"]][["Global Test"]]$Pvalue),
  OutliersDetected = !is.null(presso[["MR-PRESSO results"]][["Outlier Test"]])
)

ratio_estimates <- transform(
  analysis_data,
  WaldRatio = beta.outcome / beta.exposure,
  WaldSE = abs(se.outcome / beta.exposure)
)
ratio_estimates$Lower <- ratio_estimates$WaldRatio - stats::qnorm(0.975) * ratio_estimates$WaldSE
ratio_estimates$Upper <- ratio_estimates$WaldRatio + stats::qnorm(0.975) * ratio_estimates$WaldSE

ivw_result <- mr_results[mr_results$method == "Inverse variance weighted", ]
egger_result <- mr_results[mr_results$method == "MR Egger", ]
weighted_median <- mr_results[mr_results$method == "Weighted median", ]

assumption_audit <- data.frame(
  Assumption = c(
    "Instrument strength",
    "LD independence",
    "Allele harmonisation",
    "Directional pleiotropy",
    "Outlier pleiotropy",
    "Steiger directionality",
    "Population overlap",
    "Locus colocalisation"
  ),
  Status = c(
    "Supported",
    "Partial",
    "Partial",
    "Low power",
    "No signal",
    "Unavailable",
    "Unresolved",
    "Unavailable"
  ),
  Evidence = c(
    sprintf("Minimum F = %.1f", min(exposure$FStatistic)),
    "One SNP per chromosome; no ancestry-matched clumping",
    sprintf("%d of %d SNPs retained", nrow(analysis_data), nrow(harmonised)),
    sprintf("Egger intercept P = %.3f", egger_intercept$pval),
    sprintf("MR-PRESSO global P = %.3f", presso_global$PValue),
    "Effect-allele frequencies absent",
    "Cohort overlap cannot be excluded from compact files",
    "Top-hit exposure file lacks full locus statistics"
  ),
  stringsAsFactors = FALSE
)

write.table(
  harmonisation_ledger,
  file.path(out_dir, "harmonisation-ledger.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (object in c(
  "mr_results", "heterogeneity", "egger_intercept", "leave_one_out",
  "presso_main", "presso_global", "ratio_estimates", "assumption_audit"
)) {
  write.table(
    get(object),
    file.path(out_dir, paste0(gsub("_", "-", object), ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

scatter_x <- range(analysis_data$beta.exposure) * 1.15
line_data <- data.frame(
  Method = c("IVW", "MR Egger", "Weighted median"),
  Intercept = c(0, egger_intercept$egger_intercept, 0),
  Slope = c(ivw_result$b, egger_result$b, weighted_median$b)
)
method_colours <- c(
  IVW = pal_pub[["blue"]],
  `MR Egger` = pal_pub[["vermillion"]],
  `Weighted median` = pal_pub[["green"]]
)
figure_1 <- ggplot(
  analysis_data,
  aes(beta.exposure, beta.outcome, label = SNP)
) +
  geom_errorbar(
    aes(ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome),
    width = 0,
    colour = "#888888"
  ) +
  geom_errorbarh(
    aes(xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure),
    height = 0,
    colour = "#888888"
  ) +
  geom_point(size = 2.3, colour = "#202020") +
  ggrepel::geom_text_repel(
    seed = 20260754,
    size = 2.8,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25
  ) +
  geom_abline(
    data = line_data,
    aes(intercept = Intercept, slope = Slope, colour = Method),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = method_colours) +
  coord_cartesian(xlim = scatter_x) +
  labs(
    title = "Harmonised SNP effects do not share a stable causal slope",
    subtitle = "Five non-ambiguous instruments; error bars are 95% confidence intervals",
    x = "SNP effect on Bifidobacterium abundance",
    y = "SNP effect on log-odds of IBD",
    colour = "MR method"
  ) +
  theme_pub(10)

mr_results$method <- factor(
  mr_results$method,
  levels = rev(c(
    "Inverse variance weighted",
    "MR Egger",
    "Weighted median",
    "Simple mode",
    "Weighted mode"
  ))
)
figure_2 <- ggplot(mr_results, aes(b, method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = method == "Inverse variance weighted"), size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = pal_pub[["grey"]], `TRUE` = pal_pub[["blue"]])) +
  labs(
    title = "Sensitivity estimators all include the null",
    subtitle = "Effects are log-odds of IBD per unit higher genetically predicted abundance",
    x = "MR estimate (95% confidence interval)",
    y = NULL
  ) +
  theme_pub(10) +
  theme(legend.position = "none")

leave_one_out$Label <- ifelse(leave_one_out$SNP == "All", "All instruments", paste("Remove", leave_one_out$SNP))
leave_one_out$Label <- factor(leave_one_out$Label, levels = rev(leave_one_out$Label))
figure_3 <- ggplot(leave_one_out, aes(b, Label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, colour = pal_pub[["grey"]]) +
  geom_point(aes(colour = SNP == "All"), size = 2.3) +
  scale_colour_manual(values = c(`FALSE` = pal_pub[["orange"]], `TRUE` = pal_pub[["blue"]])) +
  labs(
    title = "Leave-one-out estimates remain imprecise",
    subtitle = "A stable sign alone would not overcome weak instrument design or pleiotropy",
    x = "IVW estimate after exclusion (95% confidence interval)",
    y = NULL
  ) +
  theme_pub(10) +
  theme(legend.position = "none")

strength_data <- merge(
  exposure[, c("SNP", "FStatistic")],
  harmonisation_ledger[, c("SNP", "mr_keep")],
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)
strength_data$Status <- ifelse(strength_data$mr_keep, "Retained", "Dropped: palindromic")
audit_a <- ggplot(strength_data, aes(FStatistic, reorder(SNP, FStatistic), colour = Status)) +
  geom_vline(xintercept = 10, linetype = 2, colour = pal_pub[["vermillion"]]) +
  geom_segment(aes(x = 0, xend = FStatistic, yend = reorder(SNP, FStatistic)), linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(
    Retained = pal_pub[["blue"]],
    `Dropped: palindromic` = pal_pub[["orange"]]
  )) +
  labs(
    title = "A. Candidate instruments",
    subtitle = "The dashed line marks F = 10",
    x = "Approximate F statistic",
    y = NULL,
    colour = "Harmonisation"
  ) +
  theme_pub(9)

assumption_audit$Assumption <- factor(
  assumption_audit$Assumption,
  levels = rev(assumption_audit$Assumption)
)
assumption_audit$Evidence <- c(
  sprintf("Minimum F = %.1f", min(exposure$FStatistic)),
  "One SNP per chromosome;\nno ancestry-matched clumping",
  sprintf("%d of %d SNPs retained", nrow(analysis_data), nrow(harmonised)),
  sprintf("Egger intercept P = %.3f", egger_intercept$pval),
  sprintf("MR-PRESSO global P = %.3f", presso_global$PValue),
  "Effect-allele frequencies absent",
  "Cohort overlap unresolved\nfrom compact files",
  "Top-hit file lacks full\nlocus statistics"
)
status_colours <- c(
  Supported = pal_pub[["green"]],
  Partial = pal_pub[["orange"]],
  `No signal` = pal_pub[["sky"]],
  `Low power` = pal_pub[["purple"]],
  Unavailable = pal_pub[["vermillion"]],
  Unresolved = pal_pub[["grey"]]
)
audit_b <- ggplot(assumption_audit, aes(1, Assumption, fill = Status)) +
  geom_tile(width = 0.28, height = 0.7, colour = "white") +
  geom_text(
    aes(x = 1.2, label = Evidence),
    hjust = 0,
    size = 2.35,
    lineheight = 0.9,
    colour = "#252525"
  ) +
  scale_fill_manual(values = status_colours) +
  coord_cartesian(xlim = c(0.82, 2.55), clip = "off") +
  labs(
    title = "B. Assumption audit",
    x = NULL,
    y = NULL,
    fill = "Status"
  ) +
  theme_pub(8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 85, 5, 5)
  )
figure_4 <- audit_a + audit_b + plot_layout(widths = c(0.85, 1.5), guides = "collect") &
  theme(legend.position = "bottom")

save_pub(figure_1, file.path(figure_dir, "54-1-harmonised-scatter"), 180, 120)
save_pub(figure_2, file.path(figure_dir, "54-2-method-forest"), 180, 110)
save_pub(figure_3, file.path(figure_dir, "54-3-leave-one-out"), 180, 110)
save_pub(figure_4, file.path(figure_dir, "54-4-assumption-audit"), 180, 130)

summary <- list(
  exposure = unique(exposure$exposure),
  outcome = unique(outcome$outcome),
  candidate_instruments = nrow(harmonised),
  retained_instruments = nrow(analysis_data),
  dropped_palindromic = harmonisation_ledger$SNP[
    harmonisation_ledger$palindromic & !harmonisation_ledger$mr_keep
  ],
  minimum_candidate_f = min(exposure$FStatistic),
  ivw_beta = unname(ivw_result$b),
  ivw_se = unname(ivw_result$se),
  ivw_p = unname(ivw_result$pval),
  ivw_or = unname(ivw_result$OR),
  ivw_heterogeneity_p = unname(
    heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"]
  ),
  egger_intercept = unname(egger_intercept$egger_intercept),
  egger_intercept_p = unname(egger_intercept$pval),
  mr_presso_global_p = unname(presso_global$PValue),
  steiger_status = "Not estimable because compact source records lack effect-allele frequencies",
  colocalisation_status = "Not estimable from a top-hit-only exposure file",
  interpretation = "This compact teaching analysis finds no robust evidence that genetically predicted Bifidobacterium abundance changes IBD risk; design limitations prevent a causal claim."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "summary.json")
)

message(
  "Article 54 complete: ", nrow(analysis_data), " instruments retained; IVW P = ",
  format(summary$ivw_p, digits = 3), "."
)
