#!/usr/bin/env Rscript

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(piecewiseSEM)
  library(plspm)
  library(svglite)
})
source("R/theme_pub.R")

set.seed(20260752)

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small")
out_dir <- file.path(root, "results", "52-structural-equation-model")
figure_dir <- file.path(root, "figures", "52-structural-equation-model")
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
environment <- read.delim(
  file.path(data_dir, "environment.tsv"),
  row.names = 1,
  check.names = FALSE
)
storage.mode(otutab) <- "numeric"
sample_ids <- intersect(
  colnames(otutab),
  intersect(rownames(metadata), rownames(environment))
)
stopifnot(
  length(sample_ids) == 90L,
  identical(rownames(otutab), rownames(taxonomy)),
  all(colSums(otutab[, sample_ids, drop = FALSE]) > 0)
)
otutab <- otutab[, sample_ids, drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
environment <- environment[sample_ids, , drop = FALSE]

relative <- sweep(otutab, 2L, colSums(otutab), "/")
shannon <- -colSums(ifelse(relative > 0, relative * log(relative), 0))
z_score <- function(value) as.numeric(scale(value))
sem_data <- data.frame(
  SampleID = sample_ids,
  Group = factor(metadata$Group, levels = c("CW", "IW", "TW")),
  Shannon = shannon[sample_ids],
  TOC = environment$TOC,
  pH = environment$pH,
  Conductivity = environment$Conductivity,
  Temperature = environment$Temperature,
  Precipitation = environment$Precipitation,
  ShannonZ = z_score(shannon[sample_ids]),
  TOCZ = z_score(log1p(environment$TOC)),
  pHZ = z_score(environment$pH),
  LogECZ = z_score(log1p(environment$Conductivity)),
  TempZ = z_score(environment$Temperature),
  PrecipZ = z_score(environment$Precipitation),
  stringsAsFactors = FALSE
)

# The directed graph is declared explicitly before model fitting. Habitat group
# is an adjustment variable in every component and is not interpreted as a
# manipulable causal exposure.
main_models <- list(
  acidity = stats::lm(pHZ ~ TempZ + PrecipZ + Group, data = sem_data),
  salinity = stats::lm(LogECZ ~ TempZ + PrecipZ + Group, data = sem_data),
  microbiome = stats::lm(
    ShannonZ ~ pHZ + LogECZ + TempZ + PrecipZ + Group,
    data = sem_data
  ),
  carbon = stats::lm(
    TOCZ ~ ShannonZ + pHZ + LogECZ + TempZ + PrecipZ + Group,
    data = sem_data
  )
)
sem_main <- piecewiseSEM::psem(
  main_models$acidity,
  main_models$salinity,
  main_models$microbiome,
  main_models$carbon,
  data = sem_data
)

# A stricter mediation-only graph forbids direct climate-to-microbiome and
# climate-to-carbon paths. It is retained as a falsification/sensitivity model.
reduced_models <- list(
  acidity = main_models$acidity,
  salinity = main_models$salinity,
  microbiome = stats::lm(
    ShannonZ ~ pHZ + LogECZ + Group,
    data = sem_data
  ),
  carbon = stats::lm(
    TOCZ ~ ShannonZ + pHZ + LogECZ + Group,
    data = sem_data
  )
)
sem_reduced <- piecewiseSEM::psem(
  reduced_models$acidity,
  reduced_models$salinity,
  reduced_models$microbiome,
  reduced_models$carbon,
  data = sem_data
)

extract_lm_paths <- function(model, response, model_name) {
  coefficients <- summary(model)$coefficients
  coefficients <- coefficients[
    !grepl("^\\(Intercept\\)$|^Group", rownames(coefficients)),
    ,
    drop = FALSE
  ]
  critical <- stats::qt(0.975, df = stats::df.residual(model))
  data.frame(
    Model = model_name,
    Response = response,
    Predictor = rownames(coefficients),
    Estimate = coefficients[, "Estimate"],
    SE = coefficients[, "Std. Error"],
    Lower = coefficients[, "Estimate"] - critical * coefficients[, "Std. Error"],
    Upper = coefficients[, "Estimate"] + critical * coefficients[, "Std. Error"],
    PValue = coefficients[, "Pr(>|t|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
main_paths <- do.call(rbind, list(
  extract_lm_paths(main_models$acidity, "pH", "Full direct-path DAG"),
  extract_lm_paths(main_models$salinity, "log Conductivity", "Full direct-path DAG"),
  extract_lm_paths(main_models$microbiome, "Shannon", "Full direct-path DAG"),
  extract_lm_paths(main_models$carbon, "log TOC", "Full direct-path DAG")
))
reduced_paths <- do.call(rbind, list(
  extract_lm_paths(reduced_models$acidity, "pH", "Mediation-only DAG"),
  extract_lm_paths(reduced_models$salinity, "log Conductivity", "Mediation-only DAG"),
  extract_lm_paths(reduced_models$microbiome, "Shannon", "Mediation-only DAG"),
  extract_lm_paths(reduced_models$carbon, "log TOC", "Mediation-only DAG")
))
sem_paths <- rbind(main_paths, reduced_paths)

fit_row <- function(model, label) {
  fisher <- piecewiseSEM::fisherC(model)
  information <- AIC(model)
  data.frame(
    Model = label,
    FisherC = fisher$Fisher.C,
    DF = fisher$df,
    PValue = fisher$P.Value,
    AIC = information$AIC,
    Parameters = information$K,
    Samples = information$n
  )
}
sem_fit <- rbind(
  fit_row(sem_main, "Full direct-path DAG"),
  fit_row(sem_reduced, "Mediation-only DAG")
)
directed_separation <- rbind(
  transform(piecewiseSEM::dSep(sem_main), Model = "Full direct-path DAG"),
  transform(piecewiseSEM::dSep(sem_reduced), Model = "Mediation-only DAG")
)

# PLS-PM is used only as an estimator-family sensitivity analysis. The latent
# Climate block contains measured temperature and precipitation; all other
# blocks are single-indicator constructs.
inner <- matrix(
  0,
  nrow = 5,
  ncol = 5,
  dimnames = list(
    c("Climate", "Acidity", "Salinity", "Microbiome", "Carbon"),
    c("Climate", "Acidity", "Salinity", "Microbiome", "Carbon")
  )
)
inner["Acidity", "Climate"] <- 1
inner["Salinity", "Climate"] <- 1
inner["Microbiome", c("Acidity", "Salinity")] <- 1
inner["Carbon", c("Acidity", "Salinity", "Microbiome")] <- 1
blocks <- list(
  c("TempZ", "PrecipZ"),
  "pHZ",
  "LogECZ",
  "ShannonZ",
  "TOCZ"
)
set.seed(20260752)
pls_fit <- plspm::plspm(
  sem_data,
  inner,
  blocks,
  modes = rep("A", 5),
  scaled = FALSE,
  boot.val = TRUE,
  br = 2000
)
pls_paths <- data.frame(
  Path = rownames(pls_fit$boot$paths),
  Estimate = pls_fit$boot$paths[, "Original"],
  BootstrapMean = pls_fit$boot$paths[, "Mean.Boot"],
  BootstrapSE = pls_fit$boot$paths[, "Std.Error"],
  Lower = pls_fit$boot$paths[, "perc.025"],
  Upper = pls_fit$boot$paths[, "perc.975"],
  row.names = NULL,
  check.names = FALSE
)

sem_data$PredictedTOCZ <- stats::predict(main_models$carbon)
sem_data$ResidualTOCZ <- stats::residuals(main_models$carbon)
write.table(
  sem_data,
  file.path(out_dir, "sample-metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (object in c(
  "sem_paths", "sem_fit", "directed_separation", "pls_paths"
)) {
  write.table(
    get(object),
    file.path(out_dir, paste0(gsub("_", "-", object), ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

microbe_carbon <- main_paths[
  main_paths$Response == "log TOC" & main_paths$Predictor == "ShannonZ",
]
node_data <- data.frame(
  Node = c("Climate", "Soil chemistry", "Microbiome", "Phenotype"),
  X = c(0, 1.4, 2.8, 4.2),
  Y = c(1.15, 1.15, 1.15, 1.15),
  Label = c(
    "Climate\nTemperature + precipitation",
    "Soil chemistry\npH + conductivity",
    "Microbiome\nShannon diversity",
    "Phenotype\nSoil organic carbon"
  )
)
edge_data <- data.frame(
  X = c(0.38, 0.35, 1.78, 1.75, 3.18),
  XEnd = c(1.02, 2.42, 2.42, 3.82, 3.82),
  Y = c(1.15, 1.28, 1.15, 1.02, 1.15),
  YEnd = c(1.15, 1.28, 1.15, 1.02, 1.15),
  Type = c("Direct", "Direct", "Direct", "Direct", "Focal")
)
figure_1 <- ggplot() +
  geom_segment(
    data = edge_data,
    aes(X, Y, xend = XEnd, yend = YEnd, colour = Type),
    linewidth = 1.0,
    arrow = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed")
  ) +
  geom_label(
    data = node_data,
    aes(X, Y, label = Label),
    size = 3.3,
    lineheight = 1.05,
    label.size = 0.35,
    fill = "white"
  ) +
  annotate(
    "label",
    x = 3.45,
    y = 1.52,
    label = sprintf(
      "Shannon -> carbon\nbeta = %.3f, P = %.3f",
      microbe_carbon$Estimate,
      microbe_carbon$PValue
    ),
    size = 3.2,
    fill = "#F7F7F7",
    label.size = 0.3
  ) +
  annotate(
    "text",
    x = 2.1,
    y = 0.68,
    label = "Wetland group adjusted in every component model",
    colour = "#666666",
    size = 3.1
  ) +
  scale_colour_manual(values = c(
    Direct = pal_pub[["grey"]],
    Focal = pal_pub[["vermillion"]]
  )) +
  coord_cartesian(xlim = c(-0.55, 4.75), ylim = c(0.55, 1.75), clip = "off") +
  labs(
    title = "A pre-specified path model separates direct and mediated routes",
    subtitle = "A compatible graph is not proof that any arrow is causal",
    colour = "Path role"
  ) +
  theme_void(base_family = font_pub) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    plot.margin = margin(10, 15, 10, 15)
  )

plot_paths <- main_paths
plot_paths$Response <- factor(
  plot_paths$Response,
  levels = c("pH", "log Conductivity", "Shannon", "log TOC")
)
plot_paths$PredictorLabel <- factor(
  plot_paths$Predictor,
  levels = rev(c("TempZ", "PrecipZ", "pHZ", "LogECZ", "ShannonZ")),
  labels = rev(c("Temperature", "Precipitation", "pH", "log Conductivity", "Shannon"))
)
figure_2 <- ggplot(
  plot_paths,
  aes(Estimate, PredictorLabel, colour = PValue < 0.05)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0) +
  geom_point(size = 2) +
  facet_wrap(~ Response, scales = "free_y", ncol = 2) +
  scale_colour_manual(
    values = c(`FALSE` = pal_pub[["grey"]], `TRUE` = pal_pub[["blue"]]),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  labs(
    title = "Standardised component-model paths",
    subtitle = "Intervals are 95% confidence intervals; habitat group coefficients are omitted",
    x = "Standardised coefficient",
    y = "Predictor",
    colour = "P < 0.05"
  ) +
  theme_pub(9)

group_colours <- c(CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]])
diagnostic_a <- ggplot(sem_data, aes(PredictedTOCZ, TOCZ, colour = Group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "A. Observed versus predicted",
    x = "Predicted log TOC (z-score)",
    y = "Observed log TOC (z-score)",
    colour = "Wetland group"
  ) +
  coord_equal() +
  theme_pub(9)
diagnostic_b <- ggplot(sem_data, aes(PredictedTOCZ, ResidualTOCZ, colour = Group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, colour = "#202020", linewidth = 0.7) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "B. Residual structure",
    x = "Predicted log TOC (z-score)",
    y = "Residual",
    colour = "Wetland group"
  ) +
  theme_pub(9)
figure_3 <- diagnostic_a + diagnostic_b + plot_layout(guides = "collect") &
  theme(legend.position = "top")

fit_plot_data <- sem_fit
fit_plot_data$Compatible <- fit_plot_data$PValue >= 0.05
audit_a <- ggplot(fit_plot_data, aes(PValue, Model, colour = Compatible)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = pal_pub[["vermillion"]]) +
  geom_segment(aes(x = 0, xend = PValue, yend = Model), linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(limits = c(0, max(0.12, fit_plot_data$PValue * 1.15))) +
  scale_colour_manual(
    values = c(`FALSE` = pal_pub[["vermillion"]], `TRUE` = pal_pub[["green"]]),
    labels = c(`FALSE` = "Rejected", `TRUE` = "Compatible")
  ) +
  labs(
    title = "A. Directed-separation test",
    subtitle = "P < 0.05 rejects the graph",
    x = "Fisher's C P value",
    y = NULL,
    colour = "Compatible"
  ) +
  theme_pub(9)

pls_paths$PathLabel <- gsub(" -> ", "  ->  ", pls_paths$Path, fixed = TRUE)
audit_b <- ggplot(pls_paths, aes(Estimate, reorder(PathLabel, Estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, colour = pal_pub[["grey"]]) +
  geom_point(size = 2, colour = pal_pub[["purple"]]) +
  labs(
    title = "B. PLS-PM estimator audit",
    subtitle = "Percentile intervals from 2,000 sample-level bootstrap draws",
    x = "PLS path coefficient (95% interval)",
    y = NULL
  ) +
  theme_pub(9)
figure_4 <- audit_a + audit_b + plot_layout(widths = c(0.8, 1.35))

save_pub(figure_1, file.path(figure_dir, "52-1-prespecified-dag"), 180, 100)
save_pub(figure_2, file.path(figure_dir, "52-2-path-coefficients"), 180, 135)
save_pub(figure_3, file.path(figure_dir, "52-3-model-diagnostics"), 180, 110)
save_pub(figure_4, file.path(figure_dir, "52-4-model-audit"), 180, 120)

summary <- list(
  samples = nrow(sem_data),
  measured_environment_variables = c(
    "pH", "conductivity", "temperature", "precipitation"
  ),
  measured_phenotype = "total organic carbon (TOC)",
  full_dag_fisher_c = unname(sem_fit$FisherC[sem_fit$Model == "Full direct-path DAG"]),
  full_dag_p = unname(sem_fit$PValue[sem_fit$Model == "Full direct-path DAG"]),
  reduced_dag_fisher_c = unname(sem_fit$FisherC[sem_fit$Model == "Mediation-only DAG"]),
  reduced_dag_p = unname(sem_fit$PValue[sem_fit$Model == "Mediation-only DAG"]),
  shannon_to_carbon_beta = unname(microbe_carbon$Estimate),
  shannon_to_carbon_p = unname(microbe_carbon$PValue),
  plspm_gof = unname(pls_fit$gof),
  interpretation = "The full graph is compatible with the observed covariance structure, whereas the mediation-only graph is rejected; neither result establishes causality."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "summary.json")
)

message(
  "Article 52 complete: full-DAG Fisher C P = ",
  format(summary$full_dag_p, digits = 3),
  "; Shannon-to-TOC P = ",
  format(summary$shannon_to_carbon_p, digits = 3),
  "."
)
