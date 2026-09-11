# 结构方程模型：从环境因子到微生物再到表型
# Run sequentially in a new working directory.
# Required packages: ggplot2, patchwork, piecewiseSEM, plspm, ragg, svglite.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/environment.tsv",
  "data/small/metadata.tsv",
  "data/small/otutab.tsv",
  "data/small/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

.libPaths(c(".r-lib", .libPaths()))
library(piecewiseSEM)
library(plspm)
library(ggplot2)
library(patchwork)

set.seed(20260752)
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
save_pub <- function(plot, file_base, width = 180, height = 120,
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

otutab <- as.matrix(read.delim(
  "data/small/otutab.tsv", row.names = 1, check.names = FALSE
))
taxonomy <- read.delim(
  "data/small/taxonomy.tsv", row.names = 1, check.names = FALSE
)
metadata <- read.delim(
  "data/small/metadata.tsv", row.names = 1, check.names = FALSE
)
environment <- read.delim(
  "data/small/environment.tsv", row.names = 1, check.names = FALSE
)
storage.mode(otutab) <- "numeric"
sample_ids <- intersect(
  colnames(otutab), intersect(rownames(metadata), rownames(environment))
)
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  length(sample_ids) == 90,
  all(colSums(otutab[, sample_ids, drop = FALSE]) > 0)
)
otutab <- otutab[, sample_ids, drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]
environment <- environment[sample_ids, , drop = FALSE]

dim(otutab)
head(taxonomy)
head(metadata)
head(environment[, c("pH", "Conductivity", "TOC", "Temperature", "Precipitation")])

relative <- sweep(otutab, 2, colSums(otutab), "/")
shannon <- -colSums(ifelse(relative > 0, relative * log(relative), 0))
z <- function(x) as.numeric(scale(x))

sem_data <- data.frame(
  SampleID = sample_ids,
  Group = factor(metadata$Group, levels = c("CW", "IW", "TW")),
  Shannon = shannon[sample_ids],
  TOC = environment$TOC,
  pH = environment$pH,
  Conductivity = environment$Conductivity,
  Temperature = environment$Temperature,
  Precipitation = environment$Precipitation,
  ShannonZ = z(shannon[sample_ids]),
  TOCZ = z(log1p(environment$TOC)),
  pHZ = z(environment$pH),
  LogECZ = z(log1p(environment$Conductivity)),
  TempZ = z(environment$Temperature),
  PrecipZ = z(environment$Precipitation)
)
summary(sem_data)

main_models <- list(
  acidity = lm(pHZ ~ TempZ + PrecipZ + Group, data = sem_data),
  salinity = lm(LogECZ ~ TempZ + PrecipZ + Group, data = sem_data),
  microbiome = lm(
    ShannonZ ~ pHZ + LogECZ + TempZ + PrecipZ + Group,
    data = sem_data
  ),
  carbon = lm(
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

piecewiseSEM::fisherC(sem_main)
piecewiseSEM::dSep(sem_main)
piecewiseSEM::coefs(sem_main)

extract_paths <- function(model, response, model_name) {
  x <- summary(model)$coefficients
  x <- x[!grepl("^\\(Intercept\\)$|^Group", rownames(x)), , drop = FALSE]
  critical <- qt(0.975, df = df.residual(model))
  data.frame(
    Model = model_name, Response = response, Predictor = rownames(x),
    Estimate = x[, "Estimate"],
    Lower = x[, "Estimate"] - critical * x[, "Std. Error"],
    Upper = x[, "Estimate"] + critical * x[, "Std. Error"],
    PValue = x[, "Pr(>|t|)"], row.names = NULL
  )
}
main_paths <- do.call(rbind, list(
  extract_paths(main_models$acidity, "pH", "Full direct-path DAG"),
  extract_paths(main_models$salinity, "log Conductivity", "Full direct-path DAG"),
  extract_paths(main_models$microbiome, "Shannon", "Full direct-path DAG"),
  extract_paths(main_models$carbon, "log TOC", "Full direct-path DAG")
))
microbe_carbon <- main_paths[
  main_paths$Response == "log TOC" & main_paths$Predictor == "ShannonZ", ]

nodes <- data.frame(
  X = c(0, 1.4, 2.8, 4.2), Y = 1.15,
  Label = c(
    "Climate\nTemperature + precipitation",
    "Soil chemistry\npH + conductivity",
    "Microbiome\nShannon diversity",
    "Phenotype\nSoil organic carbon"
  )
)
edges <- data.frame(
  X = c(0.38, 0.35, 1.78, 1.75, 3.18),
  XEnd = c(1.02, 2.42, 2.42, 3.82, 3.82),
  Y = c(1.15, 1.28, 1.15, 1.02, 1.15),
  YEnd = c(1.15, 1.28, 1.15, 1.02, 1.15),
  Type = c("Direct", "Direct", "Direct", "Direct", "Focal")
)
p_dag <- ggplot() +
  geom_segment(
    data = edges, aes(X, Y, xend = XEnd, yend = YEnd, colour = Type),
    linewidth = 1,
    arrow = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed")
  ) +
  geom_label(
    data = nodes, aes(X, Y, label = Label), size = 3.3,
    lineheight = 1.05, label.size = 0.35, fill = "white"
  ) +
  annotate(
    "label", x = 3.45, y = 1.52,
    label = sprintf("Shannon -> carbon\nbeta = %.3f, P = %.3f",
                    microbe_carbon$Estimate, microbe_carbon$PValue),
    size = 3.2, fill = "#F7F7F7", label.size = 0.3
  ) +
  annotate(
    "text", x = 2.1, y = 0.68,
    label = "Wetland group adjusted in every component model",
    colour = "#666666", size = 3.1
  ) +
  scale_colour_manual(values = c(
    Direct = pal_pub[["grey"]], Focal = pal_pub[["vermillion"]]
  )) +
  coord_cartesian(xlim = c(-0.55, 4.75), ylim = c(0.55, 1.75), clip = "off") +
  labs(
    title = "A pre-specified path model separates direct and mediated routes",
    subtitle = "A compatible graph is not proof that any arrow is causal",
    colour = "Path role"
  ) +
  theme_void(base_family = font_pub) +
  theme(
    legend.position = "top", plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    plot.margin = margin(10, 15, 10, 15)
  )
p_dag
save_pub(p_dag, "figures/52-structural-equation-model/52-1-prespecified-dag",
         height = 100)

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
p_paths <- ggplot(plot_paths, aes(Estimate, PredictorLabel, colour = PValue < 0.05)) +
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
    x = "Standardised coefficient", y = "Predictor", colour = "P < 0.05"
  ) + theme_pub(9)
p_paths
save_pub(p_paths, "figures/52-structural-equation-model/52-2-path-coefficients",
         height = 135)

reduced_models <- list(
  acidity = main_models$acidity,
  salinity = main_models$salinity,
  microbiome = lm(ShannonZ ~ pHZ + LogECZ + Group, data = sem_data),
  carbon = lm(TOCZ ~ ShannonZ + pHZ + LogECZ + Group, data = sem_data)
)
sem_reduced <- piecewiseSEM::psem(
  reduced_models$acidity, reduced_models$salinity,
  reduced_models$microbiome, reduced_models$carbon,
  data = sem_data
)

fit_row <- function(model, label) {
  f <- piecewiseSEM::fisherC(model)
  a <- AIC(model)
  data.frame(
    Model = label, FisherC = f$Fisher.C, DF = f$df, PValue = f$P.Value,
    AIC = a$AIC, Parameters = a$K, Samples = a$n
  )
}
sem_fit <- rbind(
  fit_row(sem_main, "Full direct-path DAG"),
  fit_row(sem_reduced, "Mediation-only DAG")
)
sem_fit

inner <- matrix(
  0, 5, 5,
  dimnames = list(
    c("Climate", "Acidity", "Salinity", "Microbiome", "Carbon"),
    c("Climate", "Acidity", "Salinity", "Microbiome", "Carbon")
  )
)
inner["Acidity", "Climate"] <- 1
inner["Salinity", "Climate"] <- 1
inner["Microbiome", c("Acidity", "Salinity")] <- 1
inner["Carbon", c("Acidity", "Salinity", "Microbiome")] <- 1
blocks <- list(c("TempZ", "PrecipZ"), "pHZ", "LogECZ", "ShannonZ", "TOCZ")

set.seed(20260752)
pls_fit <- plspm::plspm(
  sem_data, inner, blocks,
  modes = rep("A", 5), scaled = FALSE,
  boot.val = TRUE, br = 2000
)
pls_paths <- data.frame(
  Path = rownames(pls_fit$boot$paths),
  Estimate = pls_fit$boot$paths[, "Original"],
  Lower = pls_fit$boot$paths[, "perc.025"],
  Upper = pls_fit$boot$paths[, "perc.975"],
  row.names = NULL
)
pls_fit$gof
pls_paths

sem_data$PredictedTOCZ <- predict(main_models$carbon)
sem_data$ResidualTOCZ <- residuals(main_models$carbon)
group_colours <- c(
  CW = pal_pub[["blue"]], IW = pal_pub[["orange"]], TW = pal_pub[["green"]]
)
p_diag_a <- ggplot(sem_data, aes(PredictedTOCZ, TOCZ, colour = Group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "A. Observed versus predicted",
    x = "Predicted log TOC (z-score)", y = "Observed log TOC (z-score)",
    colour = "Wetland group"
  ) + coord_equal() + theme_pub(9)
p_diag_b <- ggplot(sem_data, aes(PredictedTOCZ, ResidualTOCZ, colour = Group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(
    inherit.aes = FALSE,
    aes(x = PredictedTOCZ, y = ResidualTOCZ),
    method = "loess", se = FALSE, colour = "#202020", linewidth = 0.7
  ) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "B. Residual structure",
    x = "Predicted log TOC (z-score)", y = "Residual",
    colour = "Wetland group"
  ) + theme_pub(9)
p_diagnostics <- p_diag_a + p_diag_b + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top")
p_diagnostics
save_pub(p_diagnostics,
         "figures/52-structural-equation-model/52-3-model-diagnostics",
         height = 110)

sem_fit$Compatible <- sem_fit$PValue >= 0.05
p_fit <- ggplot(sem_fit, aes(PValue, Model, colour = Compatible)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = pal_pub[["vermillion"]]) +
  geom_segment(aes(x = 0, xend = PValue, yend = Model), linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(limits = c(0, max(0.12, sem_fit$PValue * 1.15))) +
  scale_colour_manual(
    values = c(`FALSE` = pal_pub[["vermillion"]], `TRUE` = pal_pub[["green"]]),
    labels = c(`FALSE` = "Rejected", `TRUE` = "Compatible")
  ) +
  labs(
    title = "A. Directed-separation test", subtitle = "P < 0.05 rejects the graph",
    x = "Fisher's C P value", y = NULL, colour = "Compatible"
  ) + theme_pub(9)
pls_paths$PathLabel <- gsub(" -> ", "  ->  ", pls_paths$Path, fixed = TRUE)
p_pls <- ggplot(pls_paths, aes(Estimate, reorder(PathLabel, Estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0,
                 colour = pal_pub[["grey"]]) +
  geom_point(size = 2, colour = pal_pub[["purple"]]) +
  labs(
    title = "B. PLS-PM estimator audit",
    subtitle = "Percentile intervals from 2,000 sample-level bootstrap draws",
    x = "PLS path coefficient (95% interval)", y = NULL
  ) + theme_pub(9)
p_audit <- p_fit + p_pls + patchwork::plot_layout(widths = c(0.8, 1.35))
p_audit
save_pub(p_audit, "figures/52-structural-equation-model/52-4-model-audit",
         height = 120)
