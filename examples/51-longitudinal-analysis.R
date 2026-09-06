# 纵向：线性混合模型/轨迹/volatility
# Run sequentially in a new working directory.
# Required packages: emmeans, ggplot2, nlme, patchwork, ragg, svglite, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
  "data/small/longitudinal-dietswap/metadata.tsv",
  "data/small/longitudinal-dietswap/otutab.tsv",
  "data/small/longitudinal-dietswap/source-summary.json",
  "data/small/longitudinal-dietswap/taxonomy.tsv"
)
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

.libPaths(c(".r-lib", .libPaths()))
library(nlme)
library(emmeans)
library(vegan)
library(ggplot2)
library(patchwork)

set.seed(20260751)
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
      legend.position = "top",
      legend.key = ggplot2::element_blank()
    )
}
save_pub <- function(plot, file_base, width = 180, height = 125,
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

base <- "data/small/longitudinal-dietswap"
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
  all(otutab >= 0), all(colSums(otutab) > 0),
  !anyDuplicated(paste(metadata$SubjectID, metadata$Timepoint))
)

dim(otutab)
head(taxonomy[, c("Phylum", "Family", "Genus", "DisplayName")])
head(metadata)
with(metadata, table(Cohort, Timepoint))

relative <- sweep(otutab, 2, colSums(otutab), "/")
metadata$Shannon <- -colSums(ifelse(
  relative > 0, relative * log(relative), 0
))[rownames(metadata)]
metadata$Cohort <- factor(
  metadata$Cohort,
  levels = c("African American", "Rural African")
)
metadata$TimeFactor <- factor(metadata$Timepoint, levels = 1:6)
metadata$TimeLabel <- factor(
  paste0("T", metadata$Timepoint), levels = paste0("T", 1:6)
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

model_comparison <- anova(fit_ri_ml, fit_ar1_ml)
model_comparison
unname(coef(fit_ar1$modelStruct$corStruct, unconstrained = FALSE))
summary(fit_ri)$tTable

model_means <- as.data.frame(emmeans::emmeans(
  fit_ri, ~ Cohort * TimeFactor
))
model_means$Timepoint <- as.integer(as.character(model_means$TimeFactor))

cohort_colours <- c(
  "African American" = pal_pub[["blue"]],
  "Rural African" = pal_pub[["orange"]]
)
p_alpha <- ggplot(metadata, aes(Timepoint, Shannon, colour = Cohort)) +
  geom_line(aes(group = SubjectID), linewidth = 0.35, alpha = 0.17) +
  geom_point(size = 0.9, alpha = 0.25) +
  geom_ribbon(
    data = model_means,
    aes(x = Timepoint, ymin = lower.CL, ymax = upper.CL, fill = Cohort),
    inherit.aes = FALSE, alpha = 0.16, colour = NA
  ) +
  geom_line(
    data = model_means,
    aes(Timepoint, emmean, colour = Cohort, group = Cohort),
    linewidth = 1.05
  ) +
  geom_point(
    data = model_means, aes(Timepoint, emmean, colour = Cohort), size = 2.1
  ) +
  scale_x_continuous(breaks = 1:6, labels = paste0("T", 1:6)) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Within-subject diversity trajectories",
    subtitle = "Thin lines: participants; thick lines and bands: mixed-model means and 95% CI",
    x = "Study timepoint", y = "Shannon diversity",
    colour = "Cohort", fill = "Cohort"
  ) +
  theme_pub(10)
p_alpha
save_pub(p_alpha, "figures/51-longitudinal-analysis/51-1-alpha-trajectories")

bray <- vegan::vegdist(t(relative), method = "bray")
pcoa <- stats::cmdscale(bray, k = 5, eig = TRUE, add = TRUE)
pcoa_scores <- data.frame(
  SampleID = rownames(pcoa$points),
  Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2]
)
pcoa_scores <- cbind(
  pcoa_scores,
  metadata[pcoa_scores$SampleID,
           c("SubjectID", "Cohort", "Timepoint", "TimeLabel"), drop = FALSE]
)
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
axis_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)
centroids <- aggregate(
  cbind(Axis1, Axis2) ~ Cohort + Timepoint, pcoa_scores, mean
)

p_pcoa <- ggplot(pcoa_scores, aes(Axis1, Axis2, colour = Cohort)) +
  geom_path(
    aes(group = SubjectID), linewidth = 0.3, alpha = 0.13,
    arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  geom_point(alpha = 0.18, size = 0.8) +
  geom_path(
    data = centroids, aes(group = Cohort), linewidth = 1.2,
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
  ) +
  geom_point(
    data = centroids, aes(fill = Cohort), shape = 21, size = 3, colour = "white"
  ) +
  geom_text(
    data = centroids[centroids$Timepoint %in% c(1, 4, 6), ],
    aes(label = paste0("T", Timepoint)), colour = "#202020",
    nudge_y = 0.008, size = 2.8, show.legend = FALSE
  ) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Bray-Curtis trajectories retain subject history",
    subtitle = "Faint arrows show participants; labelled paths show cohort centroids",
    x = sprintf("PCoA 1 (%.1f%%)", axis_percent[1]),
    y = sprintf("PCoA 2 (%.1f%%)", axis_percent[2]),
    colour = "Cohort", fill = "Cohort"
  ) +
  coord_equal() + theme_pub(10)
p_pcoa
save_pub(p_pcoa, "figures/51-longitudinal-analysis/51-2-pcoa-trajectories",
         height = 135)

bray_matrix <- as.matrix(bray)
volatility <- do.call(rbind, lapply(split(metadata, metadata$SubjectID), function(x) {
  x <- x[order(x$Timepoint), , drop = FALSE]
  if (nrow(x) < 2) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(x) - 1), function(j) {
    if (x$Timepoint[j + 1] - x$Timepoint[j] != 1) return(NULL)
    from <- rownames(x)[j]
    to <- rownames(x)[j + 1]
    data.frame(
      SubjectID = x$SubjectID[j], Cohort = x$Cohort[j],
      Window = paste0("T", x$Timepoint[j], "-T", x$Timepoint[j + 1]),
      Volatility = unname(bray_matrix[from, to])
    )
  }))
}))
volatility$Cohort <- factor(volatility$Cohort, levels = levels(metadata$Cohort))
volatility$Window <- factor(
  volatility$Window, levels = paste0("T", 1:5, "-T", 2:6)
)
fit_volatility <- nlme::lme(
  Volatility ~ Cohort * Window,
  random = ~ 1 | SubjectID,
  data = volatility,
  method = "REML",
  control = nlme::lmeControl(opt = "optim")
)
head(volatility)
summary(fit_volatility)$tTable

p_volatility <- ggplot(volatility, aes(Window, Volatility, fill = Cohort)) +
  geom_boxplot(
    position = position_dodge(width = 0.72), width = 0.62,
    outlier.shape = NA, alpha = 0.72
  ) +
  geom_point(
    aes(colour = Cohort),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72),
    size = 0.8, alpha = 0.45
  ) +
  scale_colour_manual(values = cohort_colours) +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "Volatility is a within-person distance",
    subtitle = "Bray-Curtis distance between consecutive observed timepoints",
    x = "Consecutive time window", y = "Bray-Curtis volatility",
    colour = "Cohort", fill = "Cohort"
  ) +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
p_volatility
save_pub(p_volatility, "figures/51-longitudinal-analysis/51-3-volatility")

extract_fixed <- function(model, label) {
  x <- summary(model)$tTable
  critical <- qt(0.975, df = x[, "DF"])
  data.frame(
    Model = label, Term = rownames(x), Estimate = x[, "Value"],
    Lower = x[, "Value"] - critical * x[, "Std.Error"],
    Upper = x[, "Value"] + critical * x[, "Std.Error"],
    PValue = x[, "p-value"], row.names = NULL
  )
}
fixed_effects <- rbind(
  extract_fixed(fit_ri, "Random intercept"),
  extract_fixed(fit_ar1, "Random intercept + AR(1)")
)
interactions <- fixed_effects[grepl("CohortRural African:TimeFactor", fixed_effects$Term), ]
interactions$Timepoint <- sub(
  "CohortRural African:TimeFactor", "T", interactions$Term, fixed = TRUE
)

paired_change <- do.call(rbind, lapply(split(metadata, metadata$SubjectID), function(x) {
  if (!all(c(1, 4) %in% x$Timepoint)) return(NULL)
  data.frame(
    SubjectID = x$SubjectID[1], Cohort = x$Cohort[1],
    DeltaT4T1 = x$Shannon[x$Timepoint == 4][1] - x$Shannon[x$Timepoint == 1][1]
  )
}))
set.seed(20260751)
bootstrap_draws <- do.call(rbind, lapply(levels(metadata$Cohort), function(g) {
  v <- paired_change$DeltaT4T1[paired_change$Cohort == g]
  data.frame(
    Cohort = g,
    MeanChange = replicate(5000, mean(sample(v, length(v), replace = TRUE)))
  )
}))

p_audit_a <- ggplot(interactions, aes(Estimate, Timepoint, colour = Model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper),
    position = position_dodge(width = 0.45), height = 0
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8) +
  scale_colour_manual(values = c(
    "Random intercept" = pal_pub[["blue"]],
    "Random intercept + AR(1)" = pal_pub[["vermillion"]]
  )) +
  labs(
    title = "A. Correlation-structure audit",
    x = "Cohort-by-time interaction estimate (95% CI)",
    y = "Timepoint vs T1", colour = "Model"
  ) + theme_pub(9)

p_audit_b <- ggplot(bootstrap_draws, aes(MeanChange, fill = Cohort)) +
  geom_density(alpha = 0.45, colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#555555") +
  scale_fill_manual(values = cohort_colours) +
  labs(
    title = "B. Participant-level bootstrap",
    subtitle = "Mean Shannon change, T4 minus T1",
    x = "Bootstrapped mean change", y = "Density", fill = "Cohort"
  ) + theme_pub(9)

p_audit <- p_audit_a + p_audit_b + patchwork::plot_layout(widths = c(1.25, 1))
p_audit
save_pub(p_audit, "figures/51-longitudinal-analysis/51-4-model-audit",
         height = 115)
