# 约束排序：CAP/db-RDA/RDA/CCA（复现）
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggrepel, scales, vegan.

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

library(vegan)
library(ggplot2)

set.seed(20260718)

pal_pub <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442", "#999999"
)

scale_color_pub <- function(...) {
  ggplot2::scale_color_manual(values = pal_pub, ...)
}

scale_fill_pub <- function(...) {
  ggplot2::scale_fill_manual(values = pal_pub, ...)
}

theme_pub <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7),
        color = NA
      )
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 125,
  height = 100,
  units = "mm",
  dpi = 300
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"),
    plot,
    width = width,
    height = height,
    units = units,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".tiff"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(plot)
}

otutab <- read.delim(
  "data/small/otutab.tsv",
  row.names = 1,
  check.names = FALSE
)
taxonomy <- read.delim(
  "data/small/taxonomy.tsv",
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  "data/small/metadata.tsv",
  row.names = 1,
  check.names = FALSE
)
environment <- read.delim(
  "data/small/environment.tsv",
  row.names = 1,
  check.names = FALSE
)

otutab <- as.matrix(otutab)
storage.mode(otutab) <- "numeric"
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(colnames(otutab) %in% rownames(environment)),
  all(otutab >= 0),
  all(colSums(otutab) > 0)
)

environment <- environment[colnames(otutab), , drop = FALSE]

dim(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)
head(environment)

min_prevalence <- 0.05
min_samples <- ceiling(min_prevalence * ncol(otutab))
keep_feature <- rowSums(otutab > 0) >= min_samples

comm <- t(otutab[keep_feature, , drop = FALSE])
comm_rel <- sweep(comm, 1, rowSums(comm), "/")
metadata <- metadata[rownames(comm_rel), , drop = FALSE]
environment <- environment[rownames(comm_rel), , drop = FALSE]

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
metadata$Group <- factor(
  unname(group_labels[metadata$Group]),
  levels = unname(group_labels)
)

stopifnot(
  all(abs(rowSums(comm_rel) - 1) < 1e-10),
  identical(rownames(comm_rel), rownames(metadata)),
  identical(rownames(comm_rel), rownames(environment)),
  !anyNA(metadata$Group)
)

c(
  samples = nrow(comm_rel),
  retained_features = ncol(comm_rel),
  original_features = nrow(otutab)
)
table(metadata$Group)

cap_group <- vegan::capscale(
  comm_rel ~ Group,
  data = metadata,
  distance = "bray",
  add = "lingoes"
)

set.seed(20260718)
cap_group_test <- anova(
  cap_group,
  permutations = 999
)
cap_group_r2 <- vegan::RsquareAdj(cap_group)

cap_group_test
cap_group_r2

predictors <- c(
  "pH",
  "TOC",
  "NH4",
  "NO3",
  "Conductivity"
)
stopifnot(all(predictors %in% colnames(environment)))

complete_samples <- complete.cases(environment[, predictors, drop = FALSE])
comm_env <- comm_rel[complete_samples, , drop = FALSE]
metadata_env <- metadata[complete_samples, , drop = FALSE]
environment_env <- environment[
  complete_samples,
  predictors,
  drop = FALSE
]

model_data <- as.data.frame(scale(environment_env))
stopifnot(
  identical(rownames(comm_env), rownames(model_data)),
  all(vapply(model_data, is.numeric, logical(1)))
)

summary(model_data)

cap_env <- vegan::capscale(
  comm_env ~ pH + TOC + NH4 + NO3 + Conductivity,
  data = model_data,
  distance = "bray",
  add = "lingoes"
)

cap_env_r2 <- vegan::RsquareAdj(cap_env)
cap_env_vif <- vegan::vif.cca(cap_env)

cap_env_r2
cap_env_vif

set.seed(20260718)
cap_env_overall <- anova(
  cap_env,
  permutations = 999
)

set.seed(20260718)
cap_env_terms <- anova(
  cap_env,
  by = "terms",
  permutations = 999
)

set.seed(20260718)
cap_env_margin <- anova(
  cap_env,
  by = "margin",
  permutations = 999
)

set.seed(20260718)
cap_env_axes <- anova(
  cap_env,
  by = "axis",
  permutations = 999
)

cap_env_overall
cap_env_terms
cap_env_margin
cap_env_axes

cap_env_p <- unname(cap_env_overall$`Pr(>F)`[1])
cap_env_r2_raw <- unname(cap_env_r2$r.squared)
cap_env_r2_adjusted <- unname(cap_env_r2$adj.r.squared)

env_eigenvalues_mod2 <- cap_env$CCA$eig
env_axis_percent_mod2 <- 100 *
  env_eigenvalues_mod2[1:2] /
  sum(env_eigenvalues_mod2)

lab1_mod2 <- sprintf(
  "CAP1 (%.1f%% constrained)",
  env_axis_percent_mod2[1]
)
lab2_mod2 <- sprintf(
  "CAP2 (%.1f%% constrained)",
  env_axis_percent_mod2[2]
)

site_scores_mod2 <- as.data.frame(
  vegan::scores(
    cap_env,
    display = "sites",
    choices = 1:2,
    scaling = 2
  )
)
site_scores_mod2$SampleID <- rownames(site_scores_mod2)
site_scores_mod2$Group <- metadata_env[
  site_scores_mod2$SampleID,
  "Group"
]

arrow_scores_mod2 <- as.data.frame(
  vegan::scores(
    cap_env,
    display = "bp",
    choices = 1:2,
    scaling = 2
  )
)
arrow_scores_mod2$Variable <- rownames(arrow_scores_mod2)

arrow_multiplier_mod2 <- min(
  diff(range(site_scores_mod2$CAP1)) * 0.35 /
    max(abs(arrow_scores_mod2$CAP1)),
  diff(range(site_scores_mod2$CAP2)) * 0.35 /
    max(abs(arrow_scores_mod2$CAP2))
)
arrow_scores_mod2$CAP1 <- arrow_scores_mod2$CAP1 *
  arrow_multiplier_mod2
arrow_scores_mod2$CAP2 <- arrow_scores_mod2$CAP2 *
  arrow_multiplier_mod2

head(site_scores_mod2)
arrow_scores_mod2

format_p <- function(p) {
  if (is.na(p)) {
    return("= NA")
  }
  if (p < 0.001) {
    return("< 0.001")
  }
  sprintf("= %.3f", p)
}

p_cap_environment <- ggplot(
  site_scores_mod2,
  aes(CAP1, CAP2, color = Group)
) +
  stat_ellipse(
    type = "t",
    level = 0.95,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(size = 2.5, alpha = 0.86) +
  geom_segment(
    data = arrow_scores_mod2,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    inherit.aes = FALSE,
    color = "grey20",
    linewidth = 0.55,
    arrow = grid::arrow(
      length = grid::unit(2.2, "mm"),
      type = "closed"
    )
  ) +
  ggrepel::geom_text_repel(
    data = arrow_scores_mod2,
    aes(x = CAP1, y = CAP2, label = Variable),
    inherit.aes = FALSE,
    color = "grey10",
    size = 3.0,
    fontface = "bold",
    seed = 20260718,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  scale_color_pub(name = "Wetland group") +
  labs(
    title = "Environmental CAP",
    subtitle = paste0(
      "R² = ", sprintf("%.3f", cap_env_r2_raw),
      " · adjusted R² = ", sprintf("%.3f", cap_env_r2_adjusted),
      " · permutation P ", format_p(cap_env_p)
    ),
    x = lab1_mod2,
    y = lab2_mod2
  ) +
  coord_equal(clip = "off") +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9.2),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.3)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

p_cap_environment
save_pub(
  p_cap_environment,
  "figures/24-cap-environment",
  width = 150,
  height = 110
)
