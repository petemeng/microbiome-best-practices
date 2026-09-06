# 实验设计：样本量、统计效力、分组与配对
# Run sequentially in a new working directory.
# Required packages: ggplot2, lme4, permute, scales, vegan.

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- c(
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
  width = 120,
  height = 95,
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

otutab <- as.matrix(otutab)
storage.mode(otutab) <- "numeric"

stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  all(colSums(otutab) > 0),
  all(c("Group", "Type", "Saline") %in% colnames(metadata))
)

dim(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)

metadata$GroupLabel <- unname(group_labels[metadata$Group])
metadata$Salinity <- factor(
  metadata$Saline,
  levels = c("Non-saline soil", "Saline soil")
)
stopifnot(!anyNA(metadata$GroupLabel), !anyNA(metadata$Salinity))

design_counts <- as.data.frame(
  with(metadata, table(GroupLabel, Salinity)),
  responseName = "Count"
)
design_counts$Observed <- factor(
  design_counts$Count > 0,
  levels = c(FALSE, TRUE),
  labels = c("Empty design cell", "Observed samples")
)
design_counts$CellLabel <- ifelse(
  design_counts$Count == 0,
  "Not\nobserved",
  paste0("n = ", design_counts$Count)
)

design_matrix <- xtabs(Count ~ GroupLabel + Salinity, design_counts)
type_matrix <- with(metadata, table(Group, Type))

design_matrix
type_matrix

perfect_group_salinity_confounding <- all(
  rowSums(design_matrix > 0) == 1
)
type_nested_within_group <- all(
  colSums(type_matrix > 0) == 1
)

stopifnot(
  perfect_group_salinity_confounding,
  type_nested_within_group
)

library_sizes <- colSums(otutab)
rarefaction_depth <- min(library_sizes)

set.seed(20260718)
otutab_rarefied <- t(
  vegan::rrarefy(
    t(otutab),
    sample = rarefaction_depth
  )
)

stopifnot(
  all(colSums(otutab_rarefied) == rarefaction_depth),
  identical(colnames(otutab_rarefied), rownames(metadata))
)

shannon <- vegan::diversity(
  t(otutab_rarefied),
  index = "shannon"
)

alpha_df <- data.frame(
  SampleID = names(shannon),
  Shannon = as.numeric(shannon),
  Salinity = metadata[names(shannon), "Salinity"],
  Group = metadata[names(shannon), "GroupLabel"],
  stringsAsFactors = FALSE
)

alpha_summary <- do.call(
  rbind,
  lapply(split(alpha_df$Shannon, alpha_df$Salinity), function(x) {
    data.frame(
      N = length(x),
      Mean = mean(x),
      SD = stats::sd(x)
    )
  })
)

alpha_summary

cohen_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(
    ((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) /
      (nx + ny - 2)
  )
  (mean(x) - mean(y)) / pooled_sd
}

shannon_non_saline <- subset(
  alpha_df,
  Salinity == "Non-saline soil"
)$Shannon
shannon_saline <- subset(
  alpha_df,
  Salinity == "Saline soil"
)$Shannon

pilot_d <- cohen_d(
  shannon_non_saline,
  shannon_saline
)

set.seed(20260718)
bootstrap_d <- replicate(2000, {
  x_boot <- sample(
    shannon_non_saline,
    length(shannon_non_saline),
    replace = TRUE
  )
  y_boot <- sample(
    shannon_saline,
    length(shannon_saline),
    replace = TRUE
  )
  cohen_d(x_boot, y_boot)
})
bootstrap_d <- bootstrap_d[is.finite(bootstrap_d)]
pilot_ci <- unname(
  stats::quantile(
    bootstrap_d,
    probs = c(0.025, 0.975)
  )
)

pilot_test <- stats::t.test(
  shannon_non_saline,
  shannon_saline,
  var.equal = FALSE
)

data.frame(
  CohenD = pilot_d,
  BootstrapLower = pilot_ci[1],
  BootstrapUpper = pilot_ci[2],
  WelchP = pilot_test$p.value
)

pilot_effect <- abs(pilot_d)
effect_scenarios <- data.frame(
  Scenario = c(
    "Small effect",
    "Pilot estimate",
    "Moderate effect"
  ),
  Effect = c(0.20, pilot_effect, 0.50),
  stringsAsFactors = FALSE
)

n_grid <- seq(10, 450, by = 2)

power_curve <- do.call(
  rbind,
  lapply(seq_len(nrow(effect_scenarios)), function(i) {
    effect <- effect_scenarios$Effect[i]
    data.frame(
      NPerGroup = n_grid,
      Power = vapply(
        n_grid,
        function(n) {
          stats::power.t.test(
            n = n,
            delta = effect,
            sd = 1,
            sig.level = 0.05,
            type = "two.sample",
            alternative = "two.sided"
          )$power
        },
        numeric(1)
      ),
      Scenario = effect_scenarios$Scenario[i],
      Effect = effect
    )
  })
)

n80 <- effect_scenarios
n80$NPerGroup <- vapply(
  n80$Effect,
  function(effect) {
    ceiling(
      stats::power.t.test(
        delta = effect,
        sd = 1,
        sig.level = 0.05,
        power = 0.80,
        type = "two.sample",
        alternative = "two.sided"
      )$n
    )
  },
  numeric(1)
)
n80$Power <- 0.80

current_smaller_group <- min(table(alpha_df$Salinity))
current_pilot_power <- stats::power.t.test(
  n = current_smaller_group,
  delta = pilot_effect,
  sd = 1,
  sig.level = 0.05,
  type = "two.sample",
  alternative = "two.sided"
)$power

n80

design_counts$GroupLabel <- factor(
  design_counts$GroupLabel,
  levels = rev(c(
    "Inland wetland",
    "Coastal wetland",
    "Tibetan Plateau"
  ))
)

design_fill <- c(
  "Empty design cell" = "#F6D6D6",
  "Observed samples" = "#2878A8"
)

p_design_confounding <- ggplot(
  design_counts,
  aes(Salinity, GroupLabel, fill = Observed)
) +
  geom_tile(
    color = "white",
    linewidth = 1.5,
    width = 0.96,
    height = 0.90
  ) +
  geom_text(
    aes(
      label = CellLabel,
      color = Observed
    ),
    size = 4.0,
    fontface = "bold",
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = design_fill,
    name = "Design cell"
  ) +
  scale_color_manual(
    values = c(
      "Empty design cell" = "#8B1A1A",
      "Observed samples" = "white"
    ),
    guide = "none"
  ) +
  labs(
    title = "Salinity is completely confounded with wetland group",
    subtitle = paste0(
      "Only 3 of 6 design cells are observed\n",
      "Geography and salinity cannot be separated"
    ),
    x = "Salinity category",
    y = "Wetland group"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(
      color = "grey35",
      size = 9.2,
      lineheight = 0.95
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p_design_confounding
save_pub(
  p_design_confounding,
  "figures/03-design-confounding",
  width = 172,
  height = 112
)

scenario_labels <- c(
  "Small effect" = "Small effect (d = 0.20)",
  "Pilot estimate" = sprintf(
    "Pilot estimate (d = %.2f)",
    pilot_effect
  ),
  "Moderate effect" = "Moderate effect (d = 0.50)"
)

power_curve$ScenarioLabel <- factor(
  unname(scenario_labels[power_curve$Scenario]),
  levels = unname(scenario_labels)
)
n80$ScenarioLabel <- factor(
  unname(scenario_labels[n80$Scenario]),
  levels = unname(scenario_labels)
)

power_colors <- stats::setNames(
  c("#7A5195", "#0072B2", "#009E73"),
  unname(scenario_labels)
)

p_power_sensitivity <- ggplot(
  power_curve,
  aes(NPerGroup, Power, color = ScenarioLabel)
) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = current_smaller_group,
    linetype = "dotted",
    color = "#D55E00",
    linewidth = 0.7
  ) +
  geom_line(linewidth = 1.15) +
  geom_point(
    data = n80,
    aes(NPerGroup, Power, color = ScenarioLabel),
    inherit.aes = FALSE,
    size = 2.6
  ) +
  annotate(
    "text",
    x = 442,
    y = 0.845,
    label = "80% target",
    hjust = 1,
    vjust = 0,
    size = 3.1,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = current_smaller_group + 5,
    y = 0.24,
    label = paste0(
      "Current smaller group\nn = ",
      current_smaller_group
    ),
    hjust = 0,
    size = 3.0,
    color = "#B24716"
  ) +
  scale_color_manual(
    values = power_colors,
    name = "Planning scenario"
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(
      10,
      current_smaller_group,
      100,
      200,
      n80$NPerGroup,
      450
    ))),
    limits = c(10, 450),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Power depends on the effect worth detecting",
    subtitle = paste0(
      "Two equal independent groups · two-sided alpha = 0.05 · ",
      "continuous primary endpoint"
    ),
    x = "Independent samples per group",
    y = "Theoretical power"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.2),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

p_power_sensitivity
save_pub(
  p_power_sensitivity,
  "figures/03-power-sensitivity",
  width = 178,
  height = 118
)
