# Alpha 组间检验 + 复现顶刊箱线图
# Run sequentially in a new working directory.
# Required packages: digest, dplyr, emmeans, ggplot2, iNEXT, jsonlite, lmerTest, ragg, readr, scales, svglite, tidyr, vegan.

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

library(digest)
library(dplyr)
library(ggplot2)
library(iNEXT)
library(readr)
library(scales)
library(vegan)

set.seed(20260721)

font_pub <- "sans"

pal_pub <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  sky = "#56B4E9",
  yellow = "#F0E442",
  grey = "#6B7280"
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
      panel.grid.major = ggplot2::element_line(
        colour = "#E6E6E6",
        linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      axis.ticks = ggplot2::element_line(
        colour = "#1A1A1A",
        linewidth = 0.3
      ),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.15)
      ),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(
        colour = "#666666",
        hjust = 0
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2",
        colour = "#B3B3B3",
        linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 89,
  height = 70,
  units = "mm",
  dpi = 600,
  write_svg = TRUE,
  write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"),
    plot,
    width = width,
    height = height,
    units = units,
    device = grDevices::cairo_pdf,
    family = base_family,
    bg = "white"
  )
  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"),
      plot,
      width = width,
      height = height,
      units = units,
      device = svglite::svglite,
      bg = "white"
    )
  }
  ggplot2::ggsave(
    paste0(file_base, ".png"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )
  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"),
      plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw",
      bg = "white"
    )
  }
  invisible(plot)
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character()
  )
  ids <- as.character(x[[1L]])
  stopifnot(!anyDuplicated(ids))
  out <- as.data.frame(x[-1L], check.names = FALSE)
  rownames(out) <- ids
  out
}

input_paths <- c(
  otutab = "data/small/otutab.tsv",
  taxonomy = "data/small/taxonomy.tsv",
  metadata = "data/small/metadata.tsv"
)
expected_sha256 <- c(
  otutab = "76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3",
  taxonomy = "725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901",
  metadata = "df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64"
)
observed_sha256 <- vapply(
  input_paths,
  digest::digest,
  character(1),
  algo = "sha256",
  serialize = FALSE,
  file = TRUE
)
stopifnot(identical(observed_sha256, expected_sha256))

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])

otutab <- data.matrix(otutab)
storage.mode(otutab) <- "numeric"
feature_ids <- rownames(otutab)
sample_ids <- colnames(otutab)

stopifnot(
  setequal(feature_ids, rownames(taxonomy)),
  setequal(sample_ids, rownames(metadata)),
  all(is.finite(otutab)),
  all(otutab >= 0),
  all(abs(otutab - round(otutab)) < .Machine$double.eps^0.5),
  all(colSums(otutab) > 0)
)
taxonomy <- taxonomy[feature_ids, , drop = FALSE]
metadata <- metadata[sample_ids, , drop = FALSE]

dim(otutab)
sum(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

subject_candidates <- c("SubjectID", "Subject", "ParticipantID")
pair_candidates <- c("PairID", "Pair", "MatchedSet")
time_candidates <- c("Time", "Timepoint", "Visit")
block_candidates <- c("Block", "Batch", "SiteID", "Center")

present_subject <- intersect(subject_candidates, colnames(metadata))
present_pair <- intersect(pair_candidates, colnames(metadata))
present_time <- intersect(time_candidates, colnames(metadata))
present_block <- intersect(block_candidates, colnames(metadata))

design_audit <- data.frame(
  Branch = c(
    "Independent groups",
    "Paired two-condition",
    "Repeated measures",
    "Blocked or clustered"
  ),
  RequiredMetadata = c(
    "Group with independent experimental units",
    "SubjectID or PairID plus two conditions",
    "SubjectID plus Time and Condition",
    "Block, SiteID, Center, Family, Cage, or Batch"
  ),
  Available = c(
    "Group" %in% colnames(metadata),
    length(c(present_subject, present_pair)) > 0,
    length(present_subject) > 0 && length(present_time) > 0,
    length(present_block) > 0
  )
)

selected_branch <- "Independent three-group"
stopifnot(
  all(c("Group", "Type", "Saline") %in% colnames(metadata)),
  length(present_subject) == 0L,
  length(present_pair) == 0L,
  length(present_time) == 0L,
  identical(sort(unique(metadata$Group)), c("CW", "IW", "TW")),
  identical(as.integer(table(metadata$Group)[c("CW", "IW", "TW")]),
            c(30L, 30L, 30L))
)

design_audit
selected_branch
with(metadata, table(Group, Saline))

community <- t(otutab)
library_size <- rowSums(community)
observed <- vegan::specnumber(community)
shannon <- vegan::diversity(community, index = "shannon")
inverse_simpson <- vegan::diversity(
  community,
  index = "invsimpson"
)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
group_levels <- unname(group_labels[c("IW", "CW", "TW")])
group_palette <- c(
  "Inland wetland" = pal_pub[["blue"]],
  "Coastal wetland" = pal_pub[["vermillion"]],
  "Tibetan Plateau" = pal_pub[["green"]]
)
group_shapes <- c(
  "Inland wetland" = 21,
  "Coastal wetland" = 22,
  "Tibetan Plateau" = 24
)

alpha_raw <- rbind(
  data.frame(
    SampleID = sample_ids,
    Order = 0L,
    HillD = as.numeric(observed)
  ),
  data.frame(
    SampleID = sample_ids,
    Order = 1L,
    HillD = as.numeric(exp(shannon))
  ),
  data.frame(
    SampleID = sample_ids,
    Order = 2L,
    HillD = as.numeric(inverse_simpson)
  )
)
alpha_raw$Group <- metadata[alpha_raw$SampleID, "Group"]
alpha_raw$GroupLabel <- factor(
  unname(group_labels[alpha_raw$Group]),
  levels = group_levels
)
alpha_raw$Standardization <- "Raw counts"

abundance_list <- setNames(
  lapply(seq_len(nrow(community)), function(i) {
    as.numeric(community[i, community[i, ] > 0])
  }),
  sample_ids
)

size_target <- 10000L
set.seed(20260721)
size_estimate <- iNEXT::estimateD(
  abundance_list,
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "size",
  level = size_target,
  nboot = 0
)
alpha_size <- data.frame(
  SampleID = size_estimate$Assemblage,
  Group = metadata[size_estimate$Assemblage, "Group"],
  GroupLabel = unname(
    group_labels[metadata[size_estimate$Assemblage, "Group"]]
  ),
  ImpliedReads = size_estimate$m,
  Method = size_estimate$Method,
  Order = as.integer(size_estimate$Order.q),
  SampleCoverage = size_estimate$SC,
  HillD = size_estimate$qD,
  Standardization = "10,000 reads",
  stringsAsFactors = FALSE
)
alpha_size$GroupLabel <- factor(
  alpha_size$GroupLabel,
  levels = group_levels
)

coverage_target <- 0.90
set.seed(20260721)
coverage_estimate <- iNEXT::estimateD(
  abundance_list,
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = coverage_target,
  nboot = 0
)
alpha_coverage <- data.frame(
  SampleID = coverage_estimate$Assemblage,
  Group = metadata[coverage_estimate$Assemblage, "Group"],
  GroupLabel = unname(
    group_labels[metadata[coverage_estimate$Assemblage, "Group"]]
  ),
  ObservedReads = as.numeric(
    library_size[coverage_estimate$Assemblage]
  ),
  ImpliedReads = coverage_estimate$m,
  EffortRatio = coverage_estimate$m / as.numeric(
    library_size[coverage_estimate$Assemblage]
  ),
  Method = coverage_estimate$Method,
  Order = as.integer(coverage_estimate$Order.q),
  TargetCoverage = coverage_estimate$SC,
  HillD = coverage_estimate$qD,
  Standardization = "90% coverage",
  stringsAsFactors = FALSE
)
alpha_coverage$GroupLabel <- factor(
  alpha_coverage$GroupLabel,
  levels = group_levels
)

coverage_q0 <- alpha_coverage[
  alpha_coverage$Order == 0L,
  ,
  drop = FALSE
]
stopifnot(
  nrow(alpha_raw) == 270L,
  nrow(alpha_size) == 270L,
  nrow(alpha_coverage) == 270L,
  length(unique(alpha_coverage$SampleID)) == 90L,
  identical(sort(unique(alpha_coverage$Order)), 0:2),
  sum(coverage_q0$Method == "Rarefaction") == 56L,
  sum(coverage_q0$Method == "Extrapolation") == 34L,
  max(coverage_q0$EffortRatio) <= 1.5
)

c(
  samples = length(unique(alpha_coverage$SampleID)),
  rarefaction_samples = sum(coverage_q0$Method == "Rarefaction"),
  extrapolation_samples = sum(coverage_q0$Method == "Extrapolation"),
  maximum_effort_ratio = max(coverage_q0$EffortRatio)
)

sample_skewness <- function(x) {
  n <- length(x)
  s <- sd(x)
  if (n < 3L || !is.finite(s) || s == 0) {
    return(NA_real_)
  }
  n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / s)^3)
}

group_descriptive_summary <- alpha_coverage |>
  group_by(Order, Group, GroupLabel) |>
  summarise(
    Samples = n(),
    Mean = mean(HillD),
    SD = sd(HillD),
    Median = median(HillD),
    Q25 = quantile(HillD, 0.25),
    Q75 = quantile(HillD, 0.75),
    Minimum = min(HillD),
    Maximum = max(HillD),
    .groups = "drop"
  )

distribution_diagnostics <- alpha_coverage |>
  group_by(Order, Group, GroupLabel) |>
  summarise(
    Samples = n(),
    UniqueValues = n_distinct(HillD),
    Skewness = sample_skewness(HillD),
    IQR = IQR(HillD),
    LowerOutliers = sum(
      HillD < quantile(HillD, 0.25) - 1.5 * IQR(HillD)
    ),
    UpperOutliers = sum(
      HillD > quantile(HillD, 0.75) + 1.5 * IQR(HillD)
    ),
    .groups = "drop"
  )

group_descriptive_summary
distribution_diagnostics

omnibus_rows <- lapply(0:2, function(order_q) {
  x <- alpha_coverage[
    alpha_coverage$Order == order_q,
    ,
    drop = FALSE
  ]
  fit <- kruskal.test(HillD ~ GroupLabel, data = x)
  h <- unname(fit$statistic)
  k <- nlevels(droplevels(x$GroupLabel))
  n_total <- nrow(x)
  data.frame(
    Order = order_q,
    Metric = paste0("Hill q=", order_q),
    H = h,
    DF = unname(fit$parameter),
    PValue = fit$p.value,
    EpsilonSquared = max(0, (h - k + 1) / (n_total - k))
  )
})
omnibus_tests <- bind_rows(omnibus_rows) |>
  mutate(
    PAdjustedHolm = p.adjust(PValue, method = "holm"),
    RejectHolm05 = PAdjustedHolm < 0.05
  )

omnibus_tests

cliffs_delta <- function(x, y) {
  differences <- outer(x, y, "-")
  (sum(differences > 0) - sum(differences < 0)) /
    length(differences)
}

bootstrap_delta <- function(x, y, replicates, seed) {
  set.seed(seed)
  values <- vapply(
    seq_len(replicates),
    function(i) {
      xb <- sample(x, length(x), replace = TRUE)
      yb <- sample(y, length(y), replace = TRUE)
      cliffs_delta(xb, yb)
    },
    numeric(1)
  )
  c(
    lower = unname(quantile(values, 0.025, names = FALSE)),
    upper = unname(quantile(values, 0.975, names = FALSE)),
    valid = sum(is.finite(values))
  )
}

bootstrap_replicates <- 5000L
pair_definitions <- combn(group_levels, 2L, simplify = FALSE)
pairwise_rows <- list()
row_id <- 0L

for (order_q in 0:2) {
  order_data <- alpha_coverage[
    alpha_coverage$Order == order_q,
    ,
    drop = FALSE
  ]
  for (pair in pair_definitions) {
    row_id <- row_id + 1L
    x <- order_data$HillD[
      order_data$GroupLabel == pair[[1L]]
    ]
    y <- order_data$HillD[
      order_data$GroupLabel == pair[[2L]]
    ]
    wilcox_fit <- suppressWarnings(
      wilcox.test(x, y, exact = FALSE, correct = TRUE)
    )
    delta <- cliffs_delta(x, y)
    delta_ci <- bootstrap_delta(
      x,
      y,
      replicates = bootstrap_replicates,
      seed = 20260721L + 100L * order_q + row_id
    )

    pairwise_rows[[row_id]] <- data.frame(
      Order = order_q,
      GroupA = pair[[1L]],
      GroupB = pair[[2L]],
      Contrast = paste(pair[[1L]], "minus", pair[[2L]]),
      NGroupA = length(x),
      NGroupB = length(y),
      MedianGroupA = median(x),
      MedianGroupB = median(y),
      WilcoxonW = unname(wilcox_fit$statistic),
      PValue = wilcox_fit$p.value,
      CliffsDelta = delta,
      ProbabilitySuperiority = (delta + 1) / 2,
      DeltaLower95 = delta_ci[["lower"]],
      DeltaUpper95 = delta_ci[["upper"]],
      BootstrapReplicates = bootstrap_replicates,
      BootstrapValid = as.integer(delta_ci[["valid"]])
    )
  }
}

pairwise_tests <- bind_rows(pairwise_rows) |>
  mutate(
    PAdjustedHolm = p.adjust(PValue, method = "holm"),
    RejectHolm05 = PAdjustedHolm < 0.05
  )

stopifnot(
  nrow(pairwise_tests) == 9L,
  all(pairwise_tests$BootstrapValid == bootstrap_replicates),
  all(abs(pairwise_tests$CliffsDelta) <= 1),
  all(pairwise_tests$DeltaLower95 <= pairwise_tests$CliffsDelta),
  all(pairwise_tests$CliffsDelta <= pairwise_tests$DeltaUpper95)
)

pairwise_tests

welch_sensitivity <- bind_rows(lapply(0:2, function(order_q) {
  x <- alpha_coverage[
    alpha_coverage$Order == order_q,
    ,
    drop = FALSE
  ]
  fit <- oneway.test(
    HillD ~ GroupLabel,
    data = x,
    var.equal = FALSE
  )
  data.frame(
    Order = order_q,
    Statistic = unname(fit$statistic),
    DF1 = unname(fit$parameter[[1L]]),
    DF2 = unname(fit$parameter[[2L]]),
    PValue = fit$p.value
  )
})) |>
  mutate(
    PAdjustedHolm = p.adjust(PValue, method = "holm"),
    RejectHolm05 = PAdjustedHolm < 0.05
  )

standardized_hill <- bind_rows(
  alpha_raw |>
    select(SampleID, Group, GroupLabel, Order, HillD, Standardization),
  alpha_size |>
    select(SampleID, Group, GroupLabel, Order, HillD, Standardization),
  alpha_coverage |>
    select(SampleID, Group, GroupLabel, Order, HillD, Standardization)
) |>
  mutate(
    Standardization = factor(
      Standardization,
      levels = c("Raw counts", "10,000 reads", "90% coverage")
    )
  )

standardization_sensitivity <- standardized_hill |>
  group_by(Standardization, Order) |>
  group_modify(~ {
    fit <- kruskal.test(HillD ~ GroupLabel, data = .x)
    h <- unname(fit$statistic)
    k <- nlevels(droplevels(.x$GroupLabel))
    n_total <- nrow(.x)
    data.frame(
      H = h,
      DF = unname(fit$parameter),
      PValue = fit$p.value,
      EpsilonSquared = max(0, (h - k + 1) / (n_total - k))
    )
  }) |>
  ungroup() |>
  arrange(Standardization, Order) |>
  mutate(
    PAdjustedHolm9 = p.adjust(PValue, method = "holm"),
    RejectHolm9At05 = PAdjustedHolm9 < 0.05
  )

welch_sensitivity
standardization_sensitivity

design_plot_data <- data.frame(
  Branch = factor(
    c(
      "Independent groups",
      "Paired two-condition",
      "Repeated measures",
      "Blocked or clustered"
    ),
    levels = rev(c(
      "Independent groups",
      "Paired two-condition",
      "Repeated measures",
      "Blocked or clustered"
    ))
  ),
  Detail = c(
    "Group: IW / CW / TW; n = 30 each",
    "Needs SubjectID or PairID",
    "Needs SubjectID and Time",
    "Needs a declared clustering column"
  ),
  Status = c(
    "Selected",
    "Metadata absent",
    "Metadata absent",
    "Metadata absent"
  )
)

design_plot <- ggplot(
  design_plot_data,
  aes(y = Branch)
) +
  geom_segment(
    aes(x = 0.08, xend = 0.92, yend = Branch),
    colour = "#D0D0D0",
    linewidth = 1.2
  ) +
  geom_point(
    aes(x = 0.08, fill = Status),
    shape = 21,
    size = 4.2,
    stroke = 0.45,
    colour = "#1A1A1A"
  ) +
  geom_text(
    aes(x = 0.15, label = Detail),
    hjust = 0,
    size = 3.2,
    colour = "#1A1A1A"
  ) +
  geom_text(
    aes(x = 0.92, label = Status),
    hjust = 1,
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Selected" = pal_pub[["green"]],
      "Metadata absent" = "#B8B8B8"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Choose the test from the sampling design",
    subtitle = "Selected branch: independent three-group comparison",
    x = NULL,
    y = NULL,
    caption = paste(
      "No subject, pair, time, or declared block identifier is present",
      "in the demonstration metadata."
    )
  ) +
  theme_pub(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(face = "bold", size = 8.5),
    plot.margin = margin(8, 12, 8, 8)
  )

design_plot
save_pub(
  design_plot,
  "figures/20-design-selection-map",
  width = 183,
  height = 84
)

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

order_descriptions <- c(
  "0" = "q = 0 · rare features",
  "1" = "q = 1 · balanced",
  "2" = "q = 2 · dominant features"
)
panel_lookup <- setNames(
  paste0(
    order_descriptions[as.character(omnibus_tests$Order)],
    "\nKruskal-Wallis Holm P ",
    format_p(omnibus_tests$PAdjustedHolm)
  ),
  omnibus_tests$Order
)

distribution_plot_data <- alpha_coverage |>
  mutate(
    Panel = factor(
      unname(panel_lookup[as.character(Order)]),
      levels = unname(panel_lookup[as.character(0:2)])
    )
  )

distribution_plot <- ggplot(
  distribution_plot_data,
  aes(x = GroupLabel, y = HillD, fill = GroupLabel)
) +
  geom_violin(
    width = 0.88,
    alpha = 0.16,
    colour = NA,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    alpha = 0.60,
    colour = "#1A1A1A",
    linewidth = 0.4
  ) +
  geom_point(
    aes(shape = GroupLabel),
    position = position_jitter(
      width = 0.11,
      height = 0,
      seed = 20260721
    ),
    colour = "#1A1A1A",
    stroke = 0.28,
    size = 1.45,
    alpha = 0.68
  ) +
  facet_wrap(~Panel, scales = "free_y", nrow = 1L) +
  scale_fill_manual(values = group_palette, guide = "none") +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_log10(labels = label_comma()) +
  labs(
    title = "Alpha diversity by wetland group",
    subtitle = paste(
      "Hill diversity standardized to 90% sample coverage;",
      "points are independent samples"
    ),
    x = NULL,
    y = "Effective features (log scale)",
    caption = paste(
      "Boxes show medians and IQRs. Kruskal-Wallis tests rank",
      "distributions, not medians alone."
    )
  ) +
  theme_pub(base_size = 8.6) +
  theme(
    axis.text.x = element_text(angle = 24, hjust = 1, size = 7.5),
    strip.text = element_text(size = 7.5),
    legend.position = "none"
  )

distribution_plot
save_pub(
  distribution_plot,
  "figures/20-alpha-group-distributions",
  width = 183,
  height = 108
)

sensitivity_plot_data <- standardization_sensitivity |>
  mutate(
    OrderLabel = factor(
      paste0("q = ", Order),
      levels = paste0("q = ", 0:2)
    ),
    MinusLog10P = -log10(PValue),
    PLabel = paste0("P=", format_p(PValue))
  )

sensitivity_plot <- ggplot(
  sensitivity_plot_data,
  aes(
    x = Standardization,
    y = MinusLog10P,
    group = OrderLabel,
    colour = OrderLabel
  )
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    colour = "#666666",
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.7, alpha = 0.75) +
  geom_point(
    aes(size = EpsilonSquared),
    shape = 21,
    fill = "white",
    stroke = 0.75
  ) +
  geom_text(
    aes(label = PLabel),
    nudge_y = 0.10,
    size = 2.55,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "q = 0" = pal_pub[["blue"]],
      "q = 1" = pal_pub[["orange"]],
      "q = 2" = pal_pub[["purple"]]
    )
  ) +
  scale_size_continuous(
    range = c(2.2, 6.2),
    limits = c(0, max(sensitivity_plot_data$EpsilonSquared)),
    name = expression(epsilon^2)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.22)),
    breaks = c(0, -log10(0.05), 2),
    labels = c("1", "0.05", "0.01")
  ) +
  labs(
    title = paste(
      "The nominal raw-richness signal is not robust",
      "to standardization"
    ),
    subtitle = paste(
      "Kruskal-Wallis diagnostics across prespecified",
      "measurement scales"
    ),
    x = NULL,
    y = "Unadjusted P value (reverse log scale)",
    colour = "Hill order",
    caption = paste(
      "Point size is epsilon-squared. All nine sensitivity P values",
      "are non-significant after global Holm adjustment."
    )
  ) +
  theme_pub(base_size = 9) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(face = "bold")
  )

sensitivity_plot
save_pub(
  sensitivity_plot,
  "figures/20-standardization-sensitivity",
  width = 183,
  height = 92
)

forest_plot_data <- pairwise_tests |>
  mutate(
    OrderPanel = factor(
      paste0("q = ", Order),
      levels = paste0("q = ", 0:2)
    ),
    ContrastLabel = factor(
      Contrast,
      levels = rev(unique(Contrast))
    ),
    PLabel = paste0("Holm P=", format_p(PAdjustedHolm))
  )

forest_plot <- ggplot(
  forest_plot_data,
  aes(x = CliffsDelta, y = ContrastLabel)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "#666666",
    linewidth = 0.45
  ) +
  geom_errorbar(
    aes(xmin = DeltaLower95, xmax = DeltaUpper95),
    orientation = "y",
    width = 0.16,
    linewidth = 0.6,
    colour = "#3A3A3A"
  ) +
  geom_point(
    shape = 21,
    size = 2.6,
    stroke = 0.45,
    fill = pal_pub[["blue"]],
    colour = "#1A1A1A"
  ) +
  geom_text(
    aes(label = PLabel),
    nudge_y = 0.23,
    size = 2.45,
    colour = "#333333"
  ) +
  facet_wrap(~OrderPanel, nrow = 1L) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  labs(
    title = "Pairwise dominance effects remain uncertain",
    subtitle = paste(
      "Cliff's delta with stratified bootstrap",
      "95% percentile intervals"
    ),
    x = "Cliff's delta (Group A minus Group B)",
    y = NULL,
    caption = paste(
      paste(
        "Positive values mean that a random Group A sample tends",
        "to have higher diversity than Group B."
      ),
      "Holm adjustment covers all nine planned contrasts.",
      sep = "\n"
    )
  ) +
  theme_pub(base_size = 8.5) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 7.3)
  )

forest_plot
save_pub(
  forest_plot,
  "figures/20-pairwise-effect-forest",
  width = 183,
  height = 106
)

audit_summary_path <- paste0(
  "results/20-alpha-group-tests/",
  "alpha-group-tests-summary.json"
)
if (file.exists(audit_summary_path)) {
  audit_summary <- jsonlite::read_json(
    audit_summary_path,
    simplifyVector = TRUE
  )
  stopifnot(
    audit_summary$checks_total == 188L,
    audit_summary$checks_passed == 188L,
    audit_summary$checks_failed == 0L,
    audit_summary$analysis$omnibus_tests == 3L,
    audit_summary$analysis$pairwise_tests == 9L,
    audit_summary$analysis$bootstrap_replicates == 5000L,
    audit_summary$graphics$primary_figures == 4L,
    audit_summary$graphics$format_files == 16L
  )
  unlist(audit_summary[c(
    "checks_total",
    "checks_passed",
    "checks_failed"
  )])
}
