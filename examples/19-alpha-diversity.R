# Alpha 指数全家 + 稀释曲线 + 抽平之争 + Hill numbers
# Run sequentially in a new working directory.
# Required packages: dplyr, ggplot2, iNEXT, jsonlite, phyloseq, ragg, scales, svglite, systemfonts, vegan.

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

set.seed(20260721)

required_packages <- c(
  "digest", "dplyr", "ggplot2", "iNEXT", "jsonlite", "phyloseq",
  "ragg", "readr", "scales", "svglite", "systemfonts", "tibble",
  "tidyr", "vegan"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
stopifnot(length(missing_packages) == 0L)

font_pub <- "sans"
font_info <- systemfonts::font_info(font_pub)
stopifnot(
  nrow(font_info) >= 1L,
  nzchar(font_info$family[[1L]]),
  file.exists(font_info$path[[1L]]),
  isTRUE(font_info$scalable[[1L]])
)

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
        colour = "#E6E6E6", linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(colour = "#1A1A1A"),
      axis.ticks = ggplot2::element_line(
        colour = "#1A1A1A", linewidth = 0.3
      ),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.15)
      ),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D"),
      plot.caption = ggplot2::element_text(
        colour = "#666666", hjust = 0
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 183,
  height = 100,
  units = "mm",
  dpi = 600,
  write_svg = TRUE,
  write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot,
    width = width, height = height, units = units,
    device = grDevices::cairo_pdf,
    family = base_family,
    bg = "white"
  )

  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"), plot,
      width = width, height = height, units = units,
      device = svglite::svglite,
      bg = "white"
    )
  }

  ggplot2::ggsave(
    paste0(file_base, ".png"), plot,
    width = width, height = height, units = units,
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )

  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"), plot,
      width = width, height = height, units = units,
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

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(otutab, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(taxonomy)),
  phyloseq::sample_data(metadata)
)
ps

community <- t(otutab)
library_size <- rowSums(community)
observed <- vegan::specnumber(community)
frequency_one <- rowSums(community == 1)
frequency_two <- rowSums(community == 2)
goods_coverage <- 1 - frequency_one / library_size
estimate_r <- vegan::estimateR(community)
shannon <- vegan::diversity(community, index = "shannon")
simpson <- vegan::diversity(community, index = "simpson")
inverse_simpson <- vegan::diversity(community, index = "invsimpson")
pielou <- shannon / log(observed)

group_labels <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan Plateau"
)
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

alpha_raw <- data.frame(
  SampleID = sample_ids,
  Group = metadata$Group,
  GroupLabel = unname(group_labels[metadata$Group]),
  LibrarySize = as.numeric(library_size),
  Singletons = as.integer(frequency_one),
  Doubletons = as.integer(frequency_two),
  GoodsCoverage = as.numeric(goods_coverage),
  Observed = as.numeric(observed),
  Chao1 = as.numeric(estimate_r["S.chao1", ]),
  Chao1SE = as.numeric(estimate_r["se.chao1", ]),
  Shannon = as.numeric(shannon),
  Simpson = as.numeric(simpson),
  Pielou = as.numeric(pielou),
  Hill0 = as.numeric(observed),
  Hill1 = as.numeric(exp(shannon)),
  Hill2 = as.numeric(inverse_simpson),
  stringsAsFactors = FALSE
)
alpha_raw$GroupLabel <- factor(
  alpha_raw$GroupLabel,
  levels = names(group_palette)
)

phyloseq_alpha <- suppressWarnings(
  phyloseq::estimate_richness(
    ps,
    measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
  )
)
phyloseq_alpha <- phyloseq_alpha[sample_ids, , drop = FALSE]
stopifnot(
  max(abs(alpha_raw$Observed - phyloseq_alpha$Observed)) < 1e-10,
  max(abs(alpha_raw$Chao1 - phyloseq_alpha$Chao1)) < 1e-8,
  max(abs(alpha_raw$Shannon - phyloseq_alpha$Shannon)) < 1e-10,
  max(abs(alpha_raw$Simpson - phyloseq_alpha$Simpson)) < 1e-10,
  max(abs(alpha_raw$Hill2 - phyloseq_alpha$InvSimpson)) < 1e-10,
  max(abs(alpha_raw$Hill1 - exp(alpha_raw$Shannon))) < 1e-10,
  max(abs(alpha_raw$Hill2 - 1 / (1 - alpha_raw$Simpson))) < 1e-8
)

depth_observed_spearman <- suppressWarnings(
  cor(alpha_raw$LibrarySize, alpha_raw$Observed, method = "spearman")
)
alpha_raw[1:6, ]
c(
  minimum_depth = min(alpha_raw$LibrarySize),
  median_depth = median(alpha_raw$LibrarySize),
  maximum_depth = max(alpha_raw$LibrarySize),
  minimum_coverage = min(alpha_raw$GoodsCoverage),
  median_coverage = median(alpha_raw$GoodsCoverage),
  maximum_coverage = max(alpha_raw$GoodsCoverage),
  depth_observed_spearman = depth_observed_spearman
)

abundance_list <- setNames(
  lapply(seq_len(nrow(community)), function(i) {
    as.numeric(community[i, community[i, ] > 0])
  }),
  sample_ids
)

data_info <- iNEXT::DataInfo(
  abundance_list,
  datatype = "abundance"
)
data_info <- data_info[
  match(sample_ids, data_info$Assemblage),
  ,
  drop = FALSE
]
alpha_raw$iNEXTCoverage <- data_info$SC
stopifnot(
  max(abs(alpha_raw$GoodsCoverage - alpha_raw$iNEXTCoverage)) <= 0.0001
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
  TargetReads = as.numeric(size_estimate$m),
  Method = size_estimate$Method,
  Order = as.integer(size_estimate$Order.q),
  SampleCoverage = size_estimate$SC,
  HillD = size_estimate$qD,
  stringsAsFactors = FALSE
)
alpha_size$GroupLabel <- factor(
  alpha_size$GroupLabel,
  levels = names(group_palette)
)
stopifnot(
  length(unique(alpha_size$SampleID)) == 90L,
  identical(sort(unique(alpha_size$Order)), 0:2),
  identical(unique(alpha_size$Method), "Rarefaction")
)

depth_grid <- seq(1000L, size_target, by = 1000L)
rarefaction_matrix <- vapply(
  depth_grid,
  function(depth) as.numeric(
    vegan::rarefy(community, sample = depth)
  ),
  numeric(nrow(community))
)
rownames(rarefaction_matrix) <- sample_ids
colnames(rarefaction_matrix) <- as.character(depth_grid)
rarefaction_curve_data <- data.frame(
  SampleID = rep(sample_ids, times = length(depth_grid)),
  TargetReads = rep(depth_grid, each = length(sample_ids)),
  ExpectedObserved = as.numeric(rarefaction_matrix),
  stringsAsFactors = FALSE
)
rarefaction_curve_data$Group <- metadata[
  rarefaction_curve_data$SampleID,
  "Group"
]
rarefaction_curve_data$GroupLabel <- factor(
  unname(group_labels[rarefaction_curve_data$Group]),
  levels = names(group_palette)
)

size_cost <- data.frame(
  SamplesRetained = length(unique(alpha_size$SampleID)),
  ReadsAvailable = sum(community),
  ReadsRepresented = size_target * length(unique(alpha_size$SampleID)),
  FractionNotUsed = 1 -
    size_target * length(unique(alpha_size$SampleID)) / sum(community)
)
size_cost

coverage_levels <- c(0.88, 0.90, 0.92)
coverage_sensitivity <- do.call(
  rbind,
  lapply(coverage_levels, function(level) {
    set.seed(20260721)
    estimate <- iNEXT::estimateD(
      abundance_list,
      q = 0,
      datatype = "abundance",
      base = "coverage",
      level = level,
      nboot = 0
    )
    observed_reads <- as.numeric(library_size[estimate$Assemblage])
    data.frame(
      TargetCoverage = level,
      Samples = nrow(estimate),
      RarefactionSamples = sum(estimate$Method == "Rarefaction"),
      ObservedSamples = sum(estimate$Method == "Observed"),
      ExtrapolationSamples = sum(estimate$Method == "Extrapolation"),
      MedianImpliedReads = median(estimate$m),
      MaximumEffortRatio = max(estimate$m / observed_reads)
    )
  })
)
coverage_sensitivity

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
  ImpliedReads = coverage_estimate$m,
  ObservedReads = as.numeric(
    library_size[coverage_estimate$Assemblage]
  ),
  EffortRatio = coverage_estimate$m / as.numeric(
    library_size[coverage_estimate$Assemblage]
  ),
  Method = coverage_estimate$Method,
  Order = as.integer(coverage_estimate$Order.q),
  TargetCoverage = coverage_estimate$SC,
  HillD = coverage_estimate$qD,
  stringsAsFactors = FALSE
)
alpha_coverage$GroupLabel <- factor(
  alpha_coverage$GroupLabel,
  levels = names(group_palette)
)
coverage_q0 <- alpha_coverage[
  alpha_coverage$Order == 0L,
  ,
  drop = FALSE
]
stopifnot(
  sum(coverage_q0$Method == "Rarefaction") == 56L,
  sum(coverage_q0$Method == "Extrapolation") == 34L,
  max(coverage_q0$EffortRatio) <= 1.5,
  max(abs(alpha_coverage$TargetCoverage - coverage_target)) < 1e-8
)
c(
  rarefaction_samples = sum(coverage_q0$Method == "Rarefaction"),
  extrapolation_samples = sum(coverage_q0$Method == "Extrapolation"),
  median_implied_reads = median(coverage_q0$ImpliedReads),
  maximum_effort_ratio = max(coverage_q0$EffortRatio)
)

raw_hill <- rbind(
  data.frame(
    SampleID = alpha_raw$SampleID,
    Order = 0L,
    HillD = alpha_raw$Hill0
  ),
  data.frame(
    SampleID = alpha_raw$SampleID,
    Order = 1L,
    HillD = alpha_raw$Hill1
  ),
  data.frame(
    SampleID = alpha_raw$SampleID,
    Order = 2L,
    HillD = alpha_raw$Hill2
  )
)
raw_hill$Group <- metadata[raw_hill$SampleID, "Group"]
raw_hill$GroupLabel <- unname(group_labels[raw_hill$Group])

safe_spearman <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman"))
}
rank_rows <- list()
for (order_q in 0:2) {
  raw_q <- raw_hill[
    raw_hill$Order == order_q,
    c("SampleID", "HillD")
  ]
  names(raw_q)[[2L]] <- "RawHillD"
  size_q <- alpha_size[
    alpha_size$Order == order_q,
    c("SampleID", "HillD")
  ]
  names(size_q)[[2L]] <- "SizeHillD"
  coverage_q <- alpha_coverage[
    alpha_coverage$Order == order_q,
    c("SampleID", "HillD")
  ]
  names(coverage_q)[[2L]] <- "CoverageHillD"
  joined <- Reduce(
    function(x, y) merge(x, y, by = "SampleID", sort = FALSE),
    list(raw_q, size_q, coverage_q)
  )
  rank_rows[[length(rank_rows) + 1L]] <- data.frame(
    Order = order_q,
    Comparison = c(
      "Raw vs 10,000 reads",
      "Raw vs 90% coverage",
      "10,000 reads vs 90% coverage"
    ),
    SpearmanRho = c(
      safe_spearman(joined$RawHillD, joined$SizeHillD),
      safe_spearman(joined$RawHillD, joined$CoverageHillD),
      safe_spearman(joined$SizeHillD, joined$CoverageHillD)
    )
  )
}
rank_stability <- do.call(rbind, rank_rows)
rank_stability

standardized_hill <- rbind(
  transform(
    raw_hill[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "Raw counts"
  ),
  transform(
    alpha_size[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "10,000 reads"
  ),
  transform(
    alpha_coverage[c("SampleID", "Group", "GroupLabel", "Order", "HillD")],
    Standardization = "90% coverage"
  )
)
standardized_hill$Standardization <- factor(
  standardized_hill$Standardization,
  levels = c("Raw counts", "10,000 reads", "90% coverage")
)
group_descriptive_summary <- standardized_hill |>
  dplyr::group_by(Standardization, Group, GroupLabel, Order) |>
  dplyr::summarise(
    Samples = dplyr::n(),
    Median = median(HillD),
    Q25 = quantile(HillD, 0.25),
    Q75 = quantile(HillD, 0.75),
    .groups = "drop"
  )
group_descriptive_summary

depth_plot_data <- rbind(
  data.frame(
    SampleID = alpha_raw$SampleID,
    GroupLabel = alpha_raw$GroupLabel,
    LibrarySize = alpha_raw$LibrarySize,
    Metric = "Observed features (count)",
    Value = alpha_raw$Observed
  ),
  data.frame(
    SampleID = alpha_raw$SampleID,
    GroupLabel = alpha_raw$GroupLabel,
    LibrarySize = alpha_raw$LibrarySize,
    Metric = "Estimated sample coverage (%)",
    Value = 100 * alpha_raw$iNEXTCoverage
  )
)
depth_plot_data$Metric <- factor(
  depth_plot_data$Metric,
  levels = c(
    "Observed features (count)",
    "Estimated sample coverage (%)"
  )
)

depth_plot <- ggplot2::ggplot(
  depth_plot_data,
  ggplot2::aes(x = LibrarySize, y = Value)
) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "#4D4D4D",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = GroupLabel, shape = GroupLabel),
    colour = "#1A1A1A",
    stroke = 0.3,
    size = 2,
    alpha = 0.86
  ) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", nrow = 1L) +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(
    labels = scales::label_comma(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Sequencing depth and completeness are separate diagnostics",
    subtitle = paste0(
      "Raw observed richness remains depth-associated (Spearman rho = ",
      sprintf("%.3f", depth_observed_spearman), ")"
    ),
    x = "Library size (reads)",
    y = NULL,
    fill = "Wetland group",
    shape = "Wetland group",
    caption = paste(
      "Coverage describes the processed feature table;",
      "it is not absolute cell or species coverage."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

save_pub(
  depth_plot,
  "figures/19-depth-completeness-audit",
  width = 183,
  height = 92
)

curve_summary <- rarefaction_curve_data |>
  dplyr::group_by(GroupLabel, TargetReads) |>
  dplyr::summarise(
    Median = median(ExpectedObserved),
    Q25 = quantile(ExpectedObserved, 0.25),
    Q75 = quantile(ExpectedObserved, 0.75),
    .groups = "drop"
  )

rarefaction_plot <- ggplot2::ggplot(
  rarefaction_curve_data,
  ggplot2::aes(x = TargetReads, y = ExpectedObserved)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = SampleID, colour = GroupLabel),
    linewidth = 0.25,
    alpha = 0.18
  ) +
  ggplot2::geom_ribbon(
    data = curve_summary,
    ggplot2::aes(
      x = TargetReads,
      ymin = Q25,
      ymax = Q75,
      fill = GroupLabel
    ),
    inherit.aes = FALSE,
    alpha = 0.16,
    colour = NA
  ) +
  ggplot2::geom_line(
    data = curve_summary,
    ggplot2::aes(y = Median, colour = GroupLabel),
    linewidth = 0.9
  ) +
  ggplot2::geom_vline(
    xintercept = size_target,
    linetype = 2,
    colour = "#4D4D4D",
    linewidth = 0.45
  ) +
  ggplot2::facet_wrap(~GroupLabel, nrow = 1L) +
  ggplot2::scale_colour_manual(values = group_palette, guide = "none") +
  ggplot2::scale_fill_manual(values = group_palette, guide = "none") +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Sample-size rarefaction exposes unsaturated richness",
    subtitle = paste(
      "Thin lines are samples; thick lines and ribbons are",
      "descriptive medians and IQRs"
    ),
    x = "Standardized reads",
    y = "Expected observed features",
    caption = paste(
      "The dashed line marks 10,000 reads;",
      "no pooled curve is used as a replicate."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "none")

save_pub(
  rarefaction_plot,
  "figures/19-rarefaction-curves",
  width = 183,
  height = 100
)

raw_for_plot <- raw_hill[c("SampleID", "Order", "HillD")]
names(raw_for_plot)[[3L]] <- "RawHillD"
size_for_plot <- alpha_size[
  c("SampleID", "GroupLabel", "Order", "HillD")
]
size_for_plot$Standardization <- "10,000 reads"
coverage_for_plot <- alpha_coverage[
  c("SampleID", "GroupLabel", "Order", "HillD")
]
coverage_for_plot$Standardization <- "90% coverage"
comparison_plot_data <- rbind(size_for_plot, coverage_for_plot)
comparison_plot_data <- merge(
  comparison_plot_data,
  raw_for_plot,
  by = c("SampleID", "Order"),
  sort = FALSE
)
comparison_plot_data$OrderLabel <- factor(
  paste0("q = ", comparison_plot_data$Order),
  levels = paste0("q = ", 0:2)
)
comparison_plot_data$Standardization <- factor(
  comparison_plot_data$Standardization,
  levels = c("10,000 reads", "90% coverage")
)

standardization_plot <- ggplot2::ggplot(
  comparison_plot_data,
  ggplot2::aes(x = RawHillD, y = HillD)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "#7A7A7A",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = GroupLabel, shape = Standardization),
    colour = "#1A1A1A",
    stroke = 0.28,
    size = 1.65,
    alpha = 0.78
  ) +
  ggplot2::facet_wrap(~OrderLabel, scales = "free", nrow = 1L) +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_shape_manual(
    values = c("10,000 reads" = 21, "90% coverage" = 24)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(shape = 21, size = 2.4, alpha = 1)
    ),
    shape = ggplot2::guide_legend(
      override.aes = list(fill = "white", size = 2.4, alpha = 1)
    )
  ) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = paste(
      "Standardization changes rare-feature-sensitive",
      "diversity most"
    ),
    subtitle = paste(
      "Each point compares one standardized estimate with",
      "its raw plug-in estimate"
    ),
    x = "Raw effective features",
    y = "Standardized effective features",
    fill = "Wetland group",
    shape = "Standardization"
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

save_pub(
  standardization_plot,
  "figures/19-standardization-comparison",
  width = 183,
  height = 100
)

coverage_profile <- alpha_coverage
coverage_profile$OrderLabel <- coverage_profile$Order
coverage_profile_summary <- coverage_profile |>
  dplyr::group_by(GroupLabel, OrderLabel) |>
  dplyr::summarise(
    Median = median(HillD),
    Q25 = quantile(HillD, 0.25),
    Q75 = quantile(HillD, 0.75),
    .groups = "drop"
  )

hill_plot <- ggplot2::ggplot(
  coverage_profile,
  ggplot2::aes(x = OrderLabel, y = HillD)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = SampleID, colour = GroupLabel),
    linewidth = 0.25,
    alpha = 0.16
  ) +
  ggplot2::geom_ribbon(
    data = coverage_profile_summary,
    ggplot2::aes(
      x = OrderLabel,
      ymin = Q25,
      ymax = Q75,
      fill = GroupLabel,
      group = GroupLabel
    ),
    inherit.aes = FALSE,
    alpha = 0.14,
    colour = NA
  ) +
  ggplot2::geom_line(
    data = coverage_profile_summary,
    ggplot2::aes(
      y = Median,
      group = GroupLabel,
      colour = GroupLabel
    ),
    linewidth = 0.95
  ) +
  ggplot2::geom_point(
    data = coverage_profile_summary,
    ggplot2::aes(y = Median, fill = GroupLabel),
    shape = 21,
    colour = "#1A1A1A",
    stroke = 0.3,
    size = 2.2
  ) +
  ggplot2::scale_colour_manual(values = group_palette, guide = "none") +
  ggplot2::scale_fill_manual(values = group_palette) +
  ggplot2::scale_x_continuous(
    breaks = 0:2,
    labels = c(
      "q = 0\nrare taxa",
      "q = 1\nbalanced",
      "q = 2\ndominant taxa"
    )
  ) +
  ggplot2::scale_y_log10(labels = scales::label_comma()) +
  ggplot2::labs(
    title = paste(
      "Hill profiles put all three orders in",
      "effective-feature units"
    ),
    subtitle = "Diversity standardized to 90% sample coverage",
    x = "Hill order",
    y = "Effective features (log scale)",
    fill = "Wetland group",
    caption = paste(
      "Lines and IQRs are descriptive; group hypothesis tests",
      "are not performed here."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(legend.position = "top")

save_pub(
  hill_plot,
  "figures/19-hill-profile",
  width = 183,
  height = 92
)

audit_summary_path <- paste0(
  "results/19-alpha-diversity/",
  "alpha-diversity-summary.json"
)
if (file.exists(audit_summary_path)) {
  audit_summary <- jsonlite::read_json(
    audit_summary_path,
    simplifyVector = TRUE
  )
  stopifnot(
    audit_summary$checks_total == 193L,
    audit_summary$checks_passed == 193L,
    audit_summary$checks_failed == 0L,
    audit_summary$graphics$primary_figures == 4L,
    audit_summary$graphics$format_files == 16L
  )
  unlist(audit_summary[c(
    "checks_total",
    "checks_passed",
    "checks_failed"
  )])
}
