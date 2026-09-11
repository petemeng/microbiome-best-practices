# 多分类级组成（门→纲→目→科→属逐级）
# Run sequentially in a new working directory.
# Required packages: ggplot2, ragg, readr, scales, svglite.

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

library(ggplot2)

set.seed(20260727)
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
      plot.caption = ggplot2::element_text(colour = "#666666", hjust = 0),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#B3B3B3", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.position = "top"
    )
}

save_pub <- function(
  plot, file_base, width = 89, height = 70, units = "mm",
  dpi = 600, write_svg = TRUE, write_tiff = TRUE,
  base_family = font_pub
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot,
    width = width, height = height, units = units,
    device = grDevices::cairo_pdf, family = base_family, bg = "white"
  )
  if (isTRUE(write_svg)) {
    ggplot2::ggsave(
      paste0(file_base, ".svg"), plot,
      width = width, height = height, units = units,
      device = svglite::svglite, bg = "white"
    )
  }
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot,
    width = width, height = height, units = units, dpi = dpi,
    device = ragg::agg_png, bg = "white"
  )
  if (isTRUE(write_tiff)) {
    ggplot2::ggsave(
      paste0(file_base, ".tiff"), plot,
      width = width, height = height, units = units, dpi = dpi,
      device = ragg::agg_tiff, compression = "lzw", bg = "white"
    )
  }
  invisible(plot)
}

clean_taxon <- function(x, rank_name) {
  x <- trimws(as.character(x))
  x <- sub("^[A-Za-z]__", "", x)
  missing <- is.na(x) | x == "" | tolower(x) %in% c(
    "unassigned", "unclassified", "unknown", "na", "nan"
  )
  x[missing] <- paste("Unknown", tolower(rank_name))
  x
}

wrap_label <- function(x, width = 16L) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

read_keyed_tsv <- function(path) {
  x <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    na = character(),
    col_types = readr::cols(.default = readr::col_character())
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
  function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ),
  character(1)
)
stopifnot(identical(observed_sha256, expected_sha256))

otutab <- read_keyed_tsv(input_paths[["otutab"]])
taxonomy <- read_keyed_tsv(input_paths[["taxonomy"]])
metadata <- read_keyed_tsv(input_paths[["metadata"]])

counts <- as.matrix(otutab)
storage.mode(counts) <- "numeric"
stopifnot(
  identical(rownames(counts), rownames(taxonomy)),
  setequal(colnames(counts), rownames(metadata)),
  all(counts >= 0),
  all(counts == round(counts)),
  all(rowSums(counts) > 0),
  all(colSums(counts) > 0)
)
metadata <- metadata[colnames(counts), , drop = FALSE]

group_display <- c(
  IW = "Inland wetland",
  CW = "Coastal wetland",
  TW = "Tibetan wetland"
)
group_codes <- c("IW", "CW", "TW")
metadata$SampleID <- rownames(metadata)
metadata$GroupDisplay <- unname(group_display[metadata$Group])

sample_depths <- colSums(counts)
relative_feature <- sweep(counts, 2L, sample_depths, "/")
stopifnot(max(abs(colSums(relative_feature) - 1)) <= 1e-12)

data.frame(
  Features = nrow(counts),
  Samples = ncol(counts),
  Reads = sum(counts),
  MinimumDepth = min(sample_depths),
  MaximumDepth = max(sample_depths)
)

otutab[1:5, 1:5, drop = FALSE]
taxonomy[1:6, , drop = FALSE]
metadata[1:6, , drop = FALSE]

rank_names <- c("Phylum", "Class", "Order", "Family", "Genus")
rank_prefix <- c(
  Phylum = "P", Class = "C", Order = "O", Family = "F", Genus = "G"
)
unknown_labels <- setNames(paste("Unknown", tolower(rank_names)), rank_names)

taxonomy_clean <- as.data.frame(
  setNames(
    lapply(
      rank_names,
      function(rank_name) clean_taxon(taxonomy[[rank_name]], rank_name)
    ),
    rank_names
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(taxonomy_clean) <- rownames(taxonomy)

rank_relative <- list()
rank_counts <- list()
resolution_rows <- list()
for (rank_name in rank_names) {
  labels <- taxonomy_clean[[rank_name]]
  rank_counts[[rank_name]] <- rowsum(counts, labels, reorder = TRUE)
  rank_relative[[rank_name]] <- sweep(
    rank_counts[[rank_name]], 2L, sample_depths, "/"
  )
  stopifnot(
    max(abs(colSums(rank_relative[[rank_name]]) - 1)) <= 1e-12,
    !anyDuplicated(rownames(rank_relative[[rank_name]]))
  )

  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_order <- known_taxa[order(
    rowMeans(rel_rank[known_taxa, , drop = FALSE]),
    decreasing = TRUE
  )]
  top10 <- head(known_order, 10L)
  known_sample <- 1 - rel_rank[unknown_label, ]
  resolution_rows[[rank_name]] <- data.frame(
    Rank = rank_name,
    Labels = nrow(rel_rank),
    KnownLabels = length(known_taxa),
    KnownMean = mean(known_sample),
    KnownMinimum = min(known_sample),
    KnownMaximum = max(known_sample),
    UnknownMean = mean(rel_rank[unknown_label, ]),
    Top10CoverageMean = mean(colSums(rel_rank[top10, , drop = FALSE])),
    IncertaeSedisFeatures = sum(grepl(
      "Incertae Sedis", labels, fixed = TRUE
    )),
    stringsAsFactors = FALSE
  )
}
rank_resolution <- do.call(rbind, resolution_rows)
rank_resolution$ResolutionLossPP <- c(
  NA_real_,
  100 * head(rank_resolution$KnownMean, -1L) -
    100 * tail(rank_resolution$KnownMean, -1L)
)
rank_resolution

group_resolution_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  known_sample <- 1 - rel_rank[unknown_labels[[rank_name]], ]

  group_resolution_rows[[length(group_resolution_rows) + 1L]] <- data.frame(
    GroupCode = "ALL",
    Group = "All samples",
    Rank = rank_name,
    KnownMean = mean(known_sample),
    KnownMinimum = min(known_sample),
    KnownMaximum = max(known_sample),
    stringsAsFactors = FALSE
  )
  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    values <- known_sample[ids]
    group_resolution_rows[[length(group_resolution_rows) + 1L]] <- data.frame(
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      Rank = rank_name,
      KnownMean = mean(values),
      KnownMinimum = min(values),
      KnownMaximum = max(values),
      stringsAsFactors = FALSE
    )
  }
}
group_rank_resolution <- do.call(rbind, group_resolution_rows)

subset(group_rank_resolution, Rank == "Genus")

feature_equal_mean <- rowMeans(relative_feature)
known_matrix <- vapply(
  rank_names,
  function(rank_name) {
    taxonomy_clean[[rank_name]] != unknown_labels[[rank_name]]
  },
  logical(nrow(taxonomy_clean))
)
colnames(known_matrix) <- rank_names

complete_depth <- apply(known_matrix, 1L, function(x) {
  first_unknown <- which(!x)[1L]
  if (is.na(first_unknown)) length(x) else first_unknown - 1L
})
lineage_key <- apply(
  taxonomy_clean[, rank_names, drop = FALSE],
  1L,
  paste,
  collapse = " > "
)

lineage_gap_rows <- list()
for (i in 2:length(rank_names)) {
  parent_rank <- rank_names[[i - 1L]]
  child_rank <- rank_names[[i]]
  parent_unknown <- taxonomy_clean[[parent_rank]] == unknown_labels[[parent_rank]]
  child_unknown <- taxonomy_clean[[child_rank]] == unknown_labels[[child_rank]]
  gap <- parent_unknown & !child_unknown
  lineage_gap_rows[[length(lineage_gap_rows) + 1L]] <- data.frame(
    Transition = paste(parent_rank, "to", child_rank),
    ChildKnownParentUnknownFeatures = sum(gap),
    ChildKnownParentUnknownTaxa = length(unique(
      taxonomy_clean[[child_rank]][gap]
    )),
    EqualSampleMean = sum(feature_equal_mean[gap]),
    stringsAsFactors = FALSE
  )
}
lineage_gap_audit <- do.call(rbind, lineage_gap_rows)
lineage_gap_audit

collision_rows <- list()
for (i in 2:length(rank_names)) {
  rank_name <- rank_names[[i]]
  labels <- taxonomy_clean[[rank_name]]
  known_labels <- setdiff(unique(labels), unknown_labels[[rank_name]])
  for (taxon in known_labels) {
    mask <- labels == taxon
    parent_paths <- apply(
      taxonomy_clean[
        mask,
        rank_names[seq_len(i - 1L)],
        drop = FALSE
      ],
      1L,
      paste,
      collapse = " > "
    )
    parent_paths <- sort(unique(parent_paths))
    if (length(parent_paths) > 1L) {
      collision_rows[[length(collision_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        Taxon = taxon,
        ParentLineages = length(parent_paths),
        ParentPaths = paste(parent_paths, collapse = " | "),
        EqualSampleMean = sum(feature_equal_mean[mask]),
        stringsAsFactors = FALSE
      )
    }
  }
}
label_collision_audit <- do.call(rbind, collision_rows)
label_collision_audit

rank_top_rows <- list()
rank_group_rows <- list()
top_taxa_by_rank <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_means <- rowMeans(rel_rank[known_taxa, , drop = FALSE])
  top5 <- names(sort(known_means, decreasing = TRUE))[1:5]
  top_taxa_by_rank[[rank_name]] <- top5

  rank_top_rows[[rank_name]] <- data.frame(
    Rank = rank_name,
    TopRank = seq_along(top5),
    Taxon = top5,
    RankQualifiedTaxon = paste0(rank_prefix[[rank_name]], " · ", top5),
    OverallMean = unname(known_means[top5]),
    stringsAsFactors = FALSE
  )

  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    for (taxon in c(top5, unknown_label)) {
      values <- rel_rank[taxon, ids]
      rank_group_rows[[length(rank_group_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        GroupCode = group_code,
        Group = unname(group_display[group_code]),
        Taxon = taxon,
        RankQualifiedTaxon = paste0(rank_prefix[[rank_name]], " · ", taxon),
        MeanRelativeAbundance = mean(values),
        MeanPercent = 100 * mean(values),
        Prevalence = mean(values > 0),
        stringsAsFactors = FALSE
      )
    }
  }
}
rank_top_taxa <- do.call(rbind, rank_top_rows)
rank_taxon_group_summary <- do.call(rbind, rank_group_rows)
rank_top_taxa

top_n_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  known_taxa <- setdiff(rownames(rel_rank), unknown_label)
  known_means <- rowMeans(rel_rank[known_taxa, , drop = FALSE])
  known_order <- names(sort(known_means, decreasing = TRUE))
  branches <- list(
    `Top 5` = head(known_order, 5L),
    `Top 10` = head(known_order, 10L),
    `Top 20` = head(known_order, 20L),
    `Top 50` = head(known_order, 50L),
    `All known` = known_order
  )
  for (label in names(branches)) {
    selected <- branches[[label]]
    coverage <- colSums(rel_rank[selected, , drop = FALSE])
    top_n_rows[[length(top_n_rows) + 1L]] <- data.frame(
      Rank = rank_name,
      NLabel = label,
      SelectedN = length(selected),
      MeanCoverage = mean(coverage),
      MinimumCoverage = min(coverage),
      MaximumCoverage = max(coverage),
      UnknownMean = mean(rel_rank[unknown_label, ]),
      stringsAsFactors = FALSE
    )
  }
}
rank_topn_sensitivity <- do.call(rbind, top_n_rows)
subset(rank_topn_sensitivity, NLabel == "Top 10")

aggregation_rows <- list()
denominator_rows <- list()
for (rank_name in rank_names) {
  rel_rank <- rank_relative[[rank_name]]
  counts_rank <- rank_counts[[rank_name]]
  unknown_label <- unknown_labels[[rank_name]]
  top5 <- top_taxa_by_rank[[rank_name]]

  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    equal_vector <- rowMeans(rel_rank[, ids, drop = FALSE])
    pooled_vector <- rowSums(counts_rank[, ids, drop = FALSE]) /
      sum(counts_rank[, ids, drop = FALSE])
    aggregation_rows[[length(aggregation_rows) + 1L]] <- data.frame(
      Rank = rank_name,
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      TotalVariation = 0.5 * sum(abs(equal_vector - pooled_vector)),
      MaximumTaxonDifferencePP = 100 * max(abs(
        equal_vector - pooled_vector
      )),
      stringsAsFactors = FALSE
    )

    known_total <- 1 - rel_rank[unknown_label, ids]
    for (taxon in top5) {
      all_read_mean <- mean(rel_rank[taxon, ids])
      known_only_mean <- mean(rel_rank[taxon, ids] / known_total)
      denominator_rows[[length(denominator_rows) + 1L]] <- data.frame(
        Rank = rank_name,
        GroupCode = group_code,
        Group = unname(group_display[group_code]),
        Taxon = taxon,
        AllReadMean = all_read_mean,
        KnownOnlyMean = known_only_mean,
        InflationPP = 100 * (known_only_mean - all_read_mean),
        stringsAsFactors = FALSE
      )
    }
  }
}
aggregation_sensitivity <- do.call(rbind, aggregation_rows)
denominator_sensitivity <- do.call(rbind, denominator_rows)

denominator_sensitivity[
  which.max(denominator_sensitivity$InflationPP),
  ,
  drop = FALSE
]

complete_mask <- complete_depth == length(rank_names)
complete_rel <- rowsum(
  relative_feature[complete_mask, , drop = FALSE],
  lineage_key[complete_mask],
  reorder = TRUE
)
complete_means <- rowMeans(complete_rel)
complete_order <- names(sort(complete_means, decreasing = TRUE))
top12_lineages <- head(complete_order, 12L)

lineage_rows <- list()
ladder_rows <- list()
for (lineage_rank in seq_along(top12_lineages)) {
  key <- top12_lineages[[lineage_rank]]
  path_index <- match(key, lineage_key)
  path <- unlist(
    taxonomy_clean[path_index, rank_names, drop = FALSE],
    use.names = TRUE
  )
  leaf_values <- complete_rel[key, ]

  for (group_code in group_codes) {
    ids <- metadata$SampleID[metadata$Group == group_code]
    lineage_rows[[length(lineage_rows) + 1L]] <- data.frame(
      LineageRank = lineage_rank,
      LineageKey = key,
      GroupCode = group_code,
      Group = unname(group_display[group_code]),
      MeanRelativeAbundance = mean(leaf_values[ids]),
      Prevalence = mean(leaf_values[ids] > 0),
      stringsAsFactors = FALSE
    )
  }

  for (rank_i in seq_along(rank_names)) {
    rank_name <- rank_names[[rank_i]]
    prefix_match <- rep(TRUE, nrow(taxonomy_clean))
    for (prefix_i in seq_len(rank_i)) {
      prefix_rank <- rank_names[[prefix_i]]
      prefix_match <- prefix_match &
        taxonomy_clean[[prefix_rank]] == path[[prefix_rank]]
    }
    node_mean <- sum(feature_equal_mean[prefix_match])
    ladder_rows[[length(ladder_rows) + 1L]] <- data.frame(
      LineageRank = lineage_rank,
      LineageKey = key,
      LeafGenus = path[["Genus"]],
      Rank = rank_name,
      Taxon = path[[rank_name]],
      NodePercent = 100 * node_mean,
      LeafPercent = 100 * mean(leaf_values),
      stringsAsFactors = FALSE
    )
  }
}
dominant_complete_lineages <- do.call(rbind, lineage_rows)
lineage_ladder <- do.call(rbind, ladder_rows)

data.frame(
  CompleteLineages = nrow(complete_rel),
  CompleteLineageCoverage = sum(complete_means),
  Top12Coverage = sum(complete_means[top12_lineages])
)

result_dir <- "results/27-multirank-composition"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_lineage_map <- data.frame(
  FeatureID = rownames(taxonomy_clean),
  taxonomy_clean,
  LineageKey = lineage_key,
  CompleteDepth = complete_depth,
  EqualSampleMean = feature_equal_mean,
  check.names = FALSE
)

readr::write_tsv(
  taxonomy_lineage_map,
  file.path(result_dir, "article-taxonomy-lineage-map.tsv")
)
readr::write_tsv(
  rank_resolution,
  file.path(result_dir, "article-rank-resolution.tsv")
)
readr::write_tsv(
  group_rank_resolution,
  file.path(result_dir, "article-group-rank-resolution.tsv")
)
readr::write_tsv(
  lineage_gap_audit,
  file.path(result_dir, "article-lineage-gap-audit.tsv")
)
readr::write_tsv(
  label_collision_audit,
  file.path(result_dir, "article-label-collision-audit.tsv")
)
readr::write_tsv(
  rank_top_taxa,
  file.path(result_dir, "article-rank-top-taxa.tsv")
)
readr::write_tsv(
  rank_taxon_group_summary,
  file.path(result_dir, "article-rank-taxon-group-summary.tsv")
)
readr::write_tsv(
  rank_topn_sensitivity,
  file.path(result_dir, "article-rank-topn-sensitivity.tsv")
)
readr::write_tsv(
  lineage_ladder,
  file.path(result_dir, "article-lineage-ladder.tsv")
)
readr::write_tsv(
  aggregation_sensitivity,
  file.path(result_dir, "article-aggregation-sensitivity.tsv")
)
readr::write_tsv(
  denominator_sensitivity,
  file.path(result_dir, "article-denominator-sensitivity.tsv")
)

group_colours <- c(
  `All samples` = "#1F2937",
  `Inland wetland` = "#E69F00",
  `Coastal wetland` = "#0072B2",
  `Tibetan wetland` = "#009E73"
)
cascade_data <- group_rank_resolution
cascade_data$Rank <- factor(cascade_data$Rank, levels = rank_names)
cascade_data$Group <- factor(
  cascade_data$Group,
  levels = c("All samples", unname(group_display[group_codes]))
)

cascade_plot <- ggplot2::ggplot(
  cascade_data,
  ggplot2::aes(
    x = Rank,
    y = KnownMean,
    colour = Group,
    linetype = Group,
    group = Group
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(
    values = c(
      `All samples` = "22",
      `Inland wetland` = "solid",
      `Coastal wetland` = "solid",
      `Tibetan wetland` = "solid"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Taxonomic resolution falls as the requested rank becomes finer",
    subtitle = "Known fraction uses all reads and equal sample weight within each wetland group",
    x = NULL,
    y = "Reads assigned at rank",
    colour = NULL,
    linetype = NULL,
    caption = paste(
      "Genus-level resolution: inland 29.0%, coastal 38.4%, Tibetan 43.3%;\n",
      "differences in annotation completeness remain visible."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

save_pub(
  cascade_plot,
  "figures/27-rank-resolution-cascade",
  width = 150,
  height = 96
)
cascade_plot

bubble_data <- rank_taxon_group_summary
bubble_data$Rank <- factor(bubble_data$Rank, levels = rank_names)
bubble_data$Group <- factor(
  bubble_data$Group,
  levels = unname(group_display[group_codes])
)
bubble_levels <- unlist(lapply(rank_names, function(rank_name) {
  c(
    paste0(rank_prefix[[rank_name]], " · ", top_taxa_by_rank[[rank_name]]),
    paste0(rank_prefix[[rank_name]], " · ", unknown_labels[[rank_name]])
  )
}), use.names = FALSE)
bubble_data$RankQualifiedTaxon <- factor(
  bubble_data$RankQualifiedTaxon,
  levels = rev(bubble_levels)
)

bubble_plot <- ggplot2::ggplot(
  bubble_data,
  ggplot2::aes(
    x = Group,
    y = RankQualifiedTaxon,
    size = MeanPercent,
    colour = Prevalence
  )
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::facet_wrap(~Rank, ncol = 2, scales = "free_y") +
  ggplot2::scale_size_area(max_size = 11, breaks = c(1, 5, 20, 60)) +
  ggplot2::scale_colour_gradientn(
    colours = c("#D9EAF7", "#56B4E9", "#0072B2", "#003B5C"),
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      `Inland wetland` = "IW",
      `Coastal wetland` = "CW",
      `Tibetan wetland` = "TW"
    )
  ) +
  ggplot2::labs(
    title = "The same Top-5 rule yields different views at each rank",
    subtitle = "Global Top 5 named taxa plus rank-specific Unknown; one fixed list is reused across groups",
    x = "Wetland group (IW / CW / TW)",
    y = NULL,
    size = "Mean abundance (%)",
    colour = "Prevalence",
    caption = paste(
      "IW, inland; CW, coastal; TW, Tibetan wetland. P/C/O/F/G prefixes prevent",
      "equal-looking names at different ranks from being merged."
    )
  ) +
  theme_pub(base_size = 7.7) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = ggplot2::element_line(
      colour = "#ECECEC", linewidth = 0.25
    ),
    panel.spacing.x = grid::unit(1.2, "mm"),
    strip.text = ggplot2::element_text(size = 8)
  )

save_pub(
  bubble_plot,
  "figures/27-multirank-bubble",
  width = 183,
  height = 198
)
bubble_plot

ladder_data <- lineage_ladder
ladder_data$Rank <- factor(ladder_data$Rank, levels = rank_names)
leaf_order <- unique(ladder_data$LeafGenus[order(ladder_data$LineageRank)])
ladder_data$LeafGenus <- factor(
  ladder_data$LeafGenus,
  levels = rev(leaf_order)
)
ladder_data$CellLabel <- wrap_label(ladder_data$Taxon, width = 15L)
ladder_dark <- ladder_data$NodePercent >= 9

ladder_plot <- ggplot2::ggplot(
  ladder_data,
  ggplot2::aes(x = Rank, y = LeafGenus, fill = NodePercent)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(
    data = ladder_data[!ladder_dark, , drop = FALSE],
    ggplot2::aes(label = CellLabel),
    colour = "#17202A",
    family = font_pub,
    size = 2.05,
    lineheight = 0.88
  ) +
  ggplot2::geom_text(
    data = ladder_data[ladder_dark, , drop = FALSE],
    ggplot2::aes(label = CellLabel),
    colour = "white",
    family = font_pub,
    size = 2.05,
    lineheight = 0.88
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    trans = "sqrt",
    labels = scales::label_number(accuracy = 0.1, suffix = "%")
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) wrap_label(x, width = 21L)
  ) +
  ggplot2::labs(
    title = "A genus is a complete path, not an isolated display name",
    subtitle = "Top 12 complete genus lineages cover 10.19% of all reads; tile colour is prefix-node abundance",
    x = NULL,
    y = "Leaf genus",
    fill = "Node abundance",
    caption = paste(
      "Each row traces one Phylum-to-Genus path. Prefix-node percentages repeat",
      "context and must not be added across columns."
    )
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.ticks = ggplot2::element_blank()
  )

save_pub(
  ladder_plot,
  "figures/27-lineage-ladder",
  width = 183,
  height = 168
)
ladder_plot

rank_colours <- c(
  Phylum = "#0072B2",
  Class = "#E69F00",
  Order = "#009E73",
  Family = "#CC79A7",
  Genus = "#D55E00"
)
topn_plot_data <- rank_topn_sensitivity
topn_plot_data$Rank <- factor(topn_plot_data$Rank, levels = rank_names)
topn_plot_data$NLabel <- factor(
  topn_plot_data$NLabel,
  levels = c("Top 5", "Top 10", "Top 20", "Top 50", "All known")
)
top10_labels <- topn_plot_data[
  topn_plot_data$NLabel == "Top 10",
  ,
  drop = FALSE
]
top10_labels$Label <- scales::percent(
  top10_labels$MeanCoverage,
  accuracy = 0.1
)

topn_plot <- ggplot2::ggplot(
  topn_plot_data,
  ggplot2::aes(x = NLabel, y = MeanCoverage, colour = Rank, group = Rank)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::geom_text(
    data = top10_labels,
    ggplot2::aes(label = Label),
    nudge_x = 0.14,
    hjust = 0,
    size = 2.4,
    show.legend = FALSE,
    family = font_pub
  ) +
  ggplot2::scale_colour_manual(values = rank_colours, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.03),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Top 10 does not mean the same coverage at every rank",
    subtitle = "Coverage is measured against all reads after per-sample closure",
    x = "Display branch",
    y = "Mean covered abundance",
    colour = "Rank",
    caption = paste(
      "Top-10 coverage falls from 94.8% at phylum to 9.5% at genus;\n",
      "All known retains the rank-specific unclassified fraction."
    )
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

save_pub(
  topn_plot,
  "figures/27-topn-coverage",
  width = 160,
  height = 102
)
topn_plot

expected_known <- c(
  Phylum = 0.996536217607761,
  Class = 0.967978544271532,
  Order = 0.855727280469946,
  Family = 0.679627539342075,
  Genus = 0.368943886930081
)
stopifnot(
  max(abs(rank_resolution$KnownMean - expected_known)) <= 1e-12,
  nrow(label_collision_audit) == 2L,
  lineage_gap_audit$ChildKnownParentUnknownFeatures[4] == 324L,
  nrow(complete_rel) == 754L,
  abs(sum(complete_means) - 0.343761113771336) <= 1e-12,
  abs(sum(complete_means[top12_lineages]) - 0.101876289588930) <= 1e-12,
  nrow(rank_topn_sensitivity) == 25L
)

audit_path <- file.path(result_dir, "validation-checks.tsv")
if (file.exists(audit_path)) {
  validation_checks <- readr::read_tsv(
    audit_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    nrow(validation_checks) == 209L,
    all(validation_checks$status == "PASS")
  )
  data.frame(
    Validation = "Repository QA",
    Passed = sum(validation_checks$status == "PASS"),
    Total = nrow(validation_checks)
  )
} else {
  data.frame(
    Validation = "Inline reproducibility checks",
    Passed = 7L,
    Total = 7L
  )
}
