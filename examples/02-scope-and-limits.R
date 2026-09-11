# 16S 能回答什么、回答不了什么（vs 宏基因组 / 其他组学）
# Run sequentially in a new working directory.
# Required packages: ggplot2, scales.

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
  all(colSums(otutab) > 0)
)

dim(otutab)
otutab[1:5, 1:5]
taxonomy[1:5, , drop = FALSE]
head(metadata)

rank_order <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)
stopifnot(all(rank_order %in% colnames(taxonomy)))

clean_taxon <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[A-Za-z]__", "", x)
  x
}

is_assigned <- function(x) {
  !is.na(x) &
    nzchar(x) &
    !grepl(
      "^(unassigned|unclassified|uncultured|unknown|norank|NA)$",
      x,
      ignore.case = TRUE
    )
}

taxonomy_clean <- taxonomy
taxonomy_clean[rank_order] <- lapply(
  taxonomy_clean[rank_order],
  clean_taxon
)

rank_summary <- do.call(
  rbind,
  lapply(rank_order, function(rank_name) {
    values <- taxonomy_clean[[rank_name]]
    assigned <- is_assigned(values)
    data.frame(
      Rank = rank_name,
      AssignedFeatures = sum(assigned),
      TotalFeatures = length(values),
      AssignedPercent = 100 * mean(assigned),
      DistinctLabels = length(unique(values[assigned])),
      stringsAsFactors = FALSE
    )
  })
)

rank_summary$Rank <- factor(
  rank_summary$Rank,
  levels = rev(rank_order)
)

rank_summary
stopifnot(
  all(rank_summary$AssignedFeatures <= rank_summary$TotalFeatures),
  all(diff(rank_summary$AssignedPercent) <= 1e-10)
)

questions <- c(
  "Community\nstructure",
  "Species /\nstrain",
  "Gene\npotential",
  "Gene\nactivity",
  "Metabolites",
  "Absolute\nload"
)

scope_row <- function(assay, status, label) {
  stopifnot(
    length(status) == length(questions),
    length(label) == length(questions)
  )
  data.frame(
    Assay = assay,
    Question = questions,
    Status = status,
    Label = label,
    stringsAsFactors = FALSE
  )
}

scope_matrix <- do.call(
  rbind,
  list(
    scope_row(
      "16S amplicon",
      c(
        "Direct", "Partial / inferred", "Partial / inferred",
        "Not measured", "Not measured", "Not measured"
      ),
      c("Direct", "Limited", "Predicted", "No", "No", "No")
    ),
    scope_row(
      "Shotgun metagenome",
      c(
        "Direct", "Direct", "Direct",
        "Not measured", "Not measured", "Partial / inferred"
      ),
      c("Direct", "Strong", "Direct", "No", "No", "Relative")
    ),
    scope_row(
      "Metatranscriptome",
      c(
        "Partial / inferred", "Partial / inferred", "Partial / inferred",
        "Direct", "Not measured", "Partial / inferred"
      ),
      c("Biased", "Partial", "Expressed", "Direct", "No", "Relative")
    ),
    scope_row(
      "Metabolome",
      c(
        "Not measured", "Not measured", "Not measured",
        "Partial / inferred", "Direct", "Partial / inferred"
      ),
      c("No", "No", "No", "Downstream", "Direct", "Relative")
    ),
    scope_row(
      "Absolute reference\n(qPCR / flow / spike-in)",
      c(
        "Partial / inferred", "Partial / inferred", "Not measured",
        "Not measured", "Not measured", "Direct"
      ),
      c("Targeted", "Targeted", "No", "No", "No", "Direct")
    )
  )
)

scope_matrix$Question <- factor(
  scope_matrix$Question,
  levels = questions
)
scope_matrix$Assay <- factor(
  scope_matrix$Assay,
  levels = rev(unique(scope_matrix$Assay))
)
scope_matrix$Status <- factor(
  scope_matrix$Status,
  levels = c("Direct", "Partial / inferred", "Not measured")
)

scope_matrix

rank_palette <- c(
  Kingdom = "#1B4965",
  Phylum = "#245F73",
  Class = "#2C7DA0",
  Order = "#3A8D9F",
  Family = "#2A9D8F",
  Genus = "#B56576",
  Species = "#7A1F3D"
)

rank_summary$Label <- sprintf(
  "%.1f%% · n=%s · %s labels",
  rank_summary$AssignedPercent,
  format(rank_summary$AssignedFeatures, big.mark = ","),
  format(rank_summary$DistinctLabels, big.mark = ",")
)

p_taxonomy_resolution <- ggplot(
  rank_summary,
  aes(Rank, AssignedPercent, fill = Rank)
) +
  geom_col(width = 0.72) +
  geom_text(
    data = subset(rank_summary, AssignedPercent >= 50),
    aes(y = AssignedPercent - 2, label = Label),
    hjust = 1,
    size = 3.0,
    color = "white"
  ) +
  geom_text(
    data = subset(rank_summary, AssignedPercent < 50),
    aes(y = AssignedPercent + 2, label = Label),
    hjust = 0,
    size = 3.0,
    color = "grey15"
  ) +
  scale_fill_manual(values = rank_palette, guide = "none") +
  scale_y_continuous(
    limits = c(0, 112),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "Taxonomic annotation depth in a real 16S dataset",
    subtitle = paste0(
      format(nrow(taxonomy), big.mark = ","),
      " features · empty and placeholder labels excluded"
    ),
    x = NULL,
    y = "Features with a non-empty taxonomic label"
  ) +
  theme_pub(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(8, 12, 8, 8)
  )

p_taxonomy_resolution
save_pub(
  p_taxonomy_resolution,
  "figures/02-taxonomy-resolution",
  width = 165,
  height = 105
)

scope_colors <- c(
  "Direct" = "#009E73",
  "Partial / inferred" = "#E69F00",
  "Not measured" = "#D9D9D9"
)

p_assay_scope <- ggplot(
  scope_matrix,
  aes(Question, Assay, fill = Status)
) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(
    aes(label = Label),
    size = 3.0,
    color = "grey10",
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = scope_colors,
    drop = FALSE,
    name = "Measurement"
  ) +
  labs(
    title = "What each microbiome assay can measure",
    subtitle = "Direct measurement does not imply freedom from bias or confounding",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      lineheight = 0.9
    ),
    axis.text.y = element_text(
      color = "black",
      face = "bold",
      lineheight = 0.9
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p_assay_scope
save_pub(
  p_assay_scope,
  "figures/02-assay-scope-map",
  width = 180,
  height = 105
)
