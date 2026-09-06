# 群落组成图：堆叠/气泡/冲积/热图
# Run sequentially in a new working directory.
# Required packages: ggplot2, ggalluvial, ragg, svglite.

library(ggplot2)
library(ggalluvial)
set.seed(20260726)
options(timeout = 300)

data_url <- paste0(
  "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/",
  "f394362dda37f04fe62a5ba245cd00b94cbb84b9/data/small/"
)
dir.create("data/small", recursive = TRUE, showWarnings = FALSE)
for (file in c("otutab.tsv", "taxonomy.tsv", "metadata.tsv")) {
  path <- file.path("data/small", file)
  if (!file.exists(path)) download.file(paste0(data_url, file), path, mode = "wb")
}

otutab <- read.delim("data/small/otutab.tsv", row.names = 1, check.names = FALSE)
taxonomy <- read.delim("data/small/taxonomy.tsv", row.names = 1, check.names = FALSE)
metadata <- read.delim("data/small/metadata.tsv", row.names = 1, check.names = FALSE)

counts <- as.matrix(otutab)
stopifnot(
  !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
  setequal(rownames(counts), rownames(taxonomy)),
  setequal(colnames(counts), rownames(metadata)),
  all(is.finite(counts)), all(counts >= 0),
  all(counts == round(counts)), all(colSums(counts) > 0)
)
taxonomy <- taxonomy[rownames(counts), , drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]

dim(counts)
head(otutab[, 1:4], 3)
head(taxonomy[, c("Phylum", "Genus")], 3)
head(metadata[, c("Group", "Type")], 3)

group_codes <- c("IW", "CW", "TW")
group_labels <- c(IW = "Inland wetland", CW = "Coastal wetland", TW = "Tibetan wetland")
region_labels <- c(NE = "Northeast", NW = "Northwest", NC = "North China",
                   YML = "Yangtze middle-lower", SC = "Southern coast",
                   QTP = "Qinghai-Tibet Plateau")
metadata$SampleID <- rownames(metadata)
metadata$GroupDisplay <- factor(group_labels[metadata$Group], levels = group_labels)
metadata$Region <- unname(region_labels[metadata$Type])
sample_order <- metadata$SampleID[order(
  metadata$GroupDisplay, metadata$Region,
  as.integer(sub("^S", "", metadata$SampleID))
)]
table(metadata$GroupDisplay)

theme_set(theme_bw(base_size = 9, base_family = "sans") +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  plot.caption = element_text(hjust = 0)))

clean_taxon <- function(x, rank) {
  x <- sub("^[A-Za-z]__", "", trimws(as.character(x)))
  missing <- is.na(x) | x == "" | tolower(x) %in% c("unassigned", "unclassified", "unknown", "na", "nan")
  x[missing] <- paste("Unknown", tolower(rank))
  x
}
aggregate_rank <- function(counts, taxonomy, rank) {
  out <- rowsum(counts, clean_taxon(taxonomy[[rank]], rank), reorder = TRUE)
  unknown <- paste("Unknown", tolower(rank))
  if (!unknown %in% rownames(out)) {
    out <- rbind(out, matrix(0, 1, ncol(out), dimnames = list(unknown, colnames(out))))
  }
  out[order(rownames(out)), , drop = FALSE]
}
sample_depths <- colSums(counts)
phylum_counts <- aggregate_rank(counts, taxonomy, "Phylum")
phylum_rel <- sweep(phylum_counts, 2, sample_depths, "/")
stopifnot(max(abs(colSums(phylum_rel) - 1)) < 1e-12)

known_phyla <- setdiff(rownames(phylum_rel), "Unknown phylum")
known_order <- known_phyla[order(rowMeans(phylum_rel)[known_phyla], decreasing = TRUE)]
top10_phyla <- head(known_order, 10)
display_labels <- c(top10_phyla, "Other known phyla", "Unknown phylum")
display_group <- ifelse(rownames(phylum_rel) %in% top10_phyla, rownames(phylum_rel),
                        ifelse(rownames(phylum_rel) == "Unknown phylum",
                               "Unknown phylum", "Other known phyla"))
phylum_display <- rowsum(phylum_rel, display_group, reorder = FALSE)[display_labels, , drop = FALSE]

matrix_to_long <- function(x, value_name) {
  out <- data.frame(Taxon = rep(rownames(x), ncol(x)),
                    SampleID = rep(colnames(x), each = nrow(x)), Value = as.vector(x))
  names(out)[3] <- value_name
  out
}
stack_data <- matrix_to_long(phylum_display, "RelativeAbundance")
stack_data$SampleID <- factor(stack_data$SampleID, levels = sample_order)
stack_data$Taxon <- factor(stack_data$Taxon, levels = rev(display_labels))
stack_data$Group <- metadata[as.character(stack_data$SampleID), "GroupDisplay"]
taxon_palette <- c(Proteobacteria = "#0072B2", Chloroflexi = "#E69F00",
                    Bacteroidetes = "#009E73", Acidobacteria = "#D55E00",
                    Actinobacteria = "#CC79A7", Firmicutes = "#56B4E9",
                    Verrucomicrobia = "#F0E442", Planctomycetes = "#8C564B",
                    Gemmatimonadetes = "#17BECF", Nitrospirae = "#9467BD",
                    "Other known phyla" = "#B8B8B8", "Unknown phylum" = "#2F2F2F")

stack_plot <- ggplot(stack_data, aes(SampleID, RelativeAbundance, fill = Taxon)) +
  geom_col(width = 1) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = taxon_palette, drop = FALSE) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(title = "Phylum composition across 90 wetland samples",
       x = "Samples ordered by wetland group and region", y = "Relative abundance", fill = "Phylum") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        panel.spacing.x = grid::unit(1.5, "mm"))
stack_plot

group_mean <- function(x) {
  sapply(group_codes, function(g) rowMeans(x[, metadata$Group == g, drop = FALSE]))
}
phylum_group_mean <- group_mean(phylum_display)
round(100 * phylum_group_mean[c("Proteobacteria", "Acidobacteria", "Bacteroidetes"), ], 2)

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_matrices <- setNames(lapply(rank_names, function(rank) {
  sweep(aggregate_rank(counts, taxonomy, rank), 2, sample_depths, "/")
}), rank_names)
unknown_percent <- sapply(rank_names, function(rank) {
  100 * mean(rank_matrices[[rank]][paste("Unknown", tolower(rank)), ])
})
round(unknown_percent[c("Phylum", "Genus", "Species")], 2)

top_n_sensitivity <- sapply(c(5, 8, 10, 15, 20), function(n) {
  coverage <- colSums(phylum_rel[head(known_order, n), , drop = FALSE])
  100 * c(Mean = mean(coverage), Minimum = min(coverage), Maximum = max(coverage))
})
colnames(top_n_sensitivity) <- paste0("Top", c(5, 8, 10, 15, 20))
round(top_n_sensitivity, 2)

genus_rel <- rank_matrices[["Genus"]]
unknown_genus <- "Unknown genus"
known_genera <- setdiff(rownames(genus_rel), unknown_genus)
genus_order <- known_genera[order(rowMeans(genus_rel)[known_genera], decreasing = TRUE)]
top20_genera <- head(genus_order, 20)
top15_genera <- head(genus_order, 15)
ids <- metadata$SampleID[metadata$Group == "IW"]
known_fraction <- 1 - genus_rel[unknown_genus, ids]
stopifnot(all(known_fraction > 0))
denominator_sensitivity <- 100 * c(
  AllReads = mean(genus_rel["Candidatus Solibacter", ids]),
  ClassifiedReads = mean(genus_rel["Candidatus Solibacter", ids] / known_fraction)
)
round(denominator_sensitivity, 2)

pooled_phylum <- sapply(group_codes, function(g) {
  ids <- metadata$SampleID[metadata$Group == g]
  weighted <- sweep(phylum_display[, ids, drop = FALSE], 2, sample_depths[ids], "*")
  rowSums(weighted) / sum(sample_depths[ids])
})
round(100 * max(abs(phylum_group_mean - pooled_phylum)), 2)

shown_genera <- c(top15_genera, unknown_genus)
genus_means <- group_mean(genus_rel[shown_genera, , drop = FALSE])
genus_prevalence <- group_mean(genus_rel[shown_genera, , drop = FALSE] > 0)
bubble_data <- matrix_to_long(genus_means, "MeanRelativeAbundance")
bubble_data$MeanPercent <- 100 * bubble_data$MeanRelativeAbundance
bubble_data$Prevalence <- as.vector(genus_prevalence)
bubble_data$Group <- factor(group_labels[bubble_data$SampleID], levels = group_labels)
bubble_data$Taxon <- factor(bubble_data$Taxon, levels = rev(shown_genera))

bubble_plot <- ggplot(bubble_data, aes(Group, Taxon, size = MeanPercent, colour = Prevalence)) +
  geom_point(alpha = 0.9) +
  scale_size_area(max_size = 13, breaks = c(1, 5, 20, 60)) +
  scale_y_discrete(expand = expansion(add = c(1, 0.8))) +
  scale_colour_gradientn(colours = c("#D9EAF7", "#56B4E9", "#0072B2", "#003B5C"),
                         limits = c(0, 1), labels = scales::label_percent()) +
  scale_x_discrete(labels = c("Inland wetland" = "Inland\nwetland",
                              "Coastal wetland" = "Coastal\nwetland",
                              "Tibetan wetland" = "Tibetan\nwetland")) +
  labs(title = "Genus abundance and prevalence", x = NULL, y = NULL,
       size = "Mean abundance (%)", colour = "Prevalence") +
  theme(legend.position = "right")
bubble_plot

heatmap_matrix <- log1p(10000 * genus_rel[top20_genera, sample_order, drop = FALSE])
row_order <- rownames(heatmap_matrix)[hclust(dist(heatmap_matrix), method = "complete")$order]
heatmap_data <- matrix_to_long(heatmap_matrix, "Log1pAbundance")
heatmap_data$SampleID <- factor(heatmap_data$SampleID, levels = sample_order)
heatmap_data$Taxon <- factor(heatmap_data$Taxon, levels = rev(row_order))
heatmap_data$Group <- metadata[as.character(heatmap_data$SampleID), "GroupDisplay"]
heatmap_plot <- ggplot(heatmap_data, aes(SampleID, Taxon, fill = Log1pAbundance)) +
  geom_tile() +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")) +
  labs(title = "Genus composition varies within wetland groups",
       x = "Samples ordered by wetland group and region", y = NULL,
       fill = "log1p(relative\nabundance x 10,000)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), panel.spacing.x = grid::unit(1.5, "mm"))
heatmap_plot

ids <- metadata$SampleID[metadata$Group == "TW"]
round(100 * c(Mean = mean(genus_rel["Pedobacter", ids]),
               Minimum = min(genus_rel["Pedobacter", ids]),
               Maximum = max(genus_rel["Pedobacter", ids])), 2)

flow_data <- matrix_to_long(phylum_group_mean, "MeanRelativeAbundance")
flow_data$Group <- factor(group_labels[flow_data$SampleID], levels = group_labels)
flow_data$Taxon <- factor(flow_data$Taxon, levels = display_labels)
alluvial_plot <- ggplot(flow_data, aes(y = MeanRelativeAbundance, axis1 = Group, axis2 = Taxon)) +
  geom_alluvium(aes(fill = Taxon), width = 0.08, alpha = 0.8) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.08, colour = "white") +
  geom_text(stat = "stratum",
            aes(label = after_stat(ifelse(x == 1, sub(" wetland", "\nwetland", stratum), ""))),
            size = 2.5, hjust = 1, nudge_x = -0.06) +
  scale_x_discrete(limits = c("Wetland group", "Phylum"), expand = c(0.35, 0.05)) +
  scale_fill_manual(values = c(taxon_palette, setNames(rep("grey85", 3), group_labels)),
                    breaks = display_labels, name = "Phylum") +
  labs(title = "Wetland groups share dominant phyla",
       caption = "Group means; ribbons do not imply migration.") +
  theme_void(base_family = "sans", base_size = 9) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"),
        axis.text.x = element_text(), plot.margin = margin(10, 5, 10, 5))
alluvial_plot

plots <- list("26-phylum-stacked" = stack_plot, "26-genus-bubble" = bubble_plot,
              "26-genus-heatmap" = heatmap_plot, "26-group-phylum-alluvial" = alluvial_plot)
plot_sizes <- data.frame(width = c(183, 140, 183, 183), height = c(112, 118, 122, 112),
                         row.names = names(plots))
devices <- list(pdf = grDevices::cairo_pdf, svg = svglite::svglite,
                png = ragg::agg_png, tiff = ragg::agg_tiff)
dir.create("figures", showWarnings = FALSE)
for (name in names(plots)) {
  for (format in names(devices)) {
    args <- list(filename = file.path("figures", paste0(name, ".", format)),
                 plot = plots[[name]], device = devices[[format]],
                 width = plot_sizes[name, "width"], height = plot_sizes[name, "height"],
                 units = "mm", dpi = 600, bg = "white")
    if (format == "tiff") args$compression <- "lzw"
    do.call(ggsave, args)
  }
}

dir.create("results/26-community-composition", recursive = TRUE, showWarnings = FALSE)
write.table(phylum_rel, "results/26-community-composition/reader-phylum-relative.tsv",
             sep = "\t", quote = FALSE, col.names = NA)
write.table(genus_rel, "results/26-community-composition/reader-genus-relative.tsv",
             sep = "\t", quote = FALSE, col.names = NA)
