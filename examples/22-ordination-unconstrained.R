# 非约束排序：PCoA/NMDS/CLR-PCA（复现顶刊）
# Run sequentially in a new working directory.
# Required packages: ggplot2, scales, vegan.

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

min_prevalence <- 0.05
min_samples <- ceiling(min_prevalence * ncol(otutab))
keep_feature <- rowSums(otutab > 0) >= min_samples

comm <- t(otutab[keep_feature, , drop = FALSE])
comm_rel <- sweep(comm, 1, rowSums(comm), "/")
metadata <- metadata[rownames(comm_rel), , drop = FALSE]

stopifnot(
  all(abs(rowSums(comm_rel) - 1) < 1e-10),
  identical(rownames(comm_rel), rownames(metadata))
)

dist_bray <- vegan::vegdist(comm_rel, method = "bray")
c(
  samples = nrow(comm_rel),
  retained_features = ncol(comm_rel),
  original_features = nrow(otutab)
)

pcoa_fit <- stats::cmdscale(
  dist_bray,
  k = 2,
  eig = TRUE,
  add = TRUE
)

positive_eigenvalues <- pcoa_fit$eig[pcoa_fit$eig > 0]
pcoa_variance <- 100 * pcoa_fit$eig[1:2] / sum(positive_eigenvalues)

pcoa_df <- as.data.frame(pcoa_fit$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df$Salinity <- factor(
  metadata[pcoa_df$SampleID, "Saline"],
  levels = c("Non-saline soil", "Saline soil")
)

head(pcoa_df)
round(pcoa_variance, 2)

set.seed(20260718)
permanova <- vegan::adonis2(
  dist_bray ~ Salinity,
  data = pcoa_df,
  permutations = 999
)

dispersion_model <- vegan::betadisper(
  dist_bray,
  group = pcoa_df$Salinity,
  type = "median",
  bias.adjust = TRUE
)
set.seed(20260718)
dispersion_test <- vegan::permutest(
  dispersion_model,
  permutations = 999
)

permanova
dispersion_test

permanova_r2 <- unname(permanova$R2[1])
permanova_p <- unname(permanova$`Pr(>F)`[1])
dispersion_p <- unname(dispersion_test$tab$`Pr(>F)`[1])

set.seed(20260718)
nmds_fit <- vegan::metaMDS(
  dist_bray,
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

nmds_df <- as.data.frame(
  vegan::scores(nmds_fit, display = "sites")
)
nmds_df$SampleID <- rownames(nmds_df)
nmds_df$Salinity <- pcoa_df[nmds_df$SampleID, "Salinity"]

p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = Salinity)) +
  geom_point(size = 2.3, alpha = 0.85) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.7) +
  scale_color_pub(name = "Soil salinity") +
  labs(
    title = sprintf("Bray-Curtis NMDS (stress = %.3f)", nmds_fit$stress),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_pub()

p_nmds
nmds_fit$stress

clr_min_prevalence <- 0.20
keep_clr <- rowSums(otutab > 0) >=
  ceiling(clr_min_prevalence * ncol(otutab))
clr_counts <- t(otutab[keep_clr, , drop = FALSE])

zero_replacement <- 0.5
log_replaced <- log(clr_counts + zero_replacement)
clr_matrix <- log_replaced - rowMeans(log_replaced)

clr_pca_fit <- stats::prcomp(
  clr_matrix,
  center = TRUE,
  scale. = FALSE,
  rank. = 2
)
clr_variance <- 100 * summary(clr_pca_fit)$importance[2, 1:2]

clr_df <- as.data.frame(clr_pca_fit$x[, 1:2, drop = FALSE])
clr_df$SampleID <- rownames(clr_df)
clr_df$Salinity <- pcoa_df[clr_df$SampleID, "Salinity"]

p_clr <- ggplot(clr_df, aes(PC1, PC2, color = Salinity)) +
  geom_point(size = 2.3, alpha = 0.85) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.7) +
  scale_color_pub(name = "Soil salinity") +
  labs(
    title = "CLR-PCA sensitivity analysis",
    x = sprintf("PC1 (%.1f%%)", clr_variance[1]),
    y = sprintf("PC2 (%.1f%%)", clr_variance[2])
  ) +
  theme_pub()

p_clr

format_p <- function(p) {
  if (is.na(p)) {
    return("= NA")
  }
  if (p < 0.001) {
    return("< 0.001")
  }
  sprintf("= %.3f", p)
}

stat_label <- paste0(
  "PERMANOVA  R² = ", sprintf("%.3f", permanova_r2),
  ", P ", format_p(permanova_p),
  "\nDispersion  P ", format_p(dispersion_p)
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(PCoA1, PCoA2, color = Salinity)
) +
  stat_ellipse(
    type = "t",
    level = 0.95,
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  geom_point(size = 2.7, alpha = 0.88) +
  scale_color_pub(name = "Soil salinity") +
  labs(
    subtitle = stat_label,
    x = sprintf("PCoA 1 (%.1f%%)", pcoa_variance[1]),
    y = sprintf("PCoA 2 (%.1f%%)", pcoa_variance[2])
  ) +
  coord_equal() +
  theme_pub(base_size = 11) +
  theme(
    plot.subtitle = element_text(size = 9.5, lineheight = 1.05),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

p_pcoa
save_pub(
  p_pcoa,
  "figures/22-pcoa-bray",
  width = 120,
  height = 95
)
