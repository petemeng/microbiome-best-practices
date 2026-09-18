# 非约束排序：如何选择和解释 PCoA、NMDS 与 CLR-PCA
# Run sequentially in a new working directory.
# Required packages: ggplot2, scales, vegan, knitr.

pal_pub <- c("#0072B2", "#D55E00", "#009E73")
scale_color_pub <- function(...) ggplot2::scale_color_manual(values = pal_pub, ...)
scale_fill_pub <- function(...) ggplot2::scale_fill_manual(values = pal_pub, ...)
theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    )
}
save_pub <- function(plot, file_base, width = 150, height = 110) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  for (extension in c("pdf", "png", "tiff")) {
    extra <- switch(extension,
      pdf = list(device = grDevices::cairo_pdf),
      png = list(dpi = 300),
      tiff = list(dpi = 300, compression = "lzw"))
    do.call(ggplot2::ggsave, c(list(
      filename = paste0(file_base, ".", extension), plot = plot,
      width = width, height = height, units = "mm", bg = "white"
    ), extra))
  }
  invisible(plot)
}

options(timeout = 600)
data_url <- "https://raw.githubusercontent.com/petemeng/microbiome-best-practices/3cb6a817e0c73ab7ecbec1cedc1ad80cbb9ecfaa/"
input_files <- paste0("data/small/", c("otutab.tsv", "taxonomy.tsv", "metadata.tsv"))
for (path in input_files) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) download.file(paste0(data_url, path), path, mode = "wb", quiet = TRUE)
}

library(vegan)
library(ggplot2)
otutab <- as.matrix(read.delim("data/small/otutab.tsv", row.names = 1, check.names = FALSE))
taxonomy <- read.delim("data/small/taxonomy.tsv", row.names = 1, check.names = FALSE)
metadata <- read.delim("data/small/metadata.tsv", row.names = 1, check.names = FALSE)
storage.mode(otutab) <- "numeric"
stopifnot(
  !anyDuplicated(rownames(otutab)), !anyDuplicated(colnames(otutab)),
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(is.finite(otutab)), all(otutab >= 0), all(colSums(otutab) > 0)
)
otutab[1:5, 1:5]
head(taxonomy, 3)
head(metadata, 3)
metadata$Group <- factor(metadata$Group, levels = c("IW", "CW", "TW"),
                         labels = c("Inland", "Coastal", "Tibetan"))
metadata$Saline <- factor(metadata$Saline)
stopifnot(!anyNA(metadata$Group), !anyNA(metadata$Saline))

knitr::kable(table(metadata$Group, metadata$Saline),
             caption = "湿地分组与盐度标签的交叉计数")
design <- model.matrix(~ Group + Saline, metadata)
c(design_columns = ncol(design), design_rank = qr(design)$rank)

min_prevalence <- 0.05
keep_feature <- rowSums(otutab > 0) >= ceiling(min_prevalence * ncol(otutab))
comm <- t(otutab[keep_feature, , drop = FALSE])
comm_rel <- sweep(comm, 1, rowSums(comm), "/")
stopifnot(identical(rownames(comm), rownames(metadata)),
          all(abs(rowSums(comm_rel) - 1) < 1e-10))
dist_bray <- vegan::vegdist(comm_rel, method = "bray")
c(samples = nrow(comm), retained_features = ncol(comm), original_features = nrow(otutab))

pcoa_raw <- stats::cmdscale(dist_bray, k = 2, eig = TRUE, add = FALSE)
eigen_tolerance <- max(abs(pcoa_raw$eig)) * 1e-8
negative_eigenvalues <- pcoa_raw$eig[pcoa_raw$eig < -eigen_tolerance]
negative_mass_ratio <- sum(abs(negative_eigenvalues)) /
  sum(pcoa_raw$eig[pcoa_raw$eig > eigen_tolerance])
needs_correction <- length(negative_eigenvalues) > 0L
pcoa_fit <- if (needs_correction) {
  stats::cmdscale(dist_bray, k = 2, eig = TRUE, add = TRUE)
} else pcoa_raw
correction_name <- if (needs_correction) "Cailliez" else "None"
correction_constant <- if (needs_correction) pcoa_fit$ac else 0
pcoa_variance <- 100 * pcoa_fit$eig[1:2] / sum(pcoa_fit$eig[pcoa_fit$eig > eigen_tolerance])
pcoa_diagnostic <- data.frame(
  NegativeEigenvalues = length(negative_eigenvalues),
  NegativeToPositiveMass = negative_mass_ratio,
  Correction = correction_name, Constant = correction_constant,
  Axis1Percent = pcoa_variance[1], Axis2Percent = pcoa_variance[2]
)
knitr::kable(pcoa_diagnostic, digits = 4)

pcoa_df <- data.frame(pcoa_fit$points, Group = metadata$Group)
names(pcoa_df)[1:2] <- c("Axis1", "Axis2")
p_pcoa <- ggplot(pcoa_df, aes(Axis1, Axis2, colour = Group)) +
  geom_point(size = 2.4, alpha = 0.85) + scale_color_pub(name = "Wetland") +
  coord_equal() + theme_pub() +
  labs(title = "Bray-Curtis PCoA",
       subtitle = sprintf("90 samples | %.1f%% in two dimensions", sum(pcoa_variance)),
       x = sprintf("PCoA 1 (%.1f%%)", pcoa_variance[1]),
       y = sprintf("PCoA 2 (%.1f%%)", pcoa_variance[2]))
save_pub(p_pcoa, "figures/22-pcoa-bray")

set.seed(20260718)
nmds_fit <- vegan::metaMDS(dist_bray, k = 2, trymax = 200,
                          autotransform = FALSE, trace = FALSE)
set.seed(20260718)
nmds_3d <- vegan::metaMDS(dist_bray, k = 3, trymax = 200,
                         autotransform = FALSE, trace = FALSE)
nmds_diagnostic <- data.frame(
  Dimensions = c(2L, 3L), Stress = c(nmds_fit$stress, nmds_3d$stress),
  RandomStarts = c(nmds_fit$tries, nmds_3d$tries),
  RepeatedBestSolutions = c(nmds_fit$converged, nmds_3d$converged)
)
knitr::kable(nmds_diagnostic, digits = 4)
nmds_sites <- vegan::scores(nmds_fit, display = "sites")
stopifnot(identical(rownames(nmds_sites), rownames(comm)))
nmds_df <- data.frame(nmds_sites, Group = metadata$Group)
p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, colour = Group)) +
  geom_point(size = 2.4, alpha = 0.85) + scale_color_pub(name = "Wetland") +
  coord_equal() + theme_pub() +
  labs(title = "Bray-Curtis NMDS",
       subtitle = sprintf("Two-dimensional stress = %.3f", nmds_fit$stress),
       x = "NMDS1", y = "NMDS2")
save_pub(p_nmds, "figures/22-nmds-bray")

distance_fit <- rbind(
  data.frame(Method = "PCoA: absolute distances", Original = as.vector(dist_bray),
             Projected = as.vector(dist(pcoa_fit$points))),
  data.frame(Method = "NMDS: distance order", Original = as.vector(dist_bray),
             Projected = as.vector(dist(nmds_sites)))
)
p_distance <- ggplot(distance_fit, aes(Original, Projected)) +
  geom_point(size = 0.65, alpha = 0.13, colour = "#0072B2") +
  facet_wrap(~ Method, scales = "free_y", nrow = 1) + theme_pub() +
  labs(x = "Original Bray-Curtis dissimilarity", y = "Distance in two dimensions")
if (!needs_correction) p_distance <- p_distance + geom_abline(
  data = data.frame(Method = "PCoA: absolute distances", Slope = 1, Intercept = 0),
  aes(slope = Slope, intercept = Intercept), linetype = 2, colour = "#D55E00")
if (!needs_correction) stopifnot(
  length(unique(ggplot_build(p_distance)$data[[2]]$PANEL)) == 1L)
save_pub(p_distance, "figures/22-distance-fit", width = 175, height = 90)
fit_rank_cor <- c(
  PCoA = cor(as.vector(dist_bray), as.vector(dist(pcoa_fit$points)), method = "spearman"),
  NMDS = cor(as.vector(dist_bray), as.vector(dist(nmds_sites)), method = "spearman"))
round(fit_rank_cor, 3)

clr_transform <- function(counts, pseudocount) {
  log_counts <- log(counts + pseudocount)
  log_counts - rowMeans(log_counts)
}
clr_matrix <- clr_transform(comm, pseudocount = 0.5)
stopifnot(max(abs(rowSums(clr_matrix))) < 1e-7)
clr_pca_fit <- prcomp(clr_matrix, center = TRUE, scale. = FALSE, rank. = 2)
clr_variance <- 100 * clr_pca_fit$sdev[1:2]^2 / sum(clr_pca_fit$sdev^2)
clr_df <- data.frame(clr_pca_fit$x[, 1:2], Group = metadata$Group)
p_clr <- ggplot(clr_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2.4, alpha = 0.85) + scale_color_pub(name = "Wetland") +
  coord_equal() + theme_pub() +
  labs(title = "CLR-PCA", subtitle = "Same features | 0.5 added to every count",
       x = sprintf("PC1 (%.1f%%)", clr_variance[1]),
       y = sprintf("PC2 (%.1f%%)", clr_variance[2]))
save_pub(p_clr, "figures/22-clr-pca")

keep_20 <- rowSums(otutab > 0) >= ceiling(0.20 * ncol(otutab))
comm_20 <- t(otutab[keep_20, , drop = FALSE])
dist_clr <- dist(clr_matrix)
rank_agreement <- function(a, b) cor(as.vector(a), as.vector(b), method = "spearman")
sensitivity <- data.frame(
  Comparison = c("CLR 0.1 vs 0.5; prevalence 5%", "CLR 1 vs 0.5; prevalence 5%",
                 "CLR prevalence 20% vs 5%; offset 0.5", "Bray prevalence 20% vs 5%",
                 "CLR 0.5 vs Bray; prevalence 5%"),
  Features = c(ncol(comm), ncol(comm), ncol(comm_20), ncol(comm_20), ncol(comm)),
  DistanceSpearman = c(
    rank_agreement(dist(clr_transform(comm, 0.1)), dist_clr),
    rank_agreement(dist(clr_transform(comm, 1)), dist_clr),
    rank_agreement(dist(clr_transform(comm_20, 0.5)), dist_clr),
    rank_agreement(vegdist(comm_20 / rowSums(comm_20), "bray"), dist_bray),
    rank_agreement(dist_clr, dist_bray))
)
knitr::kable(sensitivity, digits = 4)

test_one_group <- function(group, label) {
  model_data <- data.frame(group = group, row.names = rownames(metadata))
  set.seed(20260718)
  fit <- vegan::adonis2(dist_bray ~ group, data = model_data, permutations = 999)
  dispersion <- vegan::betadisper(dist_bray, group = group,
                                  type = "median", bias.adjust = TRUE)
  set.seed(20260718)
  test <- vegan::permutest(dispersion, permutations = 999)
  data.frame(Grouping = label, R2 = fit$R2[1], PERMANOVA_P = fit$`Pr(>F)`[1],
             Dispersion_P = test$tab$`Pr(>F)`[1])
}
group_tests <- rbind(test_one_group(metadata$Group, "Wetland group"),
                     test_one_group(metadata$Saline, "Soil salinity"))
knitr::kable(group_tests, digits = 3)

dir.create("results/22-ordination", recursive = TRUE, showWarnings = FALSE)
for (name in c("pcoa_diagnostic", "nmds_diagnostic", "sensitivity", "group_tests")) {
  write.table(get(name), paste0("results/22-ordination/", name, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
