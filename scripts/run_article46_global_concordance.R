#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(permute)
  library(vegan)
})

set.seed(20260726)
root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "paired-ibd-multiomics")
out_dir <- file.path(root, "results", "46-integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

otutab <- as.matrix(read.delim(
  file.path(data_dir, "otutab.tsv"), row.names = 1, check.names = FALSE
))
metabolites <- as.matrix(read.delim(
  file.path(data_dir, "metabolites.tsv"), row.names = 1, check.names = FALSE
))
metadata <- read.delim(
  file.path(data_dir, "metadata.tsv"), row.names = 1, check.names = FALSE
)
storage.mode(otutab) <- "numeric"
storage.mode(metabolites) <- "numeric"
stopifnot(
  identical(colnames(otutab), rownames(metadata)),
  identical(colnames(metabolites), rownames(metadata))
)

microbe_log <- log(t(otutab) + 0.5)
microbe_clr <- microbe_log - rowMeans(microbe_log)
metabolite_log <- log1p(t(metabolites))
metabolite_scaled <- scale(metabolite_log)
metabolite_scaled[, !is.finite(colSums(metabolite_scaled))] <- 0

design <- model.matrix(~ StudyGroup, data = metadata)
residualize <- function(value) stats::lm.fit(design, value)$residuals
views <- list(
  Observed = list(Microbiome = microbe_clr, Metabolome = metabolite_scaled),
  `Group-adjusted` = list(
    Microbiome = residualize(microbe_clr),
    Metabolome = residualize(metabolite_scaled)
  )
)

test_rows <- list()
arrow_rows <- list()
for (model_index in seq_along(views)) {
  model_name <- names(views)[model_index]
  x <- views[[model_name]]$Microbiome
  y <- views[[model_name]]$Metabolome
  distance_x <- stats::dist(x, method = "euclidean")
  distance_y <- stats::dist(y, method = "euclidean")
  ordination_x <- stats::cmdscale(distance_x, k = 5, eig = TRUE)$points
  ordination_y <- stats::cmdscale(distance_y, k = 5, eig = TRUE)$points
  control <- if (model_name == "Observed") {
    permute::how(nperm = 999)
  } else {
    permute::how(nperm = 999, blocks = metadata$StudyGroup)
  }
  set.seed(20260726 + model_index)
  procrustes_test <- vegan::protest(
    ordination_x,
    ordination_y,
    permutations = control,
    symmetric = TRUE
  )
  set.seed(20260736 + model_index)
  mantel_test <- vegan::mantel(
    distance_x,
    distance_y,
    method = "spearman",
    permutations = control
  )
  test_rows[[length(test_rows) + 1L]] <- data.frame(
    Model = model_name,
    Method = "Procrustes",
    Statistic = unname(procrustes_test$t0),
    PValue = unname(procrustes_test$signif)
  )
  test_rows[[length(test_rows) + 1L]] <- data.frame(
    Model = model_name,
    Method = "Mantel",
    Statistic = unname(mantel_test$statistic),
    PValue = unname(mantel_test$signif)
  )
  fit <- vegan::procrustes(ordination_x, ordination_y, symmetric = TRUE)
  arrow_rows[[model_name]] <- data.frame(
    SampleID = rownames(metadata),
    StudyGroup = metadata$StudyGroup,
    Model = model_name,
    StartAxis1 = fit$X[, 1],
    StartAxis2 = fit$X[, 2],
    EndAxis1 = fit$Yrot[, 1],
    EndAxis2 = fit$Yrot[, 2]
  )
}

tests <- do.call(rbind, test_rows)
tests$QValue <- stats::p.adjust(tests$PValue, method = "BH")
arrows <- do.call(rbind, arrow_rows)
write.table(
  tests,
  file.path(out_dir, "global-concordance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  arrows,
  file.path(out_dir, "procrustes-arrows.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary <- list(
  status = "passed",
  seed = 20260726,
  samples = nrow(metadata),
  microbe_features = ncol(microbe_clr),
  metabolite_features = ncol(metabolite_scaled),
  tests = split(tests, seq_len(nrow(tests))),
  versions = list(
    R = as.character(getRversion()),
    vegan = as.character(packageVersion("vegan")),
    permute = as.character(packageVersion("permute"))
  )
)
write_json(
  summary,
  file.path(out_dir, "summary-global.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
