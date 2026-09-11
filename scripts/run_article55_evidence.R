#!/usr/bin/env Rscript

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(svglite)
})
source("R/theme_pub.R")

root <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root, "data", "small", "causal-evidence")
out_dir <- file.path(root, "results", "55-causal-evidence")
figure_dir <- file.path(root, "figures", "55-causal-evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cases <- read.delim(
  file.path(data_dir, "evidence-cases.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(cases) == 7L,
  identical(cases$EvidenceLevel, 1:7),
  all(grepl("^10\\.", cases$DOI))
)

level_colours <- c(
  "1" = "#D9E2E8",
  "2" = "#BFD6DF",
  "3" = "#9AC6C2",
  "4" = "#78B89B",
  "5" = "#54A477",
  "6" = "#318D5F",
  "7" = "#166B45"
)
cases$LevelLabel <- paste0("L", cases$EvidenceLevel, "  ", cases$LevelName)
cases$LevelLabel <- factor(cases$LevelLabel, levels = rev(cases$LevelLabel))
cases$CaseLabel <- paste0(cases$FirstAuthor, " ", cases$Year)

figure_1 <- ggplot(cases, aes(EvidenceLevel, LevelLabel)) +
  geom_tile(aes(fill = factor(EvidenceLevel)), width = 0.88, height = 0.78, colour = "white") +
  geom_point(shape = 21, size = 4.2, fill = "white", colour = "#202020", stroke = 0.5) +
  geom_text(aes(label = CaseLabel), nudge_x = 0.52, hjust = 0, size = 3.2) +
  geom_segment(
    data = cases[-1, ],
    aes(
      x = EvidenceLevel - 1,
      xend = EvidenceLevel,
      y = LevelLabel,
      yend = LevelLabel
    ),
    inherit.aes = FALSE,
    colour = "#777777",
    linetype = 3
  ) +
  scale_fill_manual(values = level_colours) +
  scale_x_continuous(limits = c(0.5, 8.95), breaks = 1:7) +
  labs(
    title = "Causal language should rise only with design strength",
    subtitle = "Higher tiers reduce specific alternative explanations; none makes every bias disappear",
    x = "Evidence tier",
    y = NULL
  ) +
  theme_pub(10) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

threat_matrix <- expand.grid(
  LevelName = cases$LevelName,
  Threat = c(
    "Confounding",
    "Reverse causation",
    "Compositional bias",
    "Host-context transfer",
    "Intervention specificity",
    "Mechanism ambiguity"
  ),
  stringsAsFactors = FALSE
)
threat_values <- rbind(
  c("High", "High", "High", "Not tested", "Not tested", "High"),
  c("High", "Moderate", "High", "Not tested", "Not tested", "High"),
  c("Moderate", "Lower", "Moderate", "Not tested", "Not tested", "Moderate"),
  c("Lower", "Lower", "Moderate", "Not tested", "Moderate", "Moderate"),
  c("Lower", "Lower", "Moderate", "High", "Moderate", "Moderate"),
  c("Lower", "Lower", "Moderate", "Moderate", "Lower", "Lower"),
  c("Lower", "Lower", "Moderate", "Moderate", "Lower", "Lower")
)
threat_matrix$Risk <- as.vector(threat_values)
threat_matrix$LevelName <- factor(
  threat_matrix$LevelName,
  levels = rev(cases$LevelName)
)
threat_matrix$Threat <- factor(
  threat_matrix$Threat,
  levels = c(
    "Confounding", "Reverse causation", "Compositional bias",
    "Host-context transfer", "Intervention specificity", "Mechanism ambiguity"
  )
)
risk_colours <- c(
  High = "#C74B3F",
  Moderate = "#E4A13A",
  Lower = "#66A879",
  `Not tested` = "#B9BDC4"
)
figure_2 <- ggplot(threat_matrix, aes(Threat, LevelName, fill = Risk)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = Risk), size = 2.65, colour = "white", fontface = "bold") +
  scale_fill_manual(values = risk_colours) +
  labs(
    title = "Every design leaves a different threat profile",
    subtitle = "The matrix is a reading aid, not a numerical quality score",
    x = NULL,
    y = NULL,
    fill = "Residual threat"
  ) +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(angle = 27, hjust = 1),
    panel.grid = element_blank()
  )

claim_table <- data.frame(
  Design = c(
    "Cross-sectional",
    "Longitudinal",
    "Mediation / MR",
    "Human intervention",
    "Experimental transfer",
    "Defined rescue",
    "Mechanistic chain"
  ),
  Defensible = c(
    "is associated with",
    "precedes or predicts",
    "is consistent with an indirect or genetically proxied effect",
    "changed under the assigned intervention",
    "transferred a phenotype in the model",
    "restored a phenotype under controlled conditions",
    "acts through the specified pathway in the tested system"
  ),
  Overclaim = c(
    "causes",
    "drives disease",
    "proves mediation or causality",
    "the microbiome alone caused the change",
    "must cause the same phenotype in humans",
    "is universally protective",
    "is the only mechanism"
  ),
  stringsAsFactors = FALSE
)
claim_long <- rbind(
  data.frame(Design = claim_table$Design, Column = "Defensible wording", Text = claim_table$Defensible),
  data.frame(Design = claim_table$Design, Column = "Wording to avoid", Text = claim_table$Overclaim)
)
wrap_cell <- function(value, width = 36L) {
  vapply(
    value,
    function(item) paste(strwrap(item, width = width), collapse = "\n"),
    character(1)
  )
}
claim_long$Text <- wrap_cell(claim_long$Text)
claim_long$Design <- factor(claim_long$Design, levels = rev(claim_table$Design))
claim_long$Column <- factor(claim_long$Column, levels = c("Defensible wording", "Wording to avoid"))
figure_3 <- ggplot(claim_long, aes(Column, Design, fill = Column)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = Text), size = 2.25, lineheight = 0.9, colour = "#202020") +
  scale_fill_manual(values = c(
    "Defensible wording" = "#DCEFE2",
    "Wording to avoid" = "#F5D9D4"
  )) +
  labs(
    title = "Match the verb to what the design identifies",
    subtitle = "A bounded sentence is stronger than an unsupported causal headline",
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

triangulation_nodes <- data.frame(
  Node = c(
    "Human cohort", "Genetic proxy", "Human intervention",
    "Gnotobiotic transfer", "Defined isolate", "Molecular assay",
    "Bounded causal claim"
  ),
  X = c(0, 0, 0, 2.25, 2.25, 2.25, 4.6),
  Y = c(2.2, 1.2, 0.2, 2.2, 1.2, 0.2, 1.2),
  Group = c(rep("Human evidence", 3), rep("Experimental evidence", 3), "Claim")
)
triangulation_edges <- expand.grid(
  From = 1:6,
  To = 7
)
triangulation_edges <- merge(
  triangulation_edges,
  transform(triangulation_nodes, From = seq_len(nrow(triangulation_nodes))),
  by = "From"
)
triangulation_edges$XEnd <- triangulation_nodes$X[7] - 0.32
triangulation_edges$YEnd <- triangulation_nodes$Y[7]
figure_4 <- ggplot() +
  geom_segment(
    data = triangulation_edges,
    aes(X + 0.28, Y, xend = XEnd, yend = YEnd, colour = Group),
    linewidth = 0.75,
    alpha = 0.8,
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  ) +
  geom_label(
    data = triangulation_nodes,
    aes(X, Y, label = Node, fill = Group),
    size = 3.2,
    label.size = 0.35,
    label.padding = grid::unit(0.22, "lines")
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 2.62,
    label = "Different biases should point to the same bounded conclusion",
    size = 3.4,
    colour = "#444444"
  ) +
  annotate(
    "text",
    x = 4.6,
    y = 0.68,
    label = "Specific population, exposure,\noutcome and mechanism",
    size = 3,
    colour = "#555555"
  ) +
  scale_colour_manual(values = c(
    "Human evidence" = pal_pub[["blue"]],
    "Experimental evidence" = pal_pub[["orange"]]
  )) +
  scale_fill_manual(values = c(
    "Human evidence" = "#DFEDF6",
    "Experimental evidence" = "#F8E7C8",
    Claim = "#DCEFE2"
  )) +
  coord_cartesian(xlim = c(-0.6, 5.25), ylim = c(-0.2, 2.85), clip = "off") +
  labs(
    title = "Triangulation is convergence across independent bias structures",
    subtitle = "Repeating the same observational design does not create a new evidence tier"
  ) +
  theme_void(base_family = font_pub) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "#555555", size = 10),
    plot.margin = margin(10, 15, 10, 15)
  )

save_pub(figure_1, file.path(figure_dir, "55-1-evidence-ladder"), 180, 130)
save_pub(figure_2, file.path(figure_dir, "55-2-threat-matrix"), 180, 135)
save_pub(figure_3, file.path(figure_dir, "55-3-claim-calibration"), 180, 140)
save_pub(figure_4, file.path(figure_dir, "55-4-triangulation"), 180, 115)

write.table(
  threat_matrix,
  file.path(out_dir, "threat-matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  claim_table,
  file.path(out_dir, "claim-calibration.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
summary <- list(
  primary_cases = nrow(cases),
  evidence_tiers = nrow(cases),
  statistical_analysis = FALSE,
  framework = "Design capacity, residual threats, calibrated wording, and triangulation",
  interpretation = "The highest defensible causal claim is bounded by the weakest unresolved alternative explanation, not by the most visually impressive figure."
)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out_dir, "summary.json")
)

message("Article 55 evidence framework complete: ", nrow(cases), " primary-study cases.")
