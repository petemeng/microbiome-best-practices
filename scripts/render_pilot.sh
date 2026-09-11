#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: render_pilot.sh <project-root> <run-root>" >&2
  exit 2
fi

project_root="$(cd "$1" && pwd)"
run_root="$(cd "$2" && pwd)"
stage="${run_root}/pilot"

if [[ -n "${QUARTO_BIN:-}" ]]; then
  quarto_bin="${QUARTO_BIN}"
elif command -v quarto >/dev/null 2>&1; then
  quarto_bin="$(command -v quarto)"
elif [[ -x "${project_root}/.tools/quarto-1.9.38/bin/quarto" ]]; then
  quarto_bin="${project_root}/.tools/quarto-1.9.38/bin/quarto"
else
  echo "Quarto was not found; set QUARTO_BIN or install Quarto 1.9.38." >&2
  exit 127
fi

mkdir -p \
  "${stage}/.cache" \
  "${stage}/chapters" \
  "${stage}/data/small/decontam" \
  "${stage}/data/small/fastq" \
  "${stage}/data/small/primer-trimming" \
  "${stage}/data/small/its2" \
  "${stage}/data/small/its-18s" \
  "${stage}/data/small/taxonomy-databases" \
  "${stage}/data/small/phyloseq-import" \
  "${stage}/data/small/beta-distances" \
  "${stage}/data/small/permanova-dispersion" \
  "${stage}/data/small/community-typing-dmm" \
  "${stage}/data/small/absolute-quantification" \
  "${stage}/data/small/picrust2-chemerin" \
  "${stage}/data/small/cross-cohort-crc" \
  "${stage}/data/small/survival-t1d" \
  "${stage}/data/small/paired-ibd-multiomics" \
  "${stage}/data/small/multi-kingdom-duran" \
  "${stage}/data/small/source-tracking-feast" \
  "${stage}/data/small/longitudinal-dietswap" \
  "${stage}/data/small/mr-mibiogen" \
  "${stage}/data/small/causal-evidence" \
  "${stage}/figures" \
  "${stage}/scripts" \
  "${stage}/results/07-wsl2-conda" \
  "${stage}/results/08-qiime2-install" \
  "${stage}/results/09-r-ecosystem" \
  "${stage}/results/10-import-provenance-qc" \
  "${stage}/results/11-primer-trimming" \
  "${stage}/results/12-dada2-asv" \
  "${stage}/results/13-its-18s" \
  "${stage}/results/14-taxonomy-databases" \
  "${stage}/results/15-train-classifier" \
  "${stage}/results/16-sepp-phylogeny" \
  "${stage}/results/17-phyloseq-import" \
  "${stage}/results/18-publication-graphics" \
  "${stage}/results/19-alpha-diversity" \
  "${stage}/results/20-alpha-group-tests" \
  "${stage}/results/21-beta-distances" \
  "${stage}/results/23-permanova-dispersion" \
  "${stage}/results/25-environment-variance" \
  "${stage}/results/26-community-composition" \
  "${stage}/results/27-multirank-composition" \
  "${stage}/results/28-core-rare-biosphere" \
  "${stage}/results/29-community-typing-dmm" \
  "${stage}/results/36-network-robustness" \
  "${stage}/results/37-microbial-wgcna" \
  "${stage}/results/38-community-assembly-bnti" \
  "${stage}/results/39-neutral-community-model" \
  "${stage}/results/40-niche-distance-decay" \
  "${stage}/results/41-picrust2" \
  "${stage}/results/42-functional-guilds" \
  "${stage}/results/43-random-forest" \
  "${stage}/results/44-cross-cohort-validation" \
  "${stage}/results/45-survival-analysis" \
  "${stage}/results/46-integration" \
  "${stage}/results/47-spls-diablo" \
  "${stage}/results/48-mmvec-mofa" \
  "${stage}/results/49-multi-kingdom" \
  "${stage}/results/50-source-tracking" \
  "${stage}/results/51-longitudinal-analysis" \
  "${stage}/results/52-structural-equation-model" \
  "${stage}/results/53-mediation-analysis" \
  "${stage}/results/54-mendelian-randomization" \
  "${stage}/results/55-causal-evidence"
cp "${project_root}/qa/pilot/_quarto.yml" "${stage}/_quarto.yml"
cp "${project_root}/index.qmd" "${stage}/index.qmd"
cp "${project_root}/styles.scss" "${stage}/styles.scss"
cp "${project_root}/chapters/02-scope-and-limits.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/03-study-design.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/04-contamination-controls.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/05-batch-effects.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/06-data-and-fastq.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/07-wsl2-conda.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/08-qiime2-install.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/09-r-ecosystem.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/10-import-provenance-qc.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/11-primer-trimming.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/12-dada2-asv.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/13-its-18s.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/14-taxonomy-databases.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/15-train-classifier.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/16-sepp-phylogeny.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/17-phyloseq-import.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/18-publication-graphics.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/19-alpha-diversity.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/20-alpha-group-tests.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/21-beta-distances.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/22-ordination-unconstrained.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/23-permanova-dispersion.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/24-ordination-constrained-cap.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/25-environment-variance.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/26-community-composition.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/27-multirank-composition.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/28-core-rare-biosphere.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/29-community-typing-dmm.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/30-compositional-data.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/31-da-methods.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/32-multirank-da.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/33-da-visualization.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/34-absolute-quantification.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/35-cooccurrence-networks.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/36-network-robustness.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/37-microbial-wgcna.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/38-community-assembly-bnti.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/39-neutral-community-model.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/40-niche-distance-decay.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/41-picrust2.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/42-functional-guilds.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/43-random-forest.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/44-cross-cohort-validation.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/45-survival-analysis.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/46-procrustes-mantel-halla.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/47-spls-diablo.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/48-mmvec-mofa.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/49-multi-kingdom.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/50-source-tracking.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/51-longitudinal-analysis.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/52-structural-equation-model.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/53-mediation-analysis.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/54-mendelian-randomization.qmd" "${stage}/chapters/"
cp "${project_root}/chapters/55-causal-evidence.qmd" "${stage}/chapters/"
cp "${run_root}/data/small/"*.tsv "${stage}/data/small/"
cp "${project_root}/data/small/source_summary.json" "${stage}/data/small/"
cp "${project_root}/data/small/rooted-tree.nwk.gz" "${stage}/data/small/"
cp "${run_root}/data/small/decontam/"*.tsv "${stage}/data/small/decontam/"
cp "${project_root}/data/small/fastq/"* "${stage}/data/small/fastq/"
cp "${project_root}/data/small/primer-trimming/"* \
  "${stage}/data/small/primer-trimming/"
cp "${project_root}/data/small/its2/"* \
  "${stage}/data/small/its2/"
cp "${project_root}/data/small/its-18s/"* \
  "${stage}/data/small/its-18s/"
cp "${project_root}/data/small/taxonomy-databases/"* \
  "${stage}/data/small/taxonomy-databases/"
cp "${project_root}/data/small/phyloseq-import/"* \
  "${stage}/data/small/phyloseq-import/"
cp "${project_root}/data/small/beta-distances/"* \
  "${stage}/data/small/beta-distances/"
cp "${project_root}/data/small/permanova-dispersion/"* \
  "${stage}/data/small/permanova-dispersion/"
cp "${project_root}/data/small/community-typing-dmm/"* \
  "${stage}/data/small/community-typing-dmm/"
cp "${project_root}/data/small/absolute-quantification/"* \
  "${stage}/data/small/absolute-quantification/"
cp -R "${project_root}/data/small/picrust2-chemerin/." \
  "${stage}/data/small/picrust2-chemerin/"
cp -R "${project_root}/data/small/cross-cohort-crc/." \
  "${stage}/data/small/cross-cohort-crc/"
cp -R "${project_root}/data/small/survival-t1d/." \
  "${stage}/data/small/survival-t1d/"
cp -R "${project_root}/data/small/paired-ibd-multiomics/." \
  "${stage}/data/small/paired-ibd-multiomics/"
cp -R "${project_root}/data/small/multi-kingdom-duran/." \
  "${stage}/data/small/multi-kingdom-duran/"
cp -R "${project_root}/data/small/source-tracking-feast/." \
  "${stage}/data/small/source-tracking-feast/"
cp -R "${project_root}/data/small/longitudinal-dietswap/." \
  "${stage}/data/small/longitudinal-dietswap/"
cp -R "${project_root}/data/small/mr-mibiogen/." \
  "${stage}/data/small/mr-mibiogen/"
cp -R "${project_root}/data/small/causal-evidence/." \
  "${stage}/data/small/causal-evidence/"
cp "${project_root}/scripts/run_gemelli_rpca.py" "${stage}/scripts/"
cp -R "${project_root}/results/41-picrust2/." \
  "${stage}/results/41-picrust2/"
cp "${project_root}/figures/41-"* "${stage}/figures/"
cp "${run_root}/results/07-wsl2-conda/"* "${stage}/results/07-wsl2-conda/"
cp "${run_root}/figures/07-"* "${stage}/figures/"
cp "${run_root}/results/08-qiime2-install/"* "${stage}/results/08-qiime2-install/"
cp "${run_root}/figures/08-"* "${stage}/figures/"
cp "${run_root}/results/09-r-ecosystem/"* "${stage}/results/09-r-ecosystem/"
cp "${run_root}/figures/09-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/10-import-provenance-qc/." \
  "${stage}/results/10-import-provenance-qc/"
cp "${run_root}/figures/10-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/11-primer-trimming/." \
  "${stage}/results/11-primer-trimming/"
cp "${run_root}/figures/11-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/12-dada2-asv/." \
  "${stage}/results/12-dada2-asv/"
cp "${run_root}/figures/12-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/13-its-18s/." \
  "${stage}/results/13-its-18s/"
cp "${run_root}/figures/13-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/14-taxonomy-databases/." \
  "${stage}/results/14-taxonomy-databases/"
cp "${run_root}/figures/14-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/15-train-classifier/." \
  "${stage}/results/15-train-classifier/"
cp "${run_root}/figures/15-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/16-sepp-phylogeny/." \
  "${stage}/results/16-sepp-phylogeny/"
cp "${run_root}/figures/16-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/17-phyloseq-import/." \
  "${stage}/results/17-phyloseq-import/"
cp "${run_root}/figures/17-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/18-publication-graphics/." \
  "${stage}/results/18-publication-graphics/"
cp "${run_root}/figures/18-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/19-alpha-diversity/." \
  "${stage}/results/19-alpha-diversity/"
cp "${run_root}/figures/19-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/20-alpha-group-tests/." \
  "${stage}/results/20-alpha-group-tests/"
cp "${run_root}/figures/20-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/21-beta-distances/." \
  "${stage}/results/21-beta-distances/"
cp "${run_root}/figures/21-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/23-permanova-dispersion/." \
  "${stage}/results/23-permanova-dispersion/"
cp "${run_root}/figures/23-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/25-environment-variance/." \
  "${stage}/results/25-environment-variance/"
cp "${run_root}/figures/25-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/26-community-composition/." \
  "${stage}/results/26-community-composition/"
cp "${run_root}/figures/26-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/27-multirank-composition/." \
  "${stage}/results/27-multirank-composition/"
cp "${run_root}/figures/27-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/28-core-rare-biosphere/." \
  "${stage}/results/28-core-rare-biosphere/"
cp "${run_root}/figures/28-"* "${stage}/figures/"
cp -R \
  "${run_root}/results/29-community-typing-dmm/." \
  "${stage}/results/29-community-typing-dmm/"
cp "${run_root}/figures/29-"* "${stage}/figures/"
cp -R "${project_root}/results/46-integration/." \
  "${stage}/results/46-integration/"
cp -R "${project_root}/results/47-spls-diablo/." \
  "${stage}/results/47-spls-diablo/"
cp -R "${project_root}/results/48-mmvec-mofa/." \
  "${stage}/results/48-mmvec-mofa/"
cp -R "${project_root}/results/49-multi-kingdom/." \
  "${stage}/results/49-multi-kingdom/"
cp -R "${project_root}/results/50-source-tracking/." \
  "${stage}/results/50-source-tracking/"
cp -R "${project_root}/results/51-longitudinal-analysis/." \
  "${stage}/results/51-longitudinal-analysis/"
cp -R "${project_root}/results/52-structural-equation-model/." \
  "${stage}/results/52-structural-equation-model/"
cp -R "${project_root}/results/53-mediation-analysis/." \
  "${stage}/results/53-mediation-analysis/"
cp -R "${project_root}/results/54-mendelian-randomization/." \
  "${stage}/results/54-mendelian-randomization/"
cp -R "${project_root}/results/55-causal-evidence/." \
  "${stage}/results/55-causal-evidence/"
cp -R "${project_root}/figures/55-causal-evidence" \
  "${stage}/figures/"

(
  cd "${stage}"
  export XDG_CACHE_HOME="${stage}/.cache"
  if [[ -d "${project_root}/.r-lib" ]]; then
    export R_LIBS_USER="${project_root}/.r-lib"
  fi
  export GEMELLI_PYTHON="${project_root}/.tools/gemelli-0.0.13/bin/python"
  export GEMELLI_RUNNER="${stage}/scripts/run_gemelli_rpca.py"
  "${quarto_bin}" render
)

test -f "${stage}/_site/index.html"
test -f "${stage}/_site/chapters/02-scope-and-limits.html"
test -f "${stage}/_site/chapters/03-study-design.html"
test -f "${stage}/_site/chapters/04-contamination-controls.html"
test -f "${stage}/_site/chapters/05-batch-effects.html"
test -f "${stage}/_site/chapters/06-data-and-fastq.html"
test -f "${stage}/_site/chapters/07-wsl2-conda.html"
test -f "${stage}/_site/chapters/08-qiime2-install.html"
test -f "${stage}/_site/chapters/09-r-ecosystem.html"
test -f "${stage}/_site/chapters/10-import-provenance-qc.html"
test -f "${stage}/_site/chapters/11-primer-trimming.html"
test -f "${stage}/_site/chapters/12-dada2-asv.html"
test -f "${stage}/_site/chapters/13-its-18s.html"
test -f "${stage}/_site/chapters/14-taxonomy-databases.html"
test -f "${stage}/_site/chapters/15-train-classifier.html"
test -f "${stage}/_site/chapters/16-sepp-phylogeny.html"
test -f "${stage}/_site/chapters/17-phyloseq-import.html"
test -f "${stage}/_site/chapters/18-publication-graphics.html"
test -f "${stage}/_site/chapters/19-alpha-diversity.html"
test -f "${stage}/_site/chapters/20-alpha-group-tests.html"
test -f "${stage}/_site/chapters/21-beta-distances.html"
test -f "${stage}/_site/chapters/22-ordination-unconstrained.html"
test -f "${stage}/_site/chapters/23-permanova-dispersion.html"
test -f "${stage}/_site/chapters/24-ordination-constrained-cap.html"
test -f "${stage}/_site/chapters/25-environment-variance.html"
test -f "${stage}/_site/chapters/26-community-composition.html"
test -f "${stage}/_site/chapters/27-multirank-composition.html"
test -f "${stage}/_site/chapters/28-core-rare-biosphere.html"
test -f "${stage}/_site/chapters/29-community-typing-dmm.html"
test -f "${stage}/_site/chapters/30-compositional-data.html"
test -f "${stage}/_site/chapters/31-da-methods.html"
test -f "${stage}/_site/chapters/32-multirank-da.html"
test -f "${stage}/_site/chapters/33-da-visualization.html"
test -f "${stage}/_site/chapters/34-absolute-quantification.html"
test -f "${stage}/_site/chapters/35-cooccurrence-networks.html"
test -f "${stage}/_site/chapters/36-network-robustness.html"
test -f "${stage}/_site/chapters/37-microbial-wgcna.html"
test -f "${stage}/_site/chapters/38-community-assembly-bnti.html"
test -f "${stage}/_site/chapters/39-neutral-community-model.html"
test -f "${stage}/_site/chapters/40-niche-distance-decay.html"
test -f "${stage}/_site/chapters/41-picrust2.html"
test -f "${stage}/_site/chapters/42-functional-guilds.html"
test -f "${stage}/_site/chapters/43-random-forest.html"
test -f "${stage}/_site/chapters/44-cross-cohort-validation.html"
test -f "${stage}/_site/chapters/45-survival-analysis.html"
test -f "${stage}/_site/chapters/46-procrustes-mantel-halla.html"
test -f "${stage}/_site/chapters/47-spls-diablo.html"
test -f "${stage}/_site/chapters/48-mmvec-mofa.html"
test -f "${stage}/_site/chapters/49-multi-kingdom.html"
test -f "${stage}/_site/chapters/50-source-tracking.html"
test -f "${stage}/_site/chapters/51-longitudinal-analysis.html"
test -f "${stage}/_site/chapters/52-structural-equation-model.html"
test -f "${stage}/_site/chapters/53-mediation-analysis.html"
test -f "${stage}/_site/chapters/54-mendelian-randomization.html"
test -f "${stage}/_site/chapters/55-causal-evidence.html"

for stem in \
  30-closure-artifact \
  30-measurement-scales \
  30-reference-frame \
  31-da-hit-counts \
  31-da-jaccard \
  31-da-evidence-map \
  32-rank-evidence-cascade \
  32-multirank-effect-map \
  32-family-genus-coherence \
  33-da-volcano \
  33-da-cladogram \
  33-da-manhattan \
  33-da-forest \
  34-microbial-load \
  34-relative-quantitative-effects \
  34-qmp-exemplar \
  35-network-matrix \
  35-spiec-network \
  35-network-edge-audit \
  36-role-cartography \
  36-attack-robustness \
  36-group-rewiring \
  36-topology-sensitivity \
  37-soft-threshold \
  37-module-dendrogram \
  37-module-trait \
  37-hub-sensitivity \
  38-tree-filter \
  38-bnti-rcbray \
  38-process-fractions \
  38-phylogenetic-signal \
  39-neutral-abundance-occupancy \
  39-neutral-pool-migration \
  39-neutral-classification-sensitivity \
  39-neutral-detection-sensitivity \
  40-niche-breadth \
  40-rank-abundance-models \
  40-distance-decay \
  40-distance-decay-audit \
  41-picrust2-workflow \
  41-nsti-audit \
  41-pathway-heatmap \
  41-pathway-effect \
  42-function-coverage \
  42-functional-composition \
  42-functional-heatmap \
  42-functional-sensitivity \
  43-nested-roc \
  43-calibration \
  43-permutation-importance \
  43-regression-performance \
  44-cohort-pcoa \
  44-meta-forest \
  44-heterogeneity \
  44-external-validation \
  45-treatment-km \
  45-taxon-cox \
  45-ph-diagnostics \
  45-time-dependent-roc \
  46-integration/46-1-procrustes \
  46-integration/46-2-global-tests \
  46-integration/46-3-confounding-audit \
  46-integration/46-4-halla-associations \
  47-spls-diablo/47-1-spls-generalization \
  47-spls-diablo/47-2-diablo-scores \
  47-spls-diablo/47-3-test-performance \
  47-spls-diablo/47-4-feature-stability \
  48-mmvec-mofa/48-1-cca-overfit-audit \
  48-mmvec-mofa/48-2-mmvec-loss \
  48-mmvec-mofa/48-3-mmvec-top-pairs \
  48-mmvec-mofa/48-4-mofa-summary \
  49-multi-kingdom/49-1-three-kingdom-pcoa \
  49-multi-kingdom/49-2-procrustes-audit \
  49-multi-kingdom/49-3-shared-niche-audit \
  49-multi-kingdom/49-4-stable-candidates \
  50-source-tracking/50-1-feast-main \
  50-source-tracking/50-2-method-comparison \
  50-source-tracking/50-3-loo-calibration \
  50-source-tracking/50-4-missing-source \
  51-longitudinal-analysis/51-1-alpha-trajectories \
  51-longitudinal-analysis/51-2-pcoa-trajectories \
  51-longitudinal-analysis/51-3-volatility \
  51-longitudinal-analysis/51-4-model-audit \
  52-structural-equation-model/52-1-prespecified-dag \
  52-structural-equation-model/52-2-path-coefficients \
  52-structural-equation-model/52-3-model-diagnostics \
  52-structural-equation-model/52-4-model-audit \
  53-mediation-analysis/53-1-mediation-dag \
  53-mediation-analysis/53-2-effect-decomposition \
  53-mediation-analysis/53-3-unmeasured-confounding \
  53-mediation-analysis/53-4-observed-data \
  54-mendelian-randomization/54-1-harmonised-scatter \
  54-mendelian-randomization/54-2-method-forest \
  54-mendelian-randomization/54-3-leave-one-out \
  54-mendelian-randomization/54-4-assumption-audit \
  55-causal-evidence/55-1-evidence-ladder \
  55-causal-evidence/55-2-threat-matrix \
  55-causal-evidence/55-3-claim-calibration \
  55-causal-evidence/55-4-triangulation; do
  for extension in pdf svg png tiff; do
    test -f "${stage}/figures/${stem}.${extension}"
  done
done
test -f "${stage}/figures/01-preview-pcoa.pdf"
test -f "${stage}/figures/01-preview-pcoa.png"
test -f "${stage}/figures/01-preview-pcoa.tiff"
test -f "${stage}/figures/01-preview-cap.pdf"
test -f "${stage}/figures/01-preview-cap.png"
test -f "${stage}/figures/01-preview-cap.tiff"
test -f "${stage}/figures/01-learning-map.pdf"
test -f "${stage}/figures/01-learning-map.png"
test -f "${stage}/figures/01-learning-map.tiff"
test -f "${stage}/figures/02-taxonomy-resolution.pdf"
test -f "${stage}/figures/02-taxonomy-resolution.png"
test -f "${stage}/figures/02-taxonomy-resolution.tiff"
test -f "${stage}/figures/02-assay-scope-map.pdf"
test -f "${stage}/figures/02-assay-scope-map.png"
test -f "${stage}/figures/02-assay-scope-map.tiff"
test -f "${stage}/figures/03-design-confounding.pdf"
test -f "${stage}/figures/03-design-confounding.png"
test -f "${stage}/figures/03-design-confounding.tiff"
test -f "${stage}/figures/03-power-sensitivity.pdf"
test -f "${stage}/figures/03-power-sensitivity.png"
test -f "${stage}/figures/03-power-sensitivity.tiff"
test -f "${stage}/figures/04-control-library-size.pdf"
test -f "${stage}/figures/04-control-library-size.png"
test -f "${stage}/figures/04-control-library-size.tiff"
test -f "${stage}/figures/04-contaminant-prevalence.pdf"
test -f "${stage}/figures/04-contaminant-prevalence.png"
test -f "${stage}/figures/04-contaminant-prevalence.tiff"
test -f "${stage}/figures/04-contaminant-burden.pdf"
test -f "${stage}/figures/04-contaminant-burden.png"
test -f "${stage}/figures/04-contaminant-burden.tiff"
test -f "${stage}/results/04-contamination/contaminant-classification.tsv"
test -f "${stage}/results/04-contamination/sample-contamination-burden.tsv"
test -f "${stage}/results/04-contamination/otutab.filtered.tsv"
test -f "${stage}/results/04-contamination/taxonomy.filtered.tsv"
test -f "${stage}/results/04-contamination/metadata.filtered.tsv"
test -f "${stage}/figures/05-batch-design-audit.pdf"
test -f "${stage}/figures/05-batch-design-audit.png"
test -f "${stage}/figures/05-batch-design-audit.tiff"
test -f "${stage}/figures/05-community-batch-effects.pdf"
test -f "${stage}/figures/05-community-batch-effects.png"
test -f "${stage}/figures/05-community-batch-effects.tiff"
test -f "${stage}/figures/05-feature-batch-effects.pdf"
test -f "${stage}/figures/05-feature-batch-effects.png"
test -f "${stage}/figures/05-feature-batch-effects.tiff"
test -f "${stage}/figures/05-adjustment-diagnostic.pdf"
test -f "${stage}/figures/05-adjustment-diagnostic.png"
test -f "${stage}/figures/05-adjustment-diagnostic.tiff"
test -f "${stage}/results/05-batch/batch-design-counts.tsv"
test -f "${stage}/results/05-batch/community-batch-diagnostics.tsv"
test -f "${stage}/results/05-batch/feature-batch-diagnostics.tsv"
test -f "${stage}/results/05-batch/leave-one-plate-out.tsv"
test -f "${stage}/results/05-batch/adjustment-diagnostic.tsv"
test -f "${stage}/results/05-batch/plate-adjustment-coefficients.tsv"
test -f "${stage}/results/05-batch/ordination-scores.tsv"
test -f "${stage}/figures/06-read-architecture.pdf"
test -f "${stage}/figures/06-read-architecture.png"
test -f "${stage}/figures/06-read-architecture.tiff"
test -f "${stage}/figures/06-fastq-anatomy.pdf"
test -f "${stage}/figures/06-fastq-anatomy.png"
test -f "${stage}/figures/06-fastq-anatomy.tiff"
test -f "${stage}/figures/06-quality-profile.pdf"
test -f "${stage}/figures/06-quality-profile.png"
test -f "${stage}/figures/06-quality-profile.tiff"
test -f "${stage}/figures/06-file-contract.pdf"
test -f "${stage}/figures/06-file-contract.png"
test -f "${stage}/figures/06-file-contract.tiff"
test -f "${stage}/results/06-fastq/fastq-file-audit.tsv"
test -f "${stage}/results/06-fastq/quality-by-cycle.tsv"
test -f "${stage}/results/06-fastq/fastq-record-example.tsv"
test -f "${stage}/results/06-fastq/metadata-summary.tsv"
test -f "${stage}/figures/07-wsl2-layer-map.pdf"
test -f "${stage}/figures/07-wsl2-layer-map.png"
test -f "${stage}/figures/07-wsl2-layer-map.tiff"
test -f "${stage}/figures/07-environment-validation.pdf"
test -f "${stage}/figures/07-environment-validation.png"
test -f "${stage}/figures/07-environment-validation.tiff"
test -f "${stage}/results/07-wsl2-conda/environment-audit.tsv"
test -f "${stage}/results/07-wsl2-conda/entrypoint-audit.tsv"
test -f "${stage}/results/07-wsl2-conda/fastq-smoke.tsv"
test -f "${stage}/results/07-wsl2-conda/environment-summary.json"
test -f "${stage}/results/07-wsl2-conda/environment-validation.log"
test -f "${stage}/figures/08-qiime2-install-contract.pdf"
test -f "${stage}/figures/08-qiime2-install-contract.png"
test -f "${stage}/figures/08-qiime2-install-contract.tiff"
test -f "${stage}/figures/08-artifact-smoke-audit.pdf"
test -f "${stage}/figures/08-artifact-smoke-audit.png"
test -f "${stage}/figures/08-artifact-smoke-audit.tiff"
test -f "${stage}/figures/08-error-triage.pdf"
test -f "${stage}/figures/08-error-triage.png"
test -f "${stage}/figures/08-error-triage.tiff"
test -f "${stage}/results/08-qiime2-install/atacama-emp-paired-end.qza"
test -f "${stage}/results/08-qiime2-install/qiime-info.txt"
test -f "${stage}/results/08-qiime2-install/installation-audit.tsv"
test -f "${stage}/results/08-qiime2-install/artifact-peek.tsv"
test -f "${stage}/results/08-qiime2-install/artifact-content-audit.tsv"
test -f "${stage}/results/08-qiime2-install/artifact-summary.json"
test -f "${stage}/results/08-qiime2-install/error-probes.tsv"
test -f "${stage}/results/08-qiime2-install/import-validation.log"
test -f "${stage}/figures/09-r-ecosystem-map.pdf"
test -f "${stage}/figures/09-r-ecosystem-map.png"
test -f "${stage}/figures/09-r-ecosystem-map.tiff"
test -f "${stage}/figures/09-package-version-audit.pdf"
test -f "${stage}/figures/09-package-version-audit.png"
test -f "${stage}/figures/09-package-version-audit.tiff"
test -f "${stage}/figures/09-object-parity-audit.pdf"
test -f "${stage}/figures/09-object-parity-audit.png"
test -f "${stage}/figures/09-object-parity-audit.tiff"
test -f "${stage}/results/09-r-ecosystem/environment-summary.json"
test -f "${stage}/results/09-r-ecosystem/package-audit.tsv"
test -f "${stage}/results/09-r-ecosystem/lockfile-audit.tsv"
test -f "${stage}/results/09-r-ecosystem/input-audit.tsv"
test -f "${stage}/results/09-r-ecosystem/object-parity.tsv"
test -f "${stage}/results/09-r-ecosystem/sample-library-audit.tsv"
test -f "${stage}/results/09-r-ecosystem/r-session-info.txt"
test -f "${stage}/results/09-r-ecosystem/validation.log"
test -f "${stage}/figures/10-import-provenance-map.pdf"
test -f "${stage}/figures/10-import-provenance-map.png"
test -f "${stage}/figures/10-import-provenance-map.tiff"
test -f "${stage}/figures/10-fastqc-module-audit.pdf"
test -f "${stage}/figures/10-fastqc-module-audit.png"
test -f "${stage}/figures/10-fastqc-module-audit.tiff"
test -f "${stage}/figures/10-demultiplexing-qc.pdf"
test -f "${stage}/figures/10-demultiplexing-qc.png"
test -f "${stage}/figures/10-demultiplexing-qc.tiff"
test -f "${stage}/results/10-import-provenance-qc/atacama-emp-paired-end.qza"
test -f "${stage}/results/10-import-provenance-qc/atacama-demux-paired.qza"
test -f "${stage}/results/10-import-provenance-qc/atacama-demux-details.qza"
test -f "${stage}/results/10-import-provenance-qc/atacama-demux-summary.qzv"
test -f "${stage}/results/10-import-provenance-qc/fastqc/barcodes_fastqc.html"
test -f "${stage}/results/10-import-provenance-qc/fastqc/barcodes_fastqc.zip"
test -f "${stage}/results/10-import-provenance-qc/fastqc/forward_fastqc.html"
test -f "${stage}/results/10-import-provenance-qc/fastqc/forward_fastqc.zip"
test -f "${stage}/results/10-import-provenance-qc/fastqc/reverse_fastqc.html"
test -f "${stage}/results/10-import-provenance-qc/fastqc/reverse_fastqc.zip"
test -f "${stage}/results/10-import-provenance-qc/multiqc_report.html"
test -f "${stage}/results/10-import-provenance-qc/import-provenance-audit.tsv"
test -f "${stage}/results/10-import-provenance-qc/demux-sample-counts.tsv"
test -f "${stage}/results/10-import-provenance-qc/barcode-correction-audit.tsv"
test -f "${stage}/results/10-import-provenance-qc/fastqc-module-audit.tsv"
test -f "${stage}/results/10-import-provenance-qc/multiqc-general-stats.tsv"
test -f "${stage}/results/10-import-provenance-qc/qc-summary.json"
test -f "${stage}/results/10-import-provenance-qc/validation.log"
test -f "${stage}/figures/11-primer-presence-gate.pdf"
test -f "${stage}/figures/11-primer-presence-gate.png"
test -f "${stage}/figures/11-primer-presence-gate.tiff"
test -f "${stage}/figures/11-trimming-retention-audit.pdf"
test -f "${stage}/figures/11-trimming-retention-audit.png"
test -f "${stage}/figures/11-trimming-retention-audit.tiff"
test -f "${stage}/figures/11-primer-read-architecture.pdf"
test -f "${stage}/figures/11-primer-read-architecture.png"
test -f "${stage}/figures/11-primer-read-architecture.tiff"
test -f "${stage}/figures/11-region-branch-decision.pdf"
test -f "${stage}/figures/11-region-branch-decision.png"
test -f "${stage}/figures/11-region-branch-decision.tiff"
test -f "${stage}/results/11-primer-trimming/nfcore-v4-manifest.tsv"
test -f "${stage}/results/11-primer-trimming/nfcore-v4-demux.qza"
test -f "${stage}/results/11-primer-trimming/nfcore-v4-demux-summary.qzv"
test -f "${stage}/results/11-primer-trimming/nfcore-v4-trimmed.qza"
test -f "${stage}/results/11-primer-trimming/nfcore-v4-trimmed-summary.qzv"
test -f "${stage}/results/11-primer-trimming/primer-source-audit.tsv"
test -f "${stage}/results/11-primer-trimming/primer-match-audit.tsv"
test -f "${stage}/results/11-primer-trimming/trim-retention-audit.tsv"
test -f "${stage}/results/11-primer-trimming/length-shift-audit.tsv"
test -f "${stage}/results/11-primer-trimming/parameter-sensitivity.tsv"
test -f "${stage}/results/11-primer-trimming/region-branch-decision.tsv"
test -f "${stage}/results/11-primer-trimming/primer-provenance-audit.tsv"
test -f "${stage}/results/11-primer-trimming/atacama-no-primer-probe.log"
test -f "${stage}/results/11-primer-trimming/qiime-cutadapt.log"
test -f "${stage}/results/11-primer-trimming/primer-summary.json"
test -f "${stage}/results/11-primer-trimming/validation.log"
test -f "${stage}/figures/12-quality-truncation-gate.pdf"
test -f "${stage}/figures/12-quality-truncation-gate.png"
test -f "${stage}/figures/12-quality-truncation-gate.tiff"
test -f "${stage}/figures/12-denoising-loss-ledger.pdf"
test -f "${stage}/figures/12-denoising-loss-ledger.png"
test -f "${stage}/figures/12-denoising-loss-ledger.tiff"
test -f "${stage}/figures/12-overlap-budget.pdf"
test -f "${stage}/figures/12-overlap-budget.png"
test -f "${stage}/figures/12-overlap-budget.tiff"
test -f "${stage}/figures/12-asv-output-audit.pdf"
test -f "${stage}/figures/12-asv-output-audit.png"
test -f "${stage}/figures/12-asv-output-audit.tiff"
test -f "${stage}/results/12-dada2-asv/artifact-audit.tsv"
test -f "${stage}/results/12-dada2-asv/asv-length-distribution.tsv"
test -f "${stage}/results/12-dada2-asv/dada2-base-transitions.qza"
test -f "${stage}/results/12-dada2-asv/dada2-base-transitions.qzv"
test -f "${stage}/results/12-dada2-asv/dada2-feature-frequencies.qza"
test -f "${stage}/results/12-dada2-asv/dada2-rep-seqs.qza"
test -f "${stage}/results/12-dada2-asv/dada2-rep-seqs.qzv"
test -f "${stage}/results/12-dada2-asv/dada2-sample-frequencies.qza"
test -f "${stage}/results/12-dada2-asv/dada2-stats.qza"
test -f "${stage}/results/12-dada2-asv/dada2-stats.qzv"
test -f "${stage}/results/12-dada2-asv/dada2-summary.json"
test -f "${stage}/results/12-dada2-asv/dada2-table.qza"
test -f "${stage}/results/12-dada2-asv/dada2-table-summary.qzv"
test -f "${stage}/results/12-dada2-asv/environment-audit.tsv"
test -f "${stage}/results/12-dada2-asv/input-validation.log"
test -f "${stage}/results/12-dada2-asv/overlap-audit.tsv"
test -f "${stage}/results/12-dada2-asv/overtrim-expected-failure.log"
test -f "${stage}/results/12-dada2-asv/parameter-sensitivity.tsv"
test -f "${stage}/results/12-dada2-asv/provenance-audit.tsv"
test -f "${stage}/results/12-dada2-asv/qiime-dada2.log"
test -f "${stage}/results/12-dada2-asv/quality-by-cycle.tsv"
test -f "${stage}/results/12-dada2-asv/read-length-audit.tsv"
test -f "${stage}/results/12-dada2-asv/sample-retention.tsv"
test -f "${stage}/results/12-dada2-asv/validation.log"
test -f "${stage}/figures/13-marker-branch-map.pdf"
test -f "${stage}/figures/13-marker-branch-map.png"
test -f "${stage}/figures/13-marker-branch-map.tiff"
test -f "${stage}/figures/13-itsxpress-length-audit.pdf"
test -f "${stage}/figures/13-itsxpress-length-audit.png"
test -f "${stage}/figures/13-itsxpress-length-audit.tiff"
test -f "${stage}/figures/13-its2-loss-ledger.pdf"
test -f "${stage}/figures/13-its2-loss-ledger.png"
test -f "${stage}/figures/13-its2-loss-ledger.tiff"
test -f "${stage}/figures/13-reference-database-contract.pdf"
test -f "${stage}/figures/13-reference-database-contract.png"
test -f "${stage}/figures/13-reference-database-contract.tiff"
test -f "${stage}/results/13-its-18s/artifact-audit.tsv"
test -f "${stage}/results/13-its-18s/environment-audit.tsv"
test -f "${stage}/results/13-its-18s/its-18s-summary.json"
test -f "${stage}/results/13-its-18s/its2-asv-length-distribution.tsv"
test -f "${stage}/results/13-its-18s/its2-base-transitions.qza"
test -f "${stage}/results/13-its-18s/its2-denoising-ledger.tsv"
test -f "${stage}/results/13-its-18s/its2-denoising-stats.qza"
test -f "${stage}/results/13-its-18s/its2-fastq-audit.tsv"
test -f "${stage}/results/13-its-18s/its2-feature-frequencies.qza"
test -f "${stage}/results/13-its-18s/its2-length-audit.tsv"
test -f "${stage}/results/13-its-18s/its2-manifest.tsv"
test -f "${stage}/results/13-its-18s/its2-raw-summary.qzv"
test -f "${stage}/results/13-its-18s/its2-raw.qza"
test -f "${stage}/results/13-its-18s/its2-rep-seqs.qza"
test -f "${stage}/results/13-its-18s/its2-rep-seqs.qzv"
test -f "${stage}/results/13-its-18s/its2-retention.tsv"
test -f "${stage}/results/13-its-18s/its2-sample-frequencies.qza"
test -f "${stage}/results/13-its-18s/its2-stats.qzv"
test -f "${stage}/results/13-its-18s/its2-table-summary.qzv"
test -f "${stage}/results/13-its-18s/its2-table.qza"
test -f "${stage}/results/13-its-18s/its2-trimmed-summary.qzv"
test -f "${stage}/results/13-its-18s/its2-trimmed.qza"
test -f "${stage}/results/13-its-18s/marker-branch-decision.tsv"
test -f "${stage}/results/13-its-18s/pr2-domain-audit.tsv"
test -f "${stage}/results/13-its-18s/pr2-taxonomy-depth-audit.tsv"
test -f "${stage}/results/13-its-18s/provenance-audit.tsv"
test -f "${stage}/results/13-its-18s/qiime-itsxpress.log"
test -f "${stage}/results/13-its-18s/reference-database-audit.tsv"
test -f "${stage}/results/13-its-18s/validation-checks.tsv"
test -f "${stage}/results/13-its-18s/validation.log"
test -f "${stage}/figures/14-database-decision-map.pdf"
test -f "${stage}/figures/14-database-decision-map.png"
test -f "${stage}/figures/14-database-decision-map.tiff"
test -f "${stage}/figures/14-reference-scope-audit.pdf"
test -f "${stage}/figures/14-reference-scope-audit.png"
test -f "${stage}/figures/14-reference-scope-audit.tiff"
test -f "${stage}/figures/14-rank-assignment-coverage.pdf"
test -f "${stage}/figures/14-rank-assignment-coverage.png"
test -f "${stage}/figures/14-rank-assignment-coverage.tiff"
test -f "${stage}/figures/14-taxonomy-concordance.pdf"
test -f "${stage}/figures/14-taxonomy-concordance.png"
test -f "${stage}/figures/14-taxonomy-concordance.tiff"
test -f "${stage}/results/14-taxonomy-databases/taxonomy-databases-summary.json"
test -f "${stage}/results/14-taxonomy-databases/taxonomy-assignment-summary.tsv"
test -f "${stage}/results/14-taxonomy-databases/taxonomy-rank-coverage.tsv"
test -f "${stage}/results/14-taxonomy-databases/taxonomy-crosswalk.tsv"
test -f "${stage}/results/14-taxonomy-databases/organelle-offtarget-audit.tsv"
test -f "${stage}/results/14-taxonomy-databases/provenance-audit.tsv"
test -f "${stage}/results/14-taxonomy-databases/validation-checks.tsv"
test -f "${stage}/results/14-taxonomy-databases/validation.log"
test -f "${stage}/figures/15-classifier-build-contract.pdf"
test -f "${stage}/figures/15-classifier-build-contract.png"
test -f "${stage}/figures/15-classifier-build-contract.tiff"
test -f "${stage}/figures/15-holdout-rank-performance.pdf"
test -f "${stage}/figures/15-holdout-rank-performance.png"
test -f "${stage}/figures/15-holdout-rank-performance.tiff"
test -f "${stage}/figures/15-query-confidence-sensitivity.pdf"
test -f "${stage}/figures/15-query-confidence-sensitivity.png"
test -f "${stage}/figures/15-query-confidence-sensitivity.tiff"
test -f "${stage}/figures/15-method-concordance.pdf"
test -f "${stage}/figures/15-method-concordance.png"
test -f "${stage}/figures/15-method-concordance.tiff"
test -f "${stage}/results/15-train-classifier/classifier-training-summary.json"
test -f "${stage}/results/15-train-classifier/environment-audit.tsv"
test -f "${stage}/results/15-train-classifier/reference-split-audit.tsv"
test -f "${stage}/results/15-train-classifier/reference-split-membership.tsv.gz"
test -f "${stage}/results/15-train-classifier/holdout-rank-metrics.tsv"
test -f "${stage}/results/15-train-classifier/holdout-confidence-summary.tsv"
test -f "${stage}/results/15-train-classifier/query-rank-coverage.tsv"
test -f "${stage}/results/15-train-classifier/query-method-concordance.tsv"
test -f "${stage}/results/15-train-classifier/query-taxonomy-crosswalk.tsv"
test -f "${stage}/results/15-train-classifier/silva-v4-holdout-classifier.qza"
test -f "${stage}/results/15-train-classifier/silva-v4-classifier.qza"
test -f "${stage}/results/15-train-classifier/silva-v4-holdout-sequences.qza"
test -f "${stage}/results/15-train-classifier/silva-v4-holdout-expected-taxonomy.qza"
test -f "${stage}/results/15-train-classifier/holdout-taxonomy-c050.qza"
test -f "${stage}/results/15-train-classifier/holdout-taxonomy-c070.qza"
test -f "${stage}/results/15-train-classifier/holdout-taxonomy-c090.qza"
test -f "${stage}/results/15-train-classifier/query-taxonomy-c050.qza"
test -f "${stage}/results/15-train-classifier/query-taxonomy-c070.qza"
test -f "${stage}/results/15-train-classifier/query-taxonomy-c090.qza"
test -f "${stage}/results/15-train-classifier/holdout-evaluation-c070.qzv"
test -f "${stage}/results/15-train-classifier/provenance-audit.tsv"
test -f "${stage}/results/15-train-classifier/validation-checks.tsv"
test -f "${stage}/results/15-train-classifier/qiime-classifier.log"
test -f "${stage}/results/15-train-classifier/validation.log"
test -f "${stage}/figures/16-sepp-workflow-contract.pdf"
test -f "${stage}/figures/16-sepp-workflow-contract.png"
test -f "${stage}/figures/16-sepp-workflow-contract.tiff"
test -f "${stage}/figures/16-placement-retention-audit.pdf"
test -f "${stage}/figures/16-placement-retention-audit.png"
test -f "${stage}/figures/16-placement-retention-audit.tiff"
test -f "${stage}/figures/16-placement-diagnostics.pdf"
test -f "${stage}/figures/16-placement-diagnostics.png"
test -f "${stage}/figures/16-placement-diagnostics.tiff"
test -f "${stage}/figures/16-query-placement-tree.pdf"
test -f "${stage}/figures/16-query-placement-tree.png"
test -f "${stage}/figures/16-query-placement-tree.tiff"
test -f "${stage}/results/16-sepp-phylogeny/sepp-phylogeny-summary.json"
test -f "${stage}/results/16-sepp-phylogeny/environment-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/reference-database-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/placement-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/placement-support-bins.tsv"
test -f "${stage}/results/16-sepp-phylogeny/candidate-placement-counts.tsv"
test -f "${stage}/results/16-sepp-phylogeny/tip-reconciliation.tsv"
test -f "${stage}/results/16-sepp-phylogeny/feature-retention-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/artifact-provenance-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/command-execution-audit.tsv"
test -f "${stage}/results/16-sepp-phylogeny/validation-checks.tsv"
test -f "${stage}/results/16-sepp-phylogeny/query-v4-asvs.qza"
test -f "${stage}/results/16-sepp-phylogeny/query-abundance-table.qza"
test -f "${stage}/results/16-sepp-phylogeny/sepp-insertion-tree.qza"
test -f "${stage}/results/16-sepp-phylogeny/sepp-placements.qza"
test -f "${stage}/results/16-sepp-phylogeny/sepp-filtered-table.qza"
test -f "${stage}/results/16-sepp-phylogeny/sepp-removed-table.qza"
test -f "${stage}/results/16-sepp-phylogeny/sepp-insertion-tree.nwk"
test -f "${stage}/results/16-sepp-phylogeny/sepp-placements.jplace"
test -f "${stage}/results/16-sepp-phylogeny/query-only-placement-tree.nwk"
test -f "${stage}/results/16-sepp-phylogeny/sepp-filtered-table.tsv"
test -f "${stage}/results/16-sepp-phylogeny/sepp-removed-table.tsv"
test -f "${stage}/results/16-sepp-phylogeny/qiime-sepp.log"
test -f "${stage}/results/16-sepp-phylogeny/validation.log"
test -f "${stage}/figures/17-phyloseq-object-contract.pdf"
test -f "${stage}/figures/17-phyloseq-object-contract.png"
test -f "${stage}/figures/17-phyloseq-object-contract.tiff"
test -f "${stage}/figures/17-id-reconciliation-audit.pdf"
test -f "${stage}/figures/17-id-reconciliation-audit.png"
test -f "${stage}/figures/17-id-reconciliation-audit.tiff"
test -f "${stage}/figures/17-library-feature-audit.pdf"
test -f "${stage}/figures/17-library-feature-audit.png"
test -f "${stage}/figures/17-library-feature-audit.tiff"
test -f "${stage}/figures/17-taxonomy-coverage-audit.pdf"
test -f "${stage}/figures/17-taxonomy-coverage-audit.png"
test -f "${stage}/figures/17-taxonomy-coverage-audit.tiff"
test -f "${stage}/results/17-phyloseq-import/phyloseq-import-summary.json"
test -f "${stage}/results/17-phyloseq-import/input-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/package-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/id-reconciliation.tsv"
test -f "${stage}/results/17-phyloseq-import/component-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/sample-library-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/taxonomy-coverage.tsv"
test -f "${stage}/results/17-phyloseq-import/tree-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/sequence-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/rds-roundtrip-audit.tsv"
test -f "${stage}/results/17-phyloseq-import/validation-checks.tsv"
test -f "${stage}/results/17-phyloseq-import/atacama-v4-phyloseq.rds"
test -f "${stage}/results/17-phyloseq-import/r-session-info.txt"
test -f "${stage}/results/17-phyloseq-import/validation.log"
test -f "${stage}/figures/18-encoding-redundancy-audit.pdf"
test -f "${stage}/figures/18-encoding-redundancy-audit.svg"
test -f "${stage}/figures/18-encoding-redundancy-audit.png"
test -f "${stage}/figures/18-encoding-redundancy-audit.tiff"
test -f "${stage}/figures/18-final-size-community-summary.pdf"
test -f "${stage}/figures/18-final-size-community-summary.svg"
test -f "${stage}/figures/18-final-size-community-summary.png"
test -f "${stage}/figures/18-final-size-community-summary.tiff"
test -f "${stage}/figures/18-export-decision-map.pdf"
test -f "${stage}/figures/18-export-decision-map.svg"
test -f "${stage}/figures/18-export-decision-map.png"
test -f "${stage}/figures/18-export-decision-map.tiff"
test -f "${stage}/figures/18-raster-resolution-audit.pdf"
test -f "${stage}/figures/18-raster-resolution-audit.svg"
test -f "${stage}/figures/18-raster-resolution-audit.png"
test -f "${stage}/figures/18-raster-resolution-audit.tiff"
test -f "${stage}/results/18-publication-graphics/publication-graphics-summary.json"
test -f "${stage}/results/18-publication-graphics/input-audit.tsv"
test -f "${stage}/results/18-publication-graphics/package-audit.tsv"
test -f "${stage}/results/18-publication-graphics/command-audit.tsv"
test -f "${stage}/results/18-publication-graphics/font-audit.tsv"
test -f "${stage}/results/18-publication-graphics/palette-audit.tsv"
test -f "${stage}/results/18-publication-graphics/sample-plot-data.tsv"
test -f "${stage}/results/18-publication-graphics/phylum-composition.tsv"
test -f "${stage}/results/18-publication-graphics/figure-specification.tsv"
test -f "${stage}/results/18-publication-graphics/resolution-probe-audit.tsv"
test -f "${stage}/results/18-publication-graphics/format-audit.tsv"
test -f "${stage}/results/18-publication-graphics/validation-checks.tsv"
test -f "${stage}/results/18-publication-graphics/resolution-probe-72ppi.png"
test -f "${stage}/results/18-publication-graphics/resolution-probe-300ppi.png"
test -f "${stage}/results/18-publication-graphics/resolution-probe-600ppi.png"
test -f "${stage}/results/18-publication-graphics/r-session-info.txt"
test -f "${stage}/results/18-publication-graphics/validation.log"
test -f "${stage}/figures/19-depth-completeness-audit.pdf"
test -f "${stage}/figures/19-depth-completeness-audit.svg"
test -f "${stage}/figures/19-depth-completeness-audit.png"
test -f "${stage}/figures/19-depth-completeness-audit.tiff"
test -f "${stage}/figures/19-rarefaction-curves.pdf"
test -f "${stage}/figures/19-rarefaction-curves.svg"
test -f "${stage}/figures/19-rarefaction-curves.png"
test -f "${stage}/figures/19-rarefaction-curves.tiff"
test -f "${stage}/figures/19-standardization-comparison.pdf"
test -f "${stage}/figures/19-standardization-comparison.svg"
test -f "${stage}/figures/19-standardization-comparison.png"
test -f "${stage}/figures/19-standardization-comparison.tiff"
test -f "${stage}/figures/19-hill-profile.pdf"
test -f "${stage}/figures/19-hill-profile.svg"
test -f "${stage}/figures/19-hill-profile.png"
test -f "${stage}/figures/19-hill-profile.tiff"
test -f "${stage}/results/19-alpha-diversity/alpha-diversity-summary.json"
test -f "${stage}/results/19-alpha-diversity/input-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/package-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/command-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/font-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/library-depth-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/alpha-metrics-raw.tsv"
test -f "${stage}/results/19-alpha-diversity/alpha-metrics-size-10000.tsv"
test -f "${stage}/results/19-alpha-diversity/alpha-metrics-coverage-090.tsv"
test -f "${stage}/results/19-alpha-diversity/rarefaction-curve-data.tsv"
test -f "${stage}/results/19-alpha-diversity/coverage-sensitivity.tsv"
test -f "${stage}/results/19-alpha-diversity/standardization-rank-stability.tsv"
test -f "${stage}/results/19-alpha-diversity/group-descriptive-summary.tsv"
test -f "${stage}/results/19-alpha-diversity/metric-definition.tsv"
test -f "${stage}/results/19-alpha-diversity/figure-specification.tsv"
test -f "${stage}/results/19-alpha-diversity/format-audit.tsv"
test -f "${stage}/results/19-alpha-diversity/validation-checks.tsv"
test -f "${stage}/results/19-alpha-diversity/r-session-info.txt"
test -f "${stage}/results/19-alpha-diversity/validation.log"
test -f "${stage}/figures/20-design-selection-map.pdf"
test -f "${stage}/figures/20-design-selection-map.svg"
test -f "${stage}/figures/20-design-selection-map.png"
test -f "${stage}/figures/20-design-selection-map.tiff"
test -f "${stage}/figures/20-alpha-group-distributions.pdf"
test -f "${stage}/figures/20-alpha-group-distributions.svg"
test -f "${stage}/figures/20-alpha-group-distributions.png"
test -f "${stage}/figures/20-alpha-group-distributions.tiff"
test -f "${stage}/figures/20-standardization-sensitivity.pdf"
test -f "${stage}/figures/20-standardization-sensitivity.svg"
test -f "${stage}/figures/20-standardization-sensitivity.png"
test -f "${stage}/figures/20-standardization-sensitivity.tiff"
test -f "${stage}/figures/20-pairwise-effect-forest.pdf"
test -f "${stage}/figures/20-pairwise-effect-forest.svg"
test -f "${stage}/figures/20-pairwise-effect-forest.png"
test -f "${stage}/figures/20-pairwise-effect-forest.tiff"
test -f "${stage}/results/20-alpha-group-tests/alpha-group-tests-summary.json"
test -f "${stage}/results/20-alpha-group-tests/input-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/package-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/command-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/font-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/design-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/alpha-metrics-raw.tsv"
test -f "${stage}/results/20-alpha-group-tests/alpha-metrics-size-10000.tsv"
test -f "${stage}/results/20-alpha-group-tests/alpha-metrics-coverage-090.tsv"
test -f "${stage}/results/20-alpha-group-tests/group-descriptive-summary.tsv"
test -f "${stage}/results/20-alpha-group-tests/distribution-diagnostics.tsv"
test -f "${stage}/results/20-alpha-group-tests/omnibus-tests.tsv"
test -f "${stage}/results/20-alpha-group-tests/pairwise-tests.tsv"
test -f "${stage}/results/20-alpha-group-tests/welch-sensitivity.tsv"
test -f "${stage}/results/20-alpha-group-tests/standardization-sensitivity.tsv"
test -f "${stage}/results/20-alpha-group-tests/figure-specification.tsv"
test -f "${stage}/results/20-alpha-group-tests/format-audit.tsv"
test -f "${stage}/results/20-alpha-group-tests/validation-checks.tsv"
test -f "${stage}/results/20-alpha-group-tests/r-session-info.txt"
test -f "${stage}/results/20-alpha-group-tests/validation.log"
test -f "${stage}/figures/21-distance-choice-map.pdf"
test -f "${stage}/figures/21-distance-choice-map.svg"
test -f "${stage}/figures/21-distance-choice-map.png"
test -f "${stage}/figures/21-distance-choice-map.tiff"
test -f "${stage}/figures/21-distance-correlation.pdf"
test -f "${stage}/figures/21-distance-correlation.svg"
test -f "${stage}/figures/21-distance-correlation.png"
test -f "${stage}/figures/21-distance-correlation.tiff"
test -f "${stage}/figures/21-pcoa-metric-comparison.pdf"
test -f "${stage}/figures/21-pcoa-metric-comparison.svg"
test -f "${stage}/figures/21-pcoa-metric-comparison.png"
test -f "${stage}/figures/21-pcoa-metric-comparison.tiff"
test -f "${stage}/figures/21-aitchison-zero-sensitivity.pdf"
test -f "${stage}/figures/21-aitchison-zero-sensitivity.svg"
test -f "${stage}/figures/21-aitchison-zero-sensitivity.png"
test -f "${stage}/figures/21-aitchison-zero-sensitivity.tiff"
test -f "${stage}/results/21-beta-distances/beta-distances-summary.json"
test -f "${stage}/results/21-beta-distances/input-audit.tsv"
test -f "${stage}/results/21-beta-distances/package-audit.tsv"
test -f "${stage}/results/21-beta-distances/python-package-audit.tsv"
test -f "${stage}/results/21-beta-distances/command-audit.tsv"
test -f "${stage}/results/21-beta-distances/font-audit.tsv"
test -f "${stage}/results/21-beta-distances/tree-audit.tsv"
test -f "${stage}/results/21-beta-distances/rarefaction-audit.tsv"
test -f "${stage}/results/21-beta-distances/feature-filter-audit.tsv"
test -f "${stage}/results/21-beta-distances/distance-summary.tsv"
test -f "${stage}/results/21-beta-distances/distance-correlation.tsv"
test -f "${stage}/results/21-beta-distances/pcoa-eigen-audit.tsv"
test -f "${stage}/results/21-beta-distances/pcoa-scores.tsv"
test -f "${stage}/results/21-beta-distances/pairwise-distances.tsv"
test -f "${stage}/results/21-beta-distances/aitchison-pseudocount-sensitivity.tsv"
test -f "${stage}/results/21-beta-distances/robust-aitchison-distance.tsv"
test -f "${stage}/results/21-beta-distances/robust-aitchison-sample-scores.tsv"
test -f "${stage}/results/21-beta-distances/robust-aitchison-feature-loadings.tsv"
test -f "${stage}/results/21-beta-distances/robust-aitchison-eigenvalues.tsv"
test -f "${stage}/results/21-beta-distances/robust-aitchison-ordination.txt"
test -f "${stage}/results/21-beta-distances/gemelli-summary.json"
test -f "${stage}/results/21-beta-distances/gemelli-determinism-audit.tsv"
test -f "${stage}/results/21-beta-distances/figure-specification.tsv"
test -f "${stage}/results/21-beta-distances/format-audit.tsv"
test -f "${stage}/results/21-beta-distances/validation-checks.tsv"
test -f "${stage}/results/21-beta-distances/r-session-info.txt"
test -f "${stage}/results/21-beta-distances/validation.log"
test -f "${stage}/figures/22-pcoa-bray.pdf"
test -f "${stage}/figures/22-pcoa-bray.png"
test -f "${stage}/figures/22-pcoa-bray.tiff"
test -f "${stage}/figures/23-exchangeability-design.pdf"
test -f "${stage}/figures/23-exchangeability-design.svg"
test -f "${stage}/figures/23-exchangeability-design.png"
test -f "${stage}/figures/23-exchangeability-design.tiff"
test -f "${stage}/figures/23-pcoa-centroid-dispersion.pdf"
test -f "${stage}/figures/23-pcoa-centroid-dispersion.svg"
test -f "${stage}/figures/23-pcoa-centroid-dispersion.png"
test -f "${stage}/figures/23-pcoa-centroid-dispersion.tiff"
test -f "${stage}/figures/23-pairwise-permanova.pdf"
test -f "${stage}/figures/23-pairwise-permanova.svg"
test -f "${stage}/figures/23-pairwise-permanova.png"
test -f "${stage}/figures/23-pairwise-permanova.tiff"
test -f "${stage}/figures/23-dispersion-by-treatment.pdf"
test -f "${stage}/figures/23-dispersion-by-treatment.svg"
test -f "${stage}/figures/23-dispersion-by-treatment.png"
test -f "${stage}/figures/23-dispersion-by-treatment.tiff"
test -f "${stage}/results/23-permanova-dispersion/permanova-dispersion-summary.json"
test -f "${stage}/results/23-permanova-dispersion/design-audit.tsv"
test -f "${stage}/results/23-permanova-dispersion/aggregation-audit.tsv"
test -f "${stage}/results/23-permanova-dispersion/feature-filter-audit.tsv"
test -f "${stage}/results/23-permanova-dispersion/permutation-audit.tsv"
test -f "${stage}/results/23-permanova-dispersion/permanova-omnibus.tsv"
test -f "${stage}/results/23-permanova-dispersion/permanova-factorial.tsv"
test -f "${stage}/results/23-permanova-dispersion/pairwise-permanova.tsv"
test -f "${stage}/results/23-permanova-dispersion/dispersion-global.tsv"
test -f "${stage}/results/23-permanova-dispersion/sensitivity-analysis.tsv"
test -f "${stage}/results/23-permanova-dispersion/validation-checks.tsv"
test -f "${stage}/results/23-permanova-dispersion/validation.log"
test -f "${stage}/figures/24-cap-environment.pdf"
test -f "${stage}/figures/24-cap-environment.png"
test -f "${stage}/figures/24-cap-environment.tiff"
test -f "${stage}/figures/25-environment-correlation.pdf"
test -f "${stage}/figures/25-environment-correlation.svg"
test -f "${stage}/figures/25-environment-correlation.png"
test -f "${stage}/figures/25-environment-correlation.tiff"
test -f "${stage}/figures/25-envfit-pcoa.pdf"
test -f "${stage}/figures/25-envfit-pcoa.svg"
test -f "${stage}/figures/25-envfit-pcoa.png"
test -f "${stage}/figures/25-envfit-pcoa.tiff"
test -f "${stage}/figures/25-mantel-distance-blocks.pdf"
test -f "${stage}/figures/25-mantel-distance-blocks.svg"
test -f "${stage}/figures/25-mantel-distance-blocks.png"
test -f "${stage}/figures/25-mantel-distance-blocks.tiff"
test -f "${stage}/figures/25-variation-partition.pdf"
test -f "${stage}/figures/25-variation-partition.svg"
test -f "${stage}/figures/25-variation-partition.png"
test -f "${stage}/figures/25-variation-partition.tiff"
test -f "${stage}/results/25-environment-variance/environment-variance-summary.json"
test -f "${stage}/results/25-environment-variance/environmental-correlation.tsv"
test -f "${stage}/results/25-environment-variance/envfit-results.tsv"
test -f "${stage}/results/25-environment-variance/mantel-results.tsv"
test -f "${stage}/results/25-environment-variance/partial-mantel-sensitivity.tsv"
test -f "${stage}/results/25-environment-variance/variation-partition.tsv"
test -f "${stage}/results/25-environment-variance/variation-fraction-tests.tsv"
test -f "${stage}/results/25-environment-variance/prevalence-sensitivity.tsv"
test -f "${stage}/results/25-environment-variance/conditioning-sensitivity.tsv"
test -f "${stage}/results/25-environment-variance/validation-checks.tsv"
test -f "${stage}/results/25-environment-variance/validation.log"
test -f "${stage}/figures/26-phylum-stacked.pdf"
test -f "${stage}/figures/26-phylum-stacked.svg"
test -f "${stage}/figures/26-phylum-stacked.png"
test -f "${stage}/figures/26-phylum-stacked.tiff"
test -f "${stage}/figures/26-genus-bubble.pdf"
test -f "${stage}/figures/26-genus-bubble.svg"
test -f "${stage}/figures/26-genus-bubble.png"
test -f "${stage}/figures/26-genus-bubble.tiff"
test -f "${stage}/figures/26-group-phylum-alluvial.pdf"
test -f "${stage}/figures/26-group-phylum-alluvial.svg"
test -f "${stage}/figures/26-group-phylum-alluvial.png"
test -f "${stage}/figures/26-group-phylum-alluvial.tiff"
test -f "${stage}/figures/26-genus-heatmap.pdf"
test -f "${stage}/figures/26-genus-heatmap.svg"
test -f "${stage}/figures/26-genus-heatmap.png"
test -f "${stage}/figures/26-genus-heatmap.tiff"
test -f "${stage}/results/26-community-composition/community-composition-summary.json"
test -f "${stage}/results/26-community-composition/taxonomy-rank-audit.tsv"
test -f "${stage}/results/26-community-composition/phylum-display-composition.tsv"
test -f "${stage}/results/26-community-composition/phylum-group-summary.tsv"
test -f "${stage}/results/26-community-composition/genus-bubble-summary.tsv"
test -f "${stage}/results/26-community-composition/top-n-sensitivity.tsv"
test -f "${stage}/results/26-community-composition/denominator-sensitivity.tsv"
test -f "${stage}/results/26-community-composition/aggregation-sensitivity.tsv"
test -f "${stage}/results/26-community-composition/validation-checks.tsv"
test -f "${stage}/results/26-community-composition/validation.log"
test -f "${stage}/figures/27-rank-resolution-cascade.pdf"
test -f "${stage}/figures/27-rank-resolution-cascade.svg"
test -f "${stage}/figures/27-rank-resolution-cascade.png"
test -f "${stage}/figures/27-rank-resolution-cascade.tiff"
test -f "${stage}/figures/27-multirank-bubble.pdf"
test -f "${stage}/figures/27-multirank-bubble.svg"
test -f "${stage}/figures/27-multirank-bubble.png"
test -f "${stage}/figures/27-multirank-bubble.tiff"
test -f "${stage}/figures/27-lineage-ladder.pdf"
test -f "${stage}/figures/27-lineage-ladder.svg"
test -f "${stage}/figures/27-lineage-ladder.png"
test -f "${stage}/figures/27-lineage-ladder.tiff"
test -f "${stage}/figures/27-topn-coverage.pdf"
test -f "${stage}/figures/27-topn-coverage.svg"
test -f "${stage}/figures/27-topn-coverage.png"
test -f "${stage}/figures/27-topn-coverage.tiff"
test -f "${stage}/results/27-multirank-composition/multirank-composition-summary.json"
test -f "${stage}/results/27-multirank-composition/taxonomy-lineage-map.tsv"
test -f "${stage}/results/27-multirank-composition/multirank-relative-abundance.tsv"
test -f "${stage}/results/27-multirank-composition/rank-resolution-audit.tsv"
test -f "${stage}/results/27-multirank-composition/group-rank-resolution.tsv"
test -f "${stage}/results/27-multirank-composition/lineage-gap-audit.tsv"
test -f "${stage}/results/27-multirank-composition/label-collision-audit.tsv"
test -f "${stage}/results/27-multirank-composition/rank-topn-sensitivity.tsv"
test -f "${stage}/results/27-multirank-composition/complete-lineage-abundance.tsv"
test -f "${stage}/results/27-multirank-composition/denominator-sensitivity.tsv"
test -f "${stage}/results/27-multirank-composition/validation-checks.tsv"
test -f "${stage}/results/27-multirank-composition/validation.log"
test -f "${stage}/figures/28-occupancy-abundance.pdf"
test -f "${stage}/figures/28-occupancy-abundance.svg"
test -f "${stage}/figures/28-occupancy-abundance.png"
test -f "${stage}/figures/28-occupancy-abundance.tiff"
test -f "${stage}/figures/28-group-core-membership.pdf"
test -f "${stage}/figures/28-group-core-membership.svg"
test -f "${stage}/figures/28-group-core-membership.png"
test -f "${stage}/figures/28-group-core-membership.tiff"
test -f "${stage}/figures/28-threshold-depth-sensitivity.pdf"
test -f "${stage}/figures/28-threshold-depth-sensitivity.svg"
test -f "${stage}/figures/28-threshold-depth-sensitivity.png"
test -f "${stage}/figures/28-threshold-depth-sensitivity.tiff"
test -f "${stage}/figures/28-rare-biosphere-mass.pdf"
test -f "${stage}/figures/28-rare-biosphere-mass.svg"
test -f "${stage}/figures/28-rare-biosphere-mass.png"
test -f "${stage}/figures/28-rare-biosphere-mass.tiff"
test -f "${stage}/results/28-core-rare-biosphere/core-rare-biosphere-summary.json"
test -f "${stage}/results/28-core-rare-biosphere/feature-occupancy-abundance.tsv"
test -f "${stage}/results/28-core-rare-biosphere/core-members.tsv"
test -f "${stage}/results/28-core-rare-biosphere/group-core-summary.tsv"
test -f "${stage}/results/28-core-rare-biosphere/core-membership-patterns.tsv"
test -f "${stage}/results/28-core-rare-biosphere/core-threshold-sensitivity.tsv"
test -f "${stage}/results/28-core-rare-biosphere/depth-sensitivity.tsv"
test -f "${stage}/results/28-core-rare-biosphere/feature-state-classification.tsv"
test -f "${stage}/results/28-core-rare-biosphere/rare-mass-by-sample.tsv"
test -f "${stage}/results/28-core-rare-biosphere/rare-threshold-sensitivity.tsv"
test -f "${stage}/results/28-core-rare-biosphere/validation-checks.tsv"
test -f "${stage}/results/28-core-rare-biosphere/validation.log"
test -f "${stage}/figures/29-model-selection.pdf"
test -f "${stage}/figures/29-model-selection.svg"
test -f "${stage}/figures/29-model-selection.png"
test -f "${stage}/figures/29-model-selection.tiff"
test -f "${stage}/figures/29-posterior-ordination.pdf"
test -f "${stage}/figures/29-posterior-ordination.svg"
test -f "${stage}/figures/29-posterior-ordination.png"
test -f "${stage}/figures/29-posterior-ordination.tiff"
test -f "${stage}/figures/29-component-profiles.pdf"
test -f "${stage}/figures/29-component-profiles.svg"
test -f "${stage}/figures/29-component-profiles.png"
test -f "${stage}/figures/29-component-profiles.tiff"
test -f "${stage}/figures/29-stability-audit.pdf"
test -f "${stage}/figures/29-stability-audit.svg"
test -f "${stage}/figures/29-stability-audit.png"
test -f "${stage}/figures/29-stability-audit.tiff"
test -f "${stage}/results/29-community-typing-dmm/community-typing-dmm-summary.json"
test -f "${stage}/results/29-community-typing-dmm/model-selection.tsv"
test -f "${stage}/results/29-community-typing-dmm/sample-posteriors.tsv"
test -f "${stage}/results/29-community-typing-dmm/component-profiles.tsv"
test -f "${stage}/results/29-community-typing-dmm/initialization-agreement.tsv"
test -f "${stage}/results/29-community-typing-dmm/sensitivity-agreement.tsv"
test -f "${stage}/results/29-community-typing-dmm/participant-transitions.tsv"
test -f "${stage}/results/29-community-typing-dmm/validation-checks.tsv"
test -f "${stage}/results/29-community-typing-dmm/validation.log"
