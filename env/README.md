# Reproducible environments

## R and Quarto

- R: 4.4.1
- Quarto: 1.9.38
- Quarto Linux tarball SHA-256:
  `ea8c897368791ad9f200010c087ea3111b2e556b12a960487dd4e216902aa102`
- R lockfile: `renv.lock`

The executable pilot was rendered with knitr 1.51, rmarkdown 2.31, xfun
0.60, vegan 2.6-6.1, phyloseq 1.48.0, iNEXT 3.0.2, ggplot2 3.5.2,
ggrepel 0.9.5, colorspace 2.1.0, ragg 1.3.2, svglite 2.1.3, and
systemfonts 1.1.0.

## QIIME 2

`qiime2.yml` is the unmodified official Linux environment file for the QIIME
2 2026.4 distribution:

```text
https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.4/qiime2/released/rachis-qiime2-linux-64-conda.yml
```

Stored file SHA-256:

```text
68b0e9b974937525b71cca3dec48b7dc5b6135f588a941d1313c118c2a41f4b5
```

Create and test it with:

```bash
mamba env create -n microbiome-qiime2-2026.4 \
  -f env/qiime2.yml \
  --channel-priority flexible
mamba run -n microbiome-qiime2-2026.4 qiime info
```

The environment was created and validated locally on 2026-07-18. The installed
runtime reports Python 3.12.13, QIIME 2 2026.4.0, q2cli 2026.4.0, and 27
registered plugins. The command interfaces for paired-end DADA2 denoising,
scikit-learn taxonomy classification, and SEPP fragment insertion all returned
success. See `qiime2-validation.txt` for the bounded validation record.

Article 13 adds the fixed `qiime2-itsxpress-overlay.yml` to the same named
environment:

```bash
mamba env update -n microbiome-qiime2-2026.4 \
  -f env/qiime2-itsxpress-overlay.yml \
  --channel-priority flexible
conda run -n microbiome-qiime2-2026.4 qiime itsxpress trim-pair-output-unmerged \
  --help
```

The validated runtime reports ITSxpress 2.1.4, Biopython 1.87, VSEARCH 2.30.4,
and HMMER 3.4. It completed a fresh ITSxpress plus q2-dada2 run on 5,000 real
ITS2 read pairs; the sanitized command log and environment audit are under
`../results/13-its-18s/`.

Article 14 uses the same QIIME 2 environment without loading a pretrained
scikit-learn pickle. The validated runtime reports q2-feature-classifier
2026.4.0, RESCRIPt 2026.4.0, VSEARCH 2.30.4, and scikit-learn 1.7.1. It
completed three fresh `classify-consensus-vsearch` runs against region-matched
SILVA 138.2, GTDB R11-RS232, and Greengenes2 2024.09 references; all six
taxonomy/search Artifacts passed maximum validation and provenance checks.
The 91/91 audit and command log are under
`../results/14-taxonomy-databases/`.

Article 15 uses the same runtime to fit classifiers locally rather than
loading an ABI-unknown pickle. The validated contract fixes
q2-feature-classifier 2026.4.0, scikit-learn 1.7.1, joblib 1.5.3, NumPy
2.4.2, and SciPy 1.17.1. It trains a 236,516-record holdout model, evaluates
three confidence thresholds on 26,236 records, then refits the production
classifier on all 262,752 SILVA V4 references. The 80/80 audit and sanitized
command log are under `../results/15-train-classifier/`. Its sequence-search
comparator is recomputed within the same Article 15 validation step, so this
chapter does not depend on an Article 14 taxonomy Artifact.

Article 16 uses the same QIIME 2 runtime with q2-fragment-insertion 2026.4.0,
SEPP 4.5.6, HMMER 3.4, pplacer 1.1.alpha19, biom-format 2.1.17, and DendroPy
5.0.8. The official `sepp-refs-gg-13-8.qza` Artifact is checksum-locked to
`e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360`.
The fresh 366-ASV run produced a 203,818-tip rooted insertion tree, retained
all 5,213 aggregate reads, and passed 74/74 reference, Artifact, provenance,
jplace, tip, branch-length, and feature-retention checks. The sanitized
command log and audits are under `../results/16-sepp-phylogeny/`.

Article 17 uses R 4.4.1, Bioconductor 3.19, phyloseq 1.48.0, ape 5.8,
Biostrings 2.72.1, readr 2.1.5, ggplot2 3.5.2, digest 0.6.36, and jsonlite
1.8.8 from the locked R environment. Its validator reads a checksum-locked
366 × 4 count table, 366 × 7 taxonomy, 4 × 5 metadata, 366 representative
sequences, and the complete 203,818-tip rooted SEPP tree. It explicitly
prunes the tree to 366 ASV tips, constructs all five phyloseq components, and
passes 82/82 input, ID, orientation, tree, object, figure, and RDS round-trip
checks. The session record and audits are under
`../results/17-phyloseq-import/`.

Article 18 uses the same locked R 4.4.1 environment with ggplot2 3.5.2,
colorspace 2.1.0, ragg 1.3.2, svglite 2.1.3, and systemfonts 1.1.0. The
validated font contract resolves DejaVu Sans from the host font registry.
`cairo_pdf` and `svglite` produce the vector masters; `ragg::agg_png` and
`ragg::agg_tiff(compression = "lzw")` produce 600 ppi raster copies at the
declared final physical size. Four figures yield 8 vector and 8 raster files,
all 16 dimension checks pass, all four PDFs contain an embedded font, all four
SVGs retain text nodes, all four TIFFs report LZW compression, and the complete
validator passes 171/171 checks. The session record and audits are under
`../results/18-publication-graphics/`.

Article 19 adds iNEXT 3.0.2 to the same locked environment and uses vegan
2.6-6.1 plus phyloseq 1.48.0 for independent formula and implementation
parity checks. It computes unfiltered raw richness/diversity metrics,
analytical 10,000-read rarefaction, 90% coverage interpolation/extrapolation,
coverage sensitivity, and Hill q=0/1/2 on the 13,628-feature × 90-sample
wetland table. Four figures yield 8 vector and 8 600-ppi raster files; all
193 input, package, metric, sensitivity, format, font, and dimension checks
pass. The session record and audits are under
`../results/19-alpha-diversity/`.

The host mamba configuration uses strict channel priority, which cannot solve
the official file's `deblur 1.1.1` / `sortmerna 2.0` combination. Flexible
channel priority is therefore explicit in the reproducible creation command.
Real-data execution or bounded environment validation is complete for chapters
07–50.

## Native paired multi-omics and source tracking (Articles 46–50)

Articles 46–50 use the method projects' native interfaces. The modern Python
tools share `multiomics.yml`; MMvec is isolated because its released native
runtime requires Python 3.7 and TensorFlow 1.x. SourceTracker has a separate
environment so its dependency set cannot alter the multi-omics runtime.

Create and inspect the three environments with:

```bash
CONDA_CHANNEL_PRIORITY=flexible mamba env create -f env/mmvec-native.yml
conda run -n mmvec-native python -c \
  'import pkg_resources; print(pkg_resources.get_distribution("mmvec").version)'

mamba env create -f env/multiomics.yml
R_LIBS_USER="$PWD/.r-lib" conda run -n multiomics-native python -c \
  'from importlib.metadata import version; from halla import HAllA; from mofapy2.run.entry_point import entry_point; print(version("halla"), version("mofapy2"))'

mamba env create -f env/sourcetracker2.yml
conda run -n sourcetracker2 sourcetracker2 --help
```

The locked native versions are:

- HAllA 0.8.40 and mofapy2 0.7.4 under Python 3.10.19;
- MMvec 1.0.5 and TensorFlow 1.15.0 under Python 3.7.16;
- SourceTracker 2.0.1 under Python 3.10.20;
- mixOmics 6.26.0 and FEAST 0.1.0 in the 390-package R lock.

The bounded real-data commands are:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla scripts/run_article46_global_concordance.R
R_LIBS_USER="$PWD/.r-lib" conda run -n multiomics-native python scripts/run_article46_halla.py
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla scripts/run_article47_diablo.R
conda run -n mmvec-native python scripts/run_article48_native_mmvec.py --epochs 100
conda run -n multiomics-native python scripts/run_article48_mmvec_mofa.py
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla scripts/run_article49_multi_kingdom.R
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla scripts/run_article50_feast.R
conda run -n sourcetracker2 python scripts/run_article50_sourcetracker.py
```

Article 48 records the MMvec training objective and holdout mean absolute error
separately. Its selected epoch is the lowest observed holdout error among the
fixed 100 epochs; a boundary optimum is reported as such and is not called
early stopping or convergence. Native result ledgers are stored under
`../results/46-integration/` through `../results/50-source-tracking/`.
