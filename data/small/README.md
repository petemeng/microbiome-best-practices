# Real 16S tables used by the formal articles

## microeco wetland tables

These four TSV files are generated from the example data distributed in the
CRAN source archive for `microeco 2.0.0`:

- Source:
  `https://cran.r-project.org/src/contrib/Archive/microeco/microeco_2.0.0.tar.gz`
- Source SHA-256:
  `454a3b71ceeea86bdd54f475a12b5bac9c3b124f663b8dc6e560827fca89ebfd`
- Package license: GPL-3
- Data study: An J, Liu C, Wang Q, et al. *Geoderma*. 2019;337:290–299.
  <https://doi.org/10.1016/j.geoderma.2018.09.035>

Regenerate the files with:

```bash
Rscript scripts/prepare_pilot_data.R \
  downloads/microeco_2.0.0.tar.gz \
  data/small \
  data/small/source_summary.json
```

`otutab.tsv` is features × samples. `taxonomy.tsv` is features × taxonomic
ranks. `metadata.tsv` and `environment.tsv` are samples × variables. Article
09 reads the same triad independently and verifies that its feature, sample,
total-read, and per-sample library-size contracts are identical in
`phyloseq 1.48.0` and `microeco 2.0.0`. Article 25 independently aligns all
90 community samples to the environment table, applies a prespecified 5%
prevalence filter, and audits envfit, three Mantel distance blocks, a
non-primary partial Mantel sensitivity analysis, and adjusted-R² variation
partitioning under one fixed 999-row permutation matrix. Its 221/221 audit is
written under `results/25-environment-variance/`. Article 26 independently
closes each sample to 100%, selects the phylum Top 10 by equal-sample mean,
keeps rank-specific Unknown separate from Other known taxa, and audits Top-N,
pooled-read weighting, known-only renormalization, and four composition plots.
Its 195/195 audit is written under `results/26-community-composition/`.
Article 27 independently aggregates Phylum through Genus while preserving
missing-parent gaps, uses rank-qualified display labels and full-lineage keys,
audits within-parent and cross-rank name collisions, compares rank-specific
Top-N coverage, and keeps pooled-read and known-only calculations as
sensitivity branches. Its 209/209 audit is written under
`results/27-multirank-composition/`.
Article 28 independently defines feature-level detection from sample-level
relative abundance, with a primary threshold of 0.01% and prevalence of 80%
under an all-read denominator and equal sample weights. It separates global
from IW/CW/TW group-specific core sets, audits a 4 x 4 detection-prevalence
grid and a seed-20260728 10,000-read depth branch, and assigns mutually
exclusive primary-core, conditionally-rare, sampling-limited,
persistently-rare, and intermediate states. Its 320/320 audit, including four
English-only publication figures, is written under
`results/28-core-rare-biosphere/`.

## phyloseq GlobalPatterns beta-distance tables and tree

Article 21 uses the real `GlobalPatterns` object distributed with
`phyloseq 1.48.0`, cited to Caporaso JG et al. *PNAS*. 2011;108:4516–4522
(<https://doi.org/10.1073/pnas.1000080107>). The package license is AGPL-3.
After removing 228 all-zero features, the portable bundle contains 18,988
features, 26 samples, 28,216,678 reads, seven taxonomy ranks, English habitat
metadata, and a rooted 18,988-tip tree whose tip set exactly matches the count
table.

Regenerate `data/small/beta-distances/` from the locked installed package:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/prepare_beta_distance_data.R \
  --project-root . \
  --output-dir data/small/beta-distances
```

Exported SHA-256 values:

- `otutab.tsv`: `727073a3f26aff6730ad54b55c11e1d52bd5d3d34f99dbf75d125d2e84d9722d`
- `taxonomy.tsv`: `1189e6e2a02474a665283f6b9962ec25716b6c7faf66011a6506e87a69ddf0ff`
- `metadata.tsv`: `c708b3b89edbb26a4aab26b5e575cdefef1eae28b9b7f93cea5f8f2edfc94221`
- `rooted-tree.nwk.gz`: `535f512b3018da462a416b8f04ef36369bf1c1e383510494a790508361fc8c15`

The packaged RData source hash, dimensions, study citation, license, tree
audit, and all output hashes are recorded in
`beta-distances/source-summary.json`.

## phyloseq soilrep PERMANOVA and dispersion tables

Article 23 uses the real `soilrep` object distributed with `phyloseq 1.48.0`,
cited to Zhou J et al. *ISME Journal*. 2011;5:1303–1313
(<https://doi.org/10.1038/ismej.2011.11>). The package license is AGPL-3.
The portable bundle contains 16,825 OTUs, 56 PCR/tag technical libraries,
98,022 reads, and metadata that maps those libraries to 24 biological soil
samples in a six-block warming × clipping split-plot experiment. The packaged
object has no taxonomy table, so seven rank columns are deliberately empty and
`TaxonomyStatus` records that boundary.

Regenerate `data/small/permanova-dispersion/` from the locked installed
package:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/prepare_permanova_dispersion_data.R \
  installed \
  data/small/permanova-dispersion \
  data/small/permanova-dispersion/source-summary.json
```

Source and exported SHA-256 values:

- `soilrep.RData`: `47b25a0a033b794fca1db30ae6fdb68c55a44e6c4cd0d1be76788d1d77b51f36`
- `otutab.tsv`: `de52dc7457f15f07776c216561803afd92782af7d6723ec93deb8dc0e932b900`
- `taxonomy.tsv`: `599ed73e3d1f59bc2867782c58113100def184c3ebd2cab36286b2fc96dbaee0`
- `metadata.tsv`: `5ebf2fd60582e1a0ef80fe20d251594bec3607df11d8bb74d0981236e03f982c`

Article 23 aggregates technical libraries before inference, constructs legal
split-plot permutations from the original randomization, and reports
PERMANOVA together with PERMDISP, four planned simple effects, multiplicity
control, alternative-test direction, and preprocessing/design sensitivity.
Its 214/214 audit is written under `results/23-permanova-dispersion/`.

## decontam MUClite control tables

Articles 04–05 use the real `MUClite.rds` object distributed with
Bioconductor `decontam 1.24.0`:

- Package page:
  `https://bioconductor.org/packages/3.19/bioc/html/decontam.html`
- Object SHA-256:
  `801456d5780d5f51308d04c23f4477102cc5fbe9543ffa62fad2240a99c43d62`
- Package license: Artistic-2.0
- Dimensions: 1,951 ASVs × 569 samples
- Controls: 30 negative controls across 6 plates; 539 biological oral samples
- Method/data citation: Davis NM, Proctor DM, Holmes SP, Relman DA,
  Callahan BJ. *Microbiome*. 2018;6:226.
  <https://doi.org/10.1186/s40168-018-0605-2>

Regenerate `data/small/decontam/` from the locked installed package with:

```bash
Rscript scripts/prepare_decontam_data.R \
  installed \
  data/small/decontam \
  data/small/decontam/source_summary.json
```

Exported TSV SHA-256 values:

- `otutab.tsv`:
  `45a093f93a1e9f83e341788c043174882142c91a27b409ff2b01d2696d83624f`
- `taxonomy.tsv`:
  `da5da1ec8b56b516056457edfd1758db57044bcf3cf7de2caef1bf5a809fcc0a`
- `metadata.tsv`:
  `e180f70324eb87cebdfb53b6ede10cf279a15a4e0124e82eda6d79c0d64344f2`

The original triad retains all controls. Article 04 writes classified
features, per-sample contaminant burden, and the filtered biological triad
under `results/04-contamination/`; Article 05 independently repeats the
filter and writes batch-design, community, feature, sensitivity, and
diagnostic tables under `results/05-batch/`.

## QIIME 2 Atacama paired-end FASTQ excerpt

Articles 06–08 and 10 use a deterministic excerpt of the QIIME 2 2024.5 Atacama
soils paired-end tutorial data. Article 06 audits FASTQ structure and quality;
Article 07 reuses the same locked files for the Linux-side environment smoke
test; Article 08 imports them as `EMPPairedEndSequences`, performs maximum
Artifact validation, and verifies byte-identical FASTQ export. Article 10
independently repeats the import, reverse-complements the mapping barcodes,
demultiplexes 692 of 2,000 read pairs into 56 observed samples, audits three
QIIME 2 provenance actions, and runs FastQC 0.12.1 plus MultiQC 1.33 with
45/45 checks passing. The source is the official 1% multiplexed set:

- Study: Neilson JW, Califf K, Cardona C, et al. *mSystems*.
  2017;2(3):e00195-16.
  <https://doi.org/10.1128/mSystems.00195-16>
- QIIME 2 tutorial:
  <https://docs.qiime2.org/2024.10/tutorials/atacama-soils/>
- Source records: 135,487 synchronized forward, reverse, and barcode records
- Excerpt: 2,000 evenly spaced, order-preserving synchronized records
- Metadata: 75 samples

Official 2024.5 source SHA-256 values:

- `forward.fastq.gz`:
  `392f178b5f15967fafb9332c0493a8b6ecbf4a7857538023813a630c488c2702`
- `reverse.fastq.gz`:
  `5407370313a3974e5b70878153841a61df408e49ced133b3ace44b2db286cfae`
- `barcodes.fastq.gz`:
  `c07f517e920d9e3924ce1e87102f057232b57041947523b90d560f68a0ff3c0f`
- `sample_metadata.tsv`:
  `7cff810ad86a621ebc78a16b690255d839ea58e2622ad45a11a838d0f8c5b3bd`

Bundled excerpt SHA-256 values:

- `fastq/forward.fastq.gz`:
  `e6fabdbd31db519c719447f1caf2d1cd334c6ab932353e534ea2433ead88d254`
- `fastq/reverse.fastq.gz`:
  `d2fa5841d150fead5a73067a267ec531221de11c512a38ec32ada03f1621bf65`
- `fastq/barcodes.fastq.gz`:
  `9cce1bb9a57975d14b54a5cf07c89b2b7c2095da4a8f0d1bd2720f1349e74057`
- `fastq/metadata.tsv`:
  `8d5342f0a80abf4197088bb1ce92f4903c1bd1212f6a9182e0e1c87ebc322bcb`

Regenerate the excerpt after downloading the four fixed source files:

```bash
python3 scripts/prepare_fastq_excerpt.py \
  --forward downloads/atacama-1p/forward.fastq.gz \
  --reverse downloads/atacama-1p/reverse.fastq.gz \
  --barcodes downloads/atacama-1p/barcodes.fastq.gz \
  --metadata downloads/atacama-1p/sample_metadata.tsv \
  --output-dir data/small/fastq \
  --records 2000
```

The script validates every source checksum, FASTQ record structure, equal
sequence/quality lengths, record counts, and synchronized read identifiers.
The complete provenance and output checksums are in
`fastq/source_summary.json`.

## nf-core V4 primer-retained paired-end FASTQ

Articles 11–12 use four byte-identical FASTQ pairs from the pinned
`nf-core/test-datasets` commit
`f282ad2bafb3ca90190ec87290066ee1985df655`. The files are real MiSeq V2
V4 reads retained by the `nf-core/ampliseq` test profile; each sample has
2,500 synchronized pairs, for 10,000 input pairs in total.

- Repository:
  <https://github.com/nf-core/test-datasets/tree/f282ad2bafb3ca90190ec87290066ee1985df655/testdata>
- License: MIT
- Region/layout: 16S V4, paired-end
- Forward primer: 515F (Parada), `GTGYCAGCMGCCGCGGTAA`
- Reverse primer: 806R (Apprill), `GGACTACNVGGGTWTCTAAT`
- Local files and exact SHA-256 values:
  `primer-trimming/source_summary.json`

Article 11 first proves that the primers are present, then uses anchored,
IUPAC-aware q2-cutadapt trimming and retains 8,614 pairs. Article 12 runs five
real q2-dada2 parameter profiles on that trimmed Artifact. Its selected
220/200, maxEE 2/2 profile yields 5,213 non-chimeric pairs (60.52%) and 366
ASVs; the full 122/122 audit is written under `results/12-dada2-asv/`.

## ITSxpress tutorial ITS2 paired-end FASTQ

Article 13 uses the two real paired-end samples at the pinned
`USDA-ARS-GBRU/itsxpress-tutorial` commit
`916145d3b05e20656fde30c4732dc084de72d876` (2024-04-15):

- Repository:
  <https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/tree/916145d3b05e20656fde30c4732dc084de72d876>
- Observed data: 2 samples × 2,500 synchronized pairs = 5,000 input pairs
- Region/layout: fungal ITS2, paired-end
- Prepared collection SHA-256:
  `21b351697749dab781e8396040dbde35e35b6bb3b24dd860a91d4ecfc8ac5121`
- Provenance and per-file source/prepared SHA-256 values:
  `its2/source_summary.json`

The upstream files have `.fastq.gz` names but contain plain FASTQ bytes.
`scripts/prepare_its2_tutorial_data.py` verifies those bytes and writes
deterministic gzip streams with `mtime=0` and no embedded filename:

```bash
python3 scripts/prepare_its2_tutorial_data.py \
  --source-dir downloads/itsxpress-tutorial/data \
  --output-dir data/small/its2
```

The pinned commit does not state a separate repository or data license.
Accordingly, the bundled excerpt is used only to validate the method and file
contract, not to make ecological claims. ITSxpress 2.1.4 retains 2,065 pairs;
q2-dada2 then retains 1,753 non-chimeric pairs and infers 16 ASVs spanning
143–270 nt. The complete 55/55 audit is under `results/13-its-18s/`.

## UNITE and PR2 reference audits

Article 13 keeps reference-database provenance separate from sample
classification results.

The fungal contract uses UNITE 10.0, released 2025-02-19:

- DOI: <https://doi.org/10.15156/BIO/3301241>
- License: CC BY-SA 4.0
- Raw QIIME fungi release: 20,295 reference sequences, 81,842 representative
  sequences, and 8 taxonomy labels including Species Hypothesis `sh__`
- Pinned classifier build:
  `v10.0-2025-02-19-qiime2-2026.4`, commit
  `26ca7e07979c230ea9a65565e4781b3538962352`
- Classifier SHA-256:
  `5df2370f89b7b64766c0d969ddc5374d3951fdb8bda598ae54ff15cf86b0fca0`
- Evaluation QZV SHA-256:
  `f59c439ced1622c372a647464c69c2069ac1d37402559f888c5e4a242c7d20cd`

The pinned build removes `;sh__.*` with
`q2-rescript edit-taxonomy` before dereplication and fitting, so the raw table
and trained classifier intentionally have different rank contracts. Exact
build evidence is in `its-18s/unite-v10-classifier-build-audit.json`; the
bundled evaluation artifact is
`its-18s/eval_unite_ver2025-02-19-Q2-2026.4.qzv`.

The eukaryotic SSU contract uses the PR2 5.1.1 SSU DADA2 release:

- Release DOI: <https://doi.org/10.5281/zenodo.17458343>
- Release license: MIT
- Full source: 240,201 records, all with exactly 9 taxonomy ranks and no empty
  ranks
- Source SHA-256:
  `0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d`
- Deterministic 1,002-record audit snapshot SHA-256:
  `5f8eeee6d8b76bdc951f3620d45c750b0084491cee448553ec0baa944428add4`

`scripts/prepare_pr2_reference_audit.py` streams and audits the complete source
before writing the small stratified snapshot. Full counts, release provenance,
domain composition, length summaries, and hashes are recorded in
`its-18s/pr2-v5.1.1-source-summary.json`.

## SILVA, GTDB, and Greengenes2 V4 taxonomy references

Article 14 classifies the same 366 real Article 12 V4 ASVs against three
region-matched reference contracts. The query represents 5,213 non-chimeric
reads and has a median length of 253 nt:

- Query FASTA SHA-256:
  `fa71915be8007d36ec70f85d28401f5a94d0d33eca323782cc937cfc4bcf1af4`
- Query abundance SHA-256:
  `a80e3b94d8632c48fb6c01e788e842c022c3b1ee64b79dc94456618328895e6e`
- Source and processing audit:
  `taxonomy-databases/query-source-summary.json`

All three references were derived with the same 515F
`GTGYCAGCMGCCGCGGTAA` / 806R `GGACTACNVGGGTWTCTAAT` contract, combined
primer identity 0.80, 200–400 nt retained length, bidirectional search, and
exact-sequence least-common-ancestor dereplication:

| Reference | Full records | Primer-matched V4 | Unique V4 LCA |
|---|---:|---:|---:|
| SILVA 138.2 SSURef NR99 | 510,495 | 450,770 | 262,752 |
| GTDB R11-RS232 species-representative SSU | 93,770 | 79,669 | 51,985 |
| Greengenes2 2024.09 full-length backbone | 337,506 | 335,976 | 98,185 |

Derived sequence/taxonomy SHA-256 values:

- SILVA:
  `3a251d263bd4e12c96023c84b3f2e255b8bcb0d5865dc6f299f9918e6be3a808`
  /
  `f49320aaa32aa70af5c5a48649786e3bb4177868da5da4d2143d0a54c5fa47b0`
- GTDB:
  `3c4234b4defcdc3d8c3b4d1e8956b12288bccd853a4b7f3a3003847c32e6e250`
  /
  `a23d0d4f62183e949c868c3f25aace3acdebb9362ae0bb1d0e9fe81e74c43fe0`
- Greengenes2:
  `330fe7f400c84c1e454ca9193bcc49f078b01c7880a8bfb99cbe4c0c3c1c14d2`
  /
  `9594e9ae78a568672891693207fd4dcc31f8eada83b931031642fafcff2684bd`

The immutable full releases remain ignored scratch under
`data/raw/taxonomy-databases/`. Rebuild the small region assets after
downloading and checksum-verifying those official files:

```bash
python3 scripts/prepare_taxonomy_database_assets.py \
  --project-root . \
  --raw-dir data/raw/taxonomy-databases \
  --work-dir work/14-reference-preparation \
  --output-dir data/small/taxonomy-databases \
  --qiime-env microbiome-qiime2-2026.4 \
  --threads 4
```

`scripts/validate_taxonomy_databases.py` performs fresh QIIME 2 imports,
maximum validation, and the same bounded consensus-VSEARCH classification for
all three databases. Its fixed 91/91 audit assigns 358, 356, and 362 of the
366 ASVs with SILVA, GTDB, and Greengenes2, respectively. Database species
strings remain reference-dependent labels and are not treated as proof of
species or strain.

## SILVA V4 region-specific classifier validation

Article 15 reuses the checksum-locked SILVA 138.2 V4 LCA reference and real
366-ASV query above. It groups all 262,752 reference records by complete
domain-to-genus lineage, applies a fixed SHA-256 ordering inside each group,
keeps singleton lineages in training, and creates a 236,516/26,236
training/holdout split with no shared identifier or exact sequence.

The holdout is deliberately a same-release known-lineage interpolation test:
all 5,161 holdout lineages have training representatives. It is not an
independent mock-community or novel-lineage accuracy estimate. Genus
precision/recall/F1 at confidence 0.50, 0.70, and 0.90 are
0.913/0.843/0.876, 0.938/0.800/0.864, and 0.961/0.731/0.831, respectively.

After threshold evaluation, the production classifier is refitted on all
262,752 references under QIIME 2/q2-feature-classifier 2026.4.0 and
scikit-learn 1.7.1. At confidence 0.70, 318 of 366 query ASVs representing
4,217 of 5,213 reads receive a genus label. The 80/80 audit, split membership,
metrics, model artifacts, provenance, and command log are written under
`results/15-train-classifier/` by
`scripts/validate_region_classifier.py`.

For single-article reproducibility, the validator also recomputes the bounded
SILVA consensus-VSEARCH comparator inside the Article 15 step. It does not
consume the taxonomy Artifact produced by Article 14.

## Greengenes 13_8 SEPP reference-placement validation

Article 16 reuses `taxonomy-databases/query-v4-asvs.fasta.gz` and
`taxonomy-databases/query-v4-asv-abundance.tsv`: 366 real 515F/806R V4 ASVs
representing 5,213 non-chimeric reads. The abundance file is aggregated over
the four source samples and is used only for ASV/read retention auditing; it
is not a sample-by-feature ecology table.

The external placement reference is downloaded separately to
`downloads/sepp-refs-gg-13-8.qza` from QIIME 2 public data. Its fixed contract
is:

- SHA-256:
  `e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360`
- Artifact UUID: `a14c6180-506b-4ecb-bacb-9cb30bc3044b`
- semantic type / format: `SeppReferenceDatabase` / `SeppReferenceDirFmt`
- embedded alignment/tree: 203,452 reference IDs and 1,285 alignment columns

`scripts/validate_sepp_phylogeny.py` runs a fresh QIIME 2
q2-fragment-insertion placement with fixed 1,000/5,000 alignment/placement
subsets, filters the aggregate feature table to placed tips, and reconciles
FASTA, jplace, rooted-tree, and BIOM IDs. The validated output contains all
366 query tips in a 203,818-tip rooted tree and retains all 5,213 reads. The
74/74 audit and placement diagnostics are written under
`results/16-sepp-phylogeny/`; the official reference Artifact remains an
external download and is not stored under `data/small/`.

## QIIME 2 to phyloseq import bundle

Article 17 reads only the independently frozen files under
`phyloseq-import/`. The bundle contains a 366-feature × 4-library integer
count table with 5,213 reads, a 366 × 7 SILVA taxonomy table, 4 × 5 technical
metadata, 366 representative V4 sequences, taxonomy confidence, and the
complete 203,818-tip rooted SEPP insertion tree. The tree includes 203,452
reference-only tips and 366 query tips; the R workflow explicitly prunes it
to the observed ASVs before constructing the object.

The four library names and `S103`/`S115` filename sets are technical file
labels, not biological treatment groups. No ecological group comparison is
supported by this small bundle.

Frozen file SHA-256 values:

- `otutab.tsv`:
  `6c43bb8f929f935b60e42ed32737632d3be18ddd491f5298ba3a6abd6e3f18a9`
- `taxonomy.tsv`:
  `8ef235ea19923307eef4d8e1905c2f34b78cf25cf1a83f08aa00353f2c5c0d38`
- `metadata.tsv`:
  `6628dd51051c53884e82c540ceaff8cc34f69900a78d372239a2ed41ff5b7263`
- `taxonomy-confidence.tsv`:
  `26fd4c66301f97c6feac8f3302ab19f269e9429c9ce65dea5c66ae277eca55b7`
- `representative-sequences.fasta.gz`:
  `a19e083cfb0ab757d557569bdf76aa3de78576ae35033a7e69a647976a370bc2`
- `rooted-insertion-tree.nwk.gz`:
  `7e541e2ed9e3f81d9b6c6bf98999335e7f684a37ac5dd5f2b56e61b2f4c4c47c`
- `source-summary.json`:
  `358710a72413e1bc1ff0e13cb1a3722a848583d37321e1c2442523f09d3226fd`

Regenerate the frozen bundle from the validated source Artifacts with the
locked QIIME 2 Python environment:

```bash
conda run -n microbiome-qiime2-2026.4 python \
  scripts/prepare_phyloseq_import_assets.py \
  --table-artifact results/12-dada2-asv/dada2-table.qza \
  --taxonomy-artifact results/15-train-classifier/query-taxonomy-c070.qza \
  --sequences-artifact results/12-dada2-asv/dada2-rep-seqs.qza \
  --tree-artifact results/16-sepp-phylogeny/sepp-insertion-tree.qza \
  --output-dir data/small/phyloseq-import
```

Bundle regeneration records upstream provenance; Article 17 validation does
not read those earlier result directories. It starts from the frozen files
above and runs:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_phyloseq_import.R \
  --project-root . \
  --input-dir data/small/phyloseq-import \
  --output-dir results/17-phyloseq-import \
  --figure-dir figures
```

The validator passes 82/82 checks and writes a five-component `phyloseq`
object plus readable component, ID, tree, sequence, taxonomy, library, and RDS
round-trip audits under `results/17-phyloseq-import/`.

## Article 18 publication-graphics inputs

Article 18 returns to the independently readable root triad rather than using
an earlier R object: `otutab.tsv`, `taxonomy.tsv`, and `metadata.tsv`, together
with `source_summary.json`. These are the same real wetland 16S tables from the
checksum-locked `microeco 2.0.0` source archive: 13,628 features, 90 samples,
and 1,619,670 reads, with 30 samples in each CW/IW/TW group.

The validator fixes the following input SHA-256 values:

- `otutab.tsv`:
  `76fa79c38da889f35978dc86da4641a270961746709ff38049ee5f67e3c6f7a3`
- `taxonomy.tsv`:
  `725280bb9a0cd9bda7b540022e92af945ceed52527f8b2d220b055b4489f6901`
- `metadata.tsv`:
  `df24771dccf27607ddf922c6bca2cafa876d946fbe2e09d14b601accce66ba64`
- `source_summary.json`:
  `e51bb47c6766d41ee6efe34689ce550ab50fa4ec07ca3beef5656268124bf6f6`

Run the independent graphics audit with:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_publication_graphics.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/18-publication-graphics \
  --figure-dir figures
```

It regenerates four original figures in PDF/SVG/PNG/LZW-TIFF, three raster
resolution probes, and machine-readable font, palette, dimension, format, and
package audits. The current fixed run passes 171/171 checks.

## Article 19 alpha-diversity inputs

Article 19 independently rereads the same root `otutab.tsv`, `taxonomy.tsv`,
`metadata.tsv`, and `source_summary.json`; it does not consume an RDS object or
an Article 18 result. The fixed table contains 13,628 features, 90 samples,
1,619,670 reads, and 30 samples in each CW/IW/TW group.

The analysis preserves singleton and doubleton information before calculating
Observed, bias-corrected Chao1, Shannon, Simpson, Pielou, Good's coverage, and
Hill q=0/1/2. It compares analytical diversity estimates at 10,000 reads and
at 90% sample coverage with `iNEXT 3.0.2`; it does not retain a randomly
rarefied count table for downstream models. The 90% coverage target requires
56 interpolations and 34 bounded extrapolations, with a maximum effort ratio
of 1.453.

Run the independent audit with:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_alpha_diversity.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/19-alpha-diversity \
  --figure-dir figures
```

The validator regenerates four original figures in PDF/SVG/600-ppi PNG/LZW-
TIFF and writes raw, equal-size, equal-coverage, sensitivity, rank-stability,
format, and session audits. The current fixed run passes 193/193 checks.

## Article 29 community-typing DMM inputs

Article 29 starts from the independently readable
`community-typing-dmm/otutab.tsv`, `taxonomy.tsv`, and `metadata.tsv` triad.
They are deterministically exported from `Twins.csv` and `TwinStudy.t` in
`DirichletMultinomial 1.46.0`, which distributes the real Holmes/Turnbaugh
twins gut 16S genus-count data. The fixed table contains 130 source
categories, 278 samples, and 570,851 raw reads. BMI class is retained only as
post-fit descriptive metadata; suffix-derived participant keys are audit
fields and are not asserted to be clinically validated visit labels.

Frozen file SHA-256 values:

- `otutab.tsv`:
  `28de822c434aedade101e66ac60a61df38124c9b59357095599fbd278f2ecb41`
- `taxonomy.tsv`:
  `9eaa0bf788ad1db8c314a59e774c8935bdaa33a51135eaf606c11d5a497dff08`
- `metadata.tsv`:
  `0368c9129c7d77f381688465141593a2c61a229368d2fc411c9017e1abe04409`

Regenerate the triad from the locked package and run the independent audit:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/prepare_dmm_twins_data.R \
  --output-dir data/small/community-typing-dmm

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_community_typing_dmm.R \
  --project-root . \
  --input-dir data/small/community-typing-dmm \
  --output-dir results/29-community-typing-dmm \
  --figure-dir figures
```

The validator uses unrarefied non-negative integer counts, fits candidate
`k=1–7` models, launches one fresh R process per model for the initialization
audit, and regenerates four original figures in PDF/SVG/600-ppi PNG/LZW-TIFF.
The current fixed run passes 222/222 checks. Data provenance and licensing are
recorded in `community-typing-dmm/source-summary.json`.

## Articles 30–33 and 35–40 downstream inputs

Articles 30–33 and 35–40 independently reread the root wetland
`otutab.tsv`/`taxonomy.tsv`/`metadata.tsv` triad documented under Article 18:
13,628 ASVs, 90 samples, 1,619,670 reads, and 30 samples per CW/IW/TW group.
They do not read another article's RDS or result table. Article 30 audits
closure and reference frames; Articles 31–33 aggregate by complete lineage for
differential-abundance analyses and graphics; Article 35 applies a group-blind
genus-node filter before SparCC and SPIEC-EASI inference. Articles 36–37 extend
the same group-blind universe to topology robustness and microbial WGCNA;
Articles 39–40 fit neutral and macroecological models directly from the triad.
Articles 37, 38 and 40 additionally read the matched `environment.tsv`.

Article 38 also uses `rooted-tree.nwk.gz`, deterministically extracted from the
real `phylo_tree_16S` object in `microeco 2.0.0`. The source tree is rooted and
has branch lengths, 14,096 tips and 14,095 internal nodes. It contains all
13,628 count-table OTUs plus 468 explicitly pruned extra tips. Frozen SHA-256:

- `rooted-tree.nwk.gz`:
  `05d64719bfe720714fdf03f5158893f54fc6400ab649de807cbecc5599589f61`
- updated `source_summary.json`:
  `acc15b18d3f85d6d35770d0db7580d91d0a55a838862876536500a8d7c75711b`

Regenerate the four tables and tree from the fixed source tarball with:

```bash
Rscript scripts/prepare_pilot_data.R \
  /path/to/microeco_2.0.0.tar.gz \
  data/small \
  data/small/source_summary.json
```

## Article 34 absolute-quantification inputs

`absolute-quantification/` contains an independently readable triad plus
`cell-load.tsv`. The files are deterministically exported by
`scripts/prepare_vandeputte_absolute_data.R` from the `Vandeputte` object in
Bioconductor `reconsi 1.16.0`. They contain 234 genera, 135 samples, 4,080,996
reads, and matched flow-cytometry microbial loads. The main comparison uses
only the common Disease cohort (29 Crohn's disease and 66 healthy samples);
40 healthy samples from the separate Study cohort remain available for the
cohort audit but are not pooled into the contrast.

Frozen file SHA-256 values:

- `absolute-quantification/otutab.tsv`:
  `2fb2e5042b448ea8db38559fe223b5dca6be7c1b686059451cc5e6a7127e3af0`
- `absolute-quantification/taxonomy.tsv`:
  `1518574260016b621f469abf8cfc2982f0a137194ab11a5faa6027cf26676b3d`
- `absolute-quantification/metadata.tsv`:
  `7201f8435f51833338cb8c2e660b7c5d061e5036068977b5bb8804aabeb59a81`
- `absolute-quantification/cell-load.tsv`:
  `4a8e793dcdb76a80a85e316ce4a3f810aa6be2719bc060fe431abc66946360ef`
- `absolute-quantification/source-summary.json`:
  `8ac35a7119d1772d5d1c92f30c051f2fb3bdb2b5d3dca3f83c885528352d9d94`

Regenerate the bundle with:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/prepare_vandeputte_absolute_data.R \
  data/small/absolute-quantification
```

The underlying study is Vandeputte et al., *Nature* 2017,
doi:10.1038/nature24460; package and file-level provenance are recorded in
`absolute-quantification/source-summary.json`.

## Articles 41–45 functional prediction, external validation, and survival inputs

Article 41 uses the official PICRUSt2 chemerin tutorial archive, converted to
an independently readable 37-ASV × 24-sample triad plus the exact
representative-sequence FASTA and BIOM input. The table contains 108,718
reads; all 37 sequences entered PICRUSt2 2.6.3, selected the bacterial domain,
and passed the prespecified maximum NSTI of 2.0. The upstream archive SHA-256
is `f57abdc069b6560f0ddb739cf9a341f5679b12f4947f4db2ad5e2cab893669e9`.
The official archive does not include taxonomy, so the seven taxonomy ranks
are intentionally blank and are not fabricated. File-level checksums and this
boundary are recorded in `picrust2-chemerin/source-summary.json`.

Article 42 independently rereads the root microeco wetland triad documented
above. Article 43 independently rereads the Vandeputte triad and matched
`cell-load.tsv`; its classification branch uses only the 95-sample Disease
cohort, while its regression branch uses the measured microbial load.

Article 44 uses three standardized CRC 16S cohorts from MicrobiomeHD
(Zenodo record `10.5281/zenodo.1146764`), retaining only CRC and healthy-control
stool samples. Xiang contains 43 samples (21 CRC/22 control), Zhao contains
102 (46/56), and Zackular contains 60 (30/30). Their fixed upstream archive
SHA-256 values are `c16052495cc717069c670cabc9c1028750a71b7fd04ac3d28dbc1e560521d63e`,
`2f479dab25d980f2d295244166165f5531b9d28e38896767d1acf6a08e88f91f`,
and `294305422293a2e87f9eb61d89612f49fb08b295d779b1cce7996b6d884c4498`.
The derived cohort triads and file-level checksums are under
`cross-cohort-crc/`; the source record is CC BY-NC 4.0.

Article 45 uses the MiSurv-processed NOD-mouse T1D cohort at fixed Git commit
`692a2ac7079d2dbf581893b9a820c80cde0a7e31`. The triad contains 348 OTUs,
173 baseline samples, 3,073,108 reads, 118 T1D-onset events, and 55 censored
mice. Sampling occurred at week 6 (n=6), 7 (n=164), or 8 (n=3), and each
mouse's follow-up time is calculated from its own sampling week. The event is
mouse T1D onset, not human death or overall survival. The source repository
does not state a data license; cite Zhang et al. 2018 and Gu et al. 2023 and
verify redistribution terms before republishing the source tables. Checksums
and boundaries are recorded in `survival-t1d/source-summary.json`.

Regenerate all three fixed data bundles from checksum-verified downloads with:

```bash
bash scripts/download_articles_41_45_data.sh
```

Re-run the locked PICRUSt2 analysis and its R audit figures with:

```bash
bash scripts/run_article41_picrust2.sh
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/render_article41_picrust2.R "$PWD"
```

## Articles 46–50 multi-omics, multi-kingdom, and source-tracking inputs

`paired-ibd-multiomics/` fixes the Franzosa IBD microbiome–metabolome resource
at curated-data commit `89a519d8c832008fbc6e650453e83e2f04858d02` and contains
250 microbial genera × 220 strictly paired samples, 277 named metabolites,
clinical metadata, and metabolite annotations. Articles 46–48 reread these
files independently for global concordance/HAllA, sPLS/DIABLO, and native
MMvec/MOFA analyses.

`multi-kingdom-duran/` fixes the Duran Arabidopsis supplement at commit
`6db5e85cc5d442fd95fcdcb7250b72fa9e2ff900`. Its 36 strictly paired samples
form a balanced 3-soil × 3-compartment × 4-replicate design with bacterial,
fungal, and oomycete count/taxonomy tables.

`source-tracking-feast/` fixes the official FEAST demonstration at commit
`2f8f3df8051e0e08341f597a9f4693bfb76b3bf6`: one sink, nine candidate-source
samples, four source classes, and 1,839 anonymous features. Anonymous features
remain anonymous; the tutorial does not invent taxonomy. Each directory's
`source-summary.json` records upstream URLs, checksums, transformations, and
reuse boundaries.

## Articles 51–55 longitudinal and causal-inference inputs

`longitudinal-dietswap/` contains the O'Keefe DietSwap object distributed with
Bioconductor `microbiome 1.26.0`, serialized as an independent triad with 130
features, 222 observations, 38 subjects, and six timepoints. The packaged
source SHA-256 is
`67cc6d117ce9dfbbc472fba828e1bb0d48c621fb574f4f8d3424f6cc9a77fdab`.
The assay is a HITChip 16S phylogenetic microarray rather than amplicon
sequencing; the repeated-measures workflow transfers, but assay preprocessing
does not. Cite O'Keefe et al. 2015 and Dryad `10.5061/dryad.1mn1n`.

Article 52 rereads the root wetland triad and `environment.tsv`; all SEM
variables are measured environmental or community summaries rather than
simulated covariates. Article 53 rereads `paired-ibd-multiomics/` and applies a
prespecified complete-case contract to 108 samples; no intermediate object
from Articles 46–48 is required.

`mr-mibiogen/` contains TwoSampleMR-format exposure/outcome tables and the API
selection ledger for MiBioGen `genus.Bifidobacterium.id.436` and IBD study
GCST004131. The MiBioGen source SHA-256 is
`37001a83d060596fe0b97b63d6a397f01f43a29add2925d406916b7a50b5883e`.
Six cross-chromosome candidate instruments are frozen; the lack of ancestry-
matched LD clumping, allele frequencies, Steiger inputs, and colocalisation
statistics is explicitly retained as an analysis boundary.

`causal-evidence/evidence-cases.tsv` is a seven-row, DOI-linked ledger of
primary microbiome studies used by the no-code Article 55 reading framework.
It records design capacity, bounded claims, and residual threats; it is not a
numerical evidence score. Its SHA-256 is
`a130da9f84a017a6224cf7b467aa49f898367bd0ce2886e30607f2921a333ab2`.

Regenerate or audit the fixed bundles with:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/prepare_articles_51_55_data.R
python3 scripts/prepare_article54_mr_data.py
```
