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
written under `results/25-environment-variance/`.

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
