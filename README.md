# microbiome-best-practices

一套面向“懂生物、不懂生信”科研读者的 55 篇中文 16S 微生物组分析教程。网页版使用 Quarto Book，微信公众号文章从同一份 QMD、真实数据、执行结果和原创重绘图派生。

## 发布规则

- 网页版与微信公众号版都必须单篇自足。
- 每篇 R 分析从 `otutab`、`taxonomy`、`metadata` 三件套开始；需要环境、临床或外部统计数据时，在本篇明确提供。
- 论文用于方法、数据和研究设计溯源；默认使用公开数据重新计算并重绘，不直接复制受版权保护的论文原图。
- 下游章节必须真实执行；只有 `qa_report.json.status == "passed"` 才能生成发布包。
- 默认只生成本地微信公众号 bundle；只有获得明确授权后才调用草稿 API，发布与群发仍需另行授权。
- GitHub Pages 的 deploy job 默认关闭；只有人工设置仓库变量 `ENABLE_PRODUCTION_DEPLOY=true` 后才允许部署。

## 当前状态

- `tutorial.yaml` 已固定 01–55 的篇目、文件名、数据、执行步骤和发布边界。
- 第 01–25 篇已完成；第 25 篇从真实湿地四表出发，使用同一固定置换矩阵运行 envfit、三个预设 Mantel blocks、带警告的 partial Mantel sensitivity 和 adjusted-R² VPA，并重绘四张出版图。
- Pilot QA 已通过：固定真实数据及参考资产、55 篇目录、二十五篇 HTML、273 个 PDF/SVG/PNG/TIFF 图形文件和 336 份结果/审计文件均通过，21/21 个 workflow steps 与 640/640 个发布断言全部成功；run key 为 `4343439bcad30e0a`，当前 manifest hash 为 `3397de082f938fe8679e8a5437379980666d9a39931dcc35116d0512f52a4f73`。
- 根目录 55 篇 Quarto Book 已完整渲染并产生恰好 55 个 HTML 页面；其余 30 篇仍是 `draft: true` 的目录占位，不代表正文已经完成。
- 经明确授权，第 01–25 篇已从同一份已通过 QA 的 Quarto 输出生成公众号审阅稿：102 张正文图和 25 张封面均已上传，25/25 篇草稿创建成功，草稿总数由 53 增至 78；未调用发布或群发接口。本地可追溯记录见 `rendered/wechat_review_01_25/report.json` 与 `rendered/wechat_review_01_25/live_report.json`（生成目录默认不进 Git）。
- GitHub 审阅仓库为 `petemeng/microbiome-best-practices`，第 01–25 篇位于 Draft PR [#1](https://github.com/petemeng/microbiome-best-practices/pull/1)；该 PR 尚未合并，GitHub Pages 生产部署仍关闭。
- QIIME 2 官方 2026.4 环境、固定 ITSxpress overlay 与 gemelli 0.0.13 独立环境均已在本机验证；第 07–17 篇相应环境、Artifact、导入、原始质控、去引物、DADA2、ITS/18S、数据库比较、区域 classifier、SEPP rooted phylogeny 与 R 对象验收均已通过，第 18–25 篇的出版图形、Alpha 多样性、组间检验、Beta 距离、PCoA、PERMANOVA/PERMDISP、CAP 与环境方差解构合同也已独立通过。

## 已锁定的真实数据

第 01、02、03、09、18、19、22、24、25 篇使用 CRAN 归档中的 `microeco 2.0.0`：

- 源码包：`microeco_2.0.0.tar.gz`
- SHA-256：`454a3b71ceeea86bdd54f475a12b5bac9c3b124f663b8dc6e560827fca89ebfd`
- 包许可证：GPL-3
- 表格规模：13,628 个 OTU、90 个 16S 样本、taxonomy、sample metadata 和 200 行环境观测
- 数据研究：[An et al., *Geoderma*, 2019](https://doi.org/10.1016/j.geoderma.2018.09.035)
- 软件来源：[Liu et al., *FEMS Microbiology Ecology*, 2021](https://doi.org/10.1093/femsec/fiaa255)

`data/small/` 保存可直接渲染的标准文本表；`scripts/prepare_pilot_data.R` 可从固定源码包重新生成它们。

第 04–05 篇使用 Bioconductor `decontam 1.24.0` 随包分发的真实 `MUClite.rds`：

- 对象 SHA-256：`801456d5780d5f51308d04c23f4477102cc5fbe9543ffa62fad2240a99c43d62`
- 表格规模：1,951 个 ASV、539 个口腔生物样本、30 个阴性对照、6 块实验板
- 软件/方法来源：[Davis et al., *Microbiome*, 2018](https://doi.org/10.1186/s40168-018-0605-2)
- 导出与校验：`scripts/prepare_decontam_data.R`

第 06–08、10 篇使用 QIIME 2 官方 Atacama soils 2024.5 paired-end 1% 数据的确定性真实摘录：

- 官方完整 1% 文件：135,487 组同步 forward、reverse 与 12-nt barcode records
- 本地摘录：2,000 组首尾覆盖、等间距、保序同步 records；metadata 为 75 个样本
- 数据研究：[Neilson et al., *mSystems*, 2017](https://doi.org/10.1128/mSystems.00195-16)
- 官方教程：[QIIME 2 Atacama soils](https://docs.qiime2.org/2024.10/tutorials/atacama-soils/)
- 生成与校验：`scripts/prepare_fastq_excerpt.py`
- 完整 source/output SHA-256：`data/small/fastq/source_summary.json`

第 11–12 篇使用 nf-core/test-datasets 固定提交中的真实 MiSeq V2、V4 paired-end FASTQ：

- 固定提交：`f282ad2bafb3ca90190ec87290066ee1985df655`
- 数据规模：4 个样本，每个方向每样本 2,500 条 reads，共 10,000 对输入
- 引物：515F (Parada) / 806R (Apprill)，原始 reads 保留引物
- 第 11 篇锚定去引物后保留 8,614 对；第 12 篇固定 220/200、maxEE 2/2 后保留 5,213 对 non-chimeric reads，得到 366 个 ASV
- 数据与完整 SHA-256：`data/small/primer-trimming/source_summary.json`
- 来源：[nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/f282ad2bafb3ca90190ec87290066ee1985df655/testdata)

第 13 篇使用 ITSxpress 官方教程仓库固定提交中的两份真实 ITS2 paired-end 样本：

- 固定提交：`916145d3b05e20656fde30c4732dc084de72d876`（2024-04-15）
- 数据规模：2 个样本，每样本 2,500 对 reads，共 5,000 对输入
- 固定 gzip 集合 SHA-256：`21b351697749dab781e8396040dbde35e35b6bb3b24dd860a91d4ecfc8ac5121`
- ITSxpress 保留 2,065 对，q2-dada2 最终保留 1,753 对 non-chimeric reads，得到 16 个 ASV
- 原仓库在该提交未声明独立的数据许可证，因此本地固定摘录仅用于方法与输入合同验证，不作生态结论
- 数据、逐文件 SHA-256 与 gzip 修复说明：`data/small/its2/source_summary.json`
- 来源：[ITSxpress tutorial](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/tree/916145d3b05e20656fde30c4732dc084de72d876)

同篇还锁定两套参考资产：

- UNITE 10.0（2025-02-19，DOI `10.15156/BIO/3301241`）：原始 QIIME fungi taxonomy 有 8 个标签并保留 `sh__`，固定的 QIIME 2 2026.4 分类器在训练前移除 `sh__`，分类器 SHA-256 为 `5df2370f89b7b64766c0d969ddc5374d3951fdb8bda598ae54ff15cf86b0fca0`
- PR2 5.1.1（DOI `10.5281/zenodo.17458343`）：完整 SSU DADA2 FASTA 有 240,201 条记录，全部具有 9 层 taxonomy；源文件 SHA-256 为 `0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d`
- 审计资产：`data/small/its-18s/`

第 14 篇使用第 12 篇真实 DADA2 输出独立打包的 366 条 V4 ASV 与 5,213 条 reads 频数，并锁定三套同区域参考：

- query FASTA SHA-256：`fa71915be8007d36ec70f85d28401f5a94d0d33eca323782cc937cfc4bcf1af4`
- SILVA 138.2 SSURef NR99：510,495 条完整记录，450,770 条 primer-matched V4，262,752 条 unique V4 LCA
- GTDB R11-RS232 species-representative SSU：93,770 条完整记录，79,669 条 primer-matched V4，51,985 条 unique V4 LCA
- Greengenes2 2024.09 full-length backbone：337,506 条完整记录，335,976 条 primer-matched V4，98,185 条 unique V4 LCA
- 三库统一使用 515F/806R、primer identity 0.80、200–400 nt、both orientation、exact-sequence LCA 和相同 consensus-VSEARCH 参数
- 91/91 验收分别注释 358、356、362 个 ASV；species 字符串只作为数据库条件标签，不作为物种或菌株证明
- 来源、许可证、完整/区域 checksum 与 rank 审计：`data/small/taxonomy-databases/`

第 15 篇复用其中 checksum-locked SILVA 138.2 V4 LCA reference 与真实 query：

- 262,752 条 unique V4 reference 按完整 genus lineage 和固定 SHA-256 顺序拆为 236,516 条 training 与 26,236 条 holdout；711 个 singleton lineage 只进入 training
- 训练集与 holdout 没有共享 ID 或 exact sequence；留出集是同一 release 的 known-lineage interpolation，不是外部 mock-community accuracy
- confidence 0.50/0.70/0.90 的 genus precision/recall/F1 分别为 0.913/0.843/0.876、0.938/0.800/0.864、0.961/0.731/0.831
- 用全部 262,752 条 reference 重训生产 classifier；0.70 下 318/366 个 ASV、4,217/5,213 条 reads 获得 genus 标签
- 与 SILVA sequence-search 在 comparable genus reads 中 exact-name agreement 为 97.74%；只作为方法敏感性
- sequence-search 对照在第 15 篇验证步骤内用同一 SILVA 区域参考重新计算，不读取第 14 篇的分类 Artifact
- 80/80 验收覆盖环境版本、输入 checksum、拆分、训练、六次分类、11 个 maximum-validated QIIME 对象、provenance 和四张英文证据图

第 16 篇复用同一组 366 条真实 V4 ASV 与 5,213 条 reads 聚合频数，并锁定 QIIME 2 官方 Greengenes 13_8 99% SEPP reference：

- reference Artifact：`sepp-refs-gg-13-8.qza`，50,161,069 bytes，SHA-256 `e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360`
- Artifact UUID：`a14c6180-506b-4ecb-bacb-9cb30bc3044b`；语义类型 `SeppReferenceDatabase`；内含 203,452 条 1,285-column aligned SSU 与 203,452-tip rooted tree
- 当前运行时固定 QIIME 2/q2-fragment-insertion 2026.4.0、SEPP 4.5.6、HMMER 3.4、pplacer 1.1.alpha19，alignment/placement subset 为 1,000/5,000
- 输出 rooted tree 有 203,818 tips，完整保留 203,452 个 reference tips 与 366 个 query tips；无缺失或负 branch lengths
- filter-features 保留 366/366 ASV 与 5,213/5,213 reads；131 个 ASV 的 top likelihood weight 低于 0.50，因此“全部插入”不解释为“全部可靠定位”
- 74/74 验收覆盖 reference contents、环境、maximum validation、provenance、jplace、tip reconciliation、branch lengths、feature/read retention 和四张英文证据图
- 官方 reference：[QIIME 2 public data](https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza)；方法：[Janssen et al., *mSystems*, 2018](https://doi.org/10.1128/mSystems.00021-18)

第 17 篇使用独立冻结的 QIIME 2 → R 输入包 `data/small/phyloseq-import/`：

- `otutab.tsv` 为 366 features × 4 technical libraries 的非负整数矩阵，总计 5,213 reads；四个文件标签不作为生物学分组
- `taxonomy.tsv` 为 366 × 7；SILVA 138.2 confidence 0.70 在 Genus 层覆盖 318 个 ASV、4,217 条 reads，Species 保持 0
- 366 条代表序列长度为 252–303 nt；完整 rooted SEPP insertion tree 有 203,818 tips，其中 203,452 个 reference-only tips 与 366 个 query tips
- R 中先核对方向、唯一 ID 和集合差，再用 `ape::keep.tip()` 显式裁剪为 366-tip 对象树，构建 `otu_table`、`tax_table`、`sample_data`、`phy_tree`、`refseq` 五组件对象
- 82/82 检查及 RDS 往返全部通过；对象保留 366 features、4 samples 和 5,213 reads
- 七份输入的逐文件 SHA-256、来源与许可证记录在 `data/small/phyloseq-import/source-summary.json`

第 18 篇重新从根目录真实湿地 `otutab`、`taxonomy`、`metadata` 三件套开始：

- 输入为 13,628 features × 90 samples、总计 1,619,670 reads，CW/IW/TW 各 30 个样本；四份输入及来源摘要均通过 SHA-256 校验
- 固定 DejaVu Sans、三组色盲友好颜色与形状冗余编码，并审计 standard/deutan/protan 三种视图
- 四张原创图各输出 PDF、SVG、PNG、LZW TIFF；矢量文件保留文本或嵌入字体，位图按终稿 89/183 mm 物理尺寸以 600 ppi 生成
- 89 × 70 mm 分辨率探针实际得到 72/300/600 ppi 对应的 252 × 198、1,051 × 826、2,102 × 1,653 pixels
- 171/171 检查覆盖输入、包版本、字体、调色板、物理尺寸、像素尺寸、有效 ppi、格式签名与 TIFF 压缩；明细在 `results/18-publication-graphics/`

第 19 篇从同一 checksum-locked 三件套独立计算 Alpha 多样性：

- 原始 library size 为 10,364–37,374 reads，Observed richness 与 depth 的 Spearman 相关为 0.781；sample coverage 为 86.25%–95.32%
- 10,000-read 解析稀释保留 90/90 个样本，但相当于不利用 44.4% 的已测 reads；不保存随机抽平 count table
- 0.90 coverage 标准化包含 56 个 interpolation 与 34 个短程 extrapolation 样本，最大 effort ratio 为 1.453
- Hill q=0/1/2 分别报告 observed richness、exp(Shannon) 与 inverse Simpson 的有效 feature 数；组别汇总仅作描述
- 193/193 检查覆盖输入、包版本、指标公式、`phyloseq`/`vegan` 一致性、coverage sensitivity、图形格式与物理尺寸；明细在 `results/19-alpha-diversity/`

第 20 篇再次从三件套独立重算 Alpha 指标，并把设计选择、检验族和效应量固定为可审计合同：

- metadata 含 `Group`、`Type`、`Saline`，但没有受试者、配对、时间或区组标识；当前数据只执行 30 + 30 + 30 的独立三组分支，同时明确保留配对和重复测量模板
- 主分析使用 0.90 coverage Hill q=0/1/2；三个 Kruskal–Wallis 检验经 Holm 校正后均不显著，最大计划成对效应为 Inland − Tibetan 的 q=0 Cliff's delta 0.258，95% bootstrap CI 为 −0.033–0.547
- 原始 q=0 的名义 P=0.0346，但 10,000-read 与 0.90 coverage 口径分别为 P=0.511 和 0.216；九项标准化敏感性经统一 Holm 校正后均不显著
- 188/188 检查覆盖输入、设计门禁、三种标准化、主检验、多重比较、5,000 次 bootstrap、Welch 敏感性与四张图的 PDF/SVG/PNG/LZW-TIFF 文件；明细在 `results/20-alpha-group-tests/`

第 21 篇使用 `phyloseq 1.48.0` 随包分发的真实 `GlobalPatterns` 数据，并导出单篇可读的三件套与根树：

- 删除 228 个全零 feature 后，输入为 18,988 features × 26 samples、28,216,678 reads；rooted tree 有 18,988 个唯一 tips，tip 集与计数表精确相等且 branch length 有效
- Bray-Curtis、binary Jaccard 与 weighted/unweighted UniFrac 使用 seed 20260722 的每样本 50,000-read 稀释表；26 个样本全部保留，总计 1,300,000 reads
- Aitchison 与 Robust Aitchison 从未稀释计数开始，统一使用 total count > 10、prevalence > 10% 的 11,063-feature 集；经典 CLR 明报 pseudocount 0.5，并审计 0.1/0.5/1.0
- Robust Aitchison 固定 Python 3.12.13、gemelli 0.0.13、3 components、5 iterations 与 ARPACK seed；两次独立运行的五个核心输出逐字节一致
- weighted UniFrac 的 3 个小负特征值只在该面板使用 Lingoes 校正；其余五类距离保留原始 PCoA 口径
- 239/239 检查覆盖输入哈希、树、标准化、零值、六个距离矩阵、PCoA、双跑确定性和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/21-beta-distances/`

第 23 篇使用 `phyloseq 1.48.0` 随包分发的真实 `soilrep` 数据，并导出独立三件套：

- 源 `soilrep.RData` SHA-256 为 `47b25a0a033b794fca1db30ae6fdb68c55a44e6c4cd0d1be76788d1d77b51f36`；数据研究为 [Zhou et al., *ISME Journal*, 2011](https://doi.org/10.1038/ismej.2011.11)
- 输入含 16,825 个 OTU、56 个 PCR/tag 技术文库和 98,022 reads；推断前按 `BiologicalSampleID` 聚合为 24 个土壤样本、6 个区组、12 个 main plots 和四种处理组合
- 源对象没有 taxonomy；七个 rank 列明确留空并记录 `TaxonomyStatus`，本篇不虚构也不使用分类注释
- 主分析从 262,143 个非恒等合法 split-plot assignments 中以 seed 20260723 无放回抽取 9,999 个；四个 simple effects 各穷举 63 个非恒等 paired swaps，并对四项 P 值作 Holm 校正
- 研究锚定 feature filter 下，treatment omnibus total/partial R² 为 0.0981/0.1847、P=0.0065，PERMDISP P=0.3714；三个 factorial terms 与四个 simple effects 经 Holm 校正后均不显著
- 不过滤 feature 的合法置换敏感性分支 P=0.7642，因此正文把结果限定为 filter-sensitive 的条件性整体位置差异，不宣称稳健的普遍处理效应
- 214/214 检查覆盖三件套、技术重复聚合、交换性、置换矩阵、PERMANOVA/PERMDISP、ANOSIM/MRPP 方向、敏感性分支和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/23-permanova-dispersion/`

第 25 篇从 checksum-locked 的真实湿地三件套与环境表独立重算环境方差解构：

- 5% prevalence filter 保留 11,113/13,628 个 OTU；envfit 与 Mantel 使用 relative-abundance Bray–Curtis，VPA 使用 Hellinger-transformed community matrix
- metadata 没有 subject/pair/block/site/time/batch 标识，因此示例用 seed 20260725 生成 999 个唯一自由置换；该边界不作为样本独立性的额外证明
- 11 个连续变量的 envfit 检验经 BH 校正后 10 个达到 0.05；其 r² 只表示与当前两条 PCoA 轴的拟合，不是完整群落解释率
- edaphic、climate、geographic 三个预设 Mantel blocks 的 Spearman r 为 0.332、0.395、0.587，三项 Holm-adjusted P 均为 0.003；Pearson partial Mantel 只保存为非主敏感性结果
- Hellinger VPA combined adjusted R² 为 0.2817；pure edaphic、pure wetland group、shared 和 unexplained fractions 为 0.0967、0.1140、0.0710、0.7183，只有两个 pure fractions 直接检验
- Group 与 Saline 完全混杂，TOC 与 TN 的 Pearson r 为 0.988；正文不把 group fraction 写成盐度因果贡献，也不把 shared fraction 归因给单一变量集
- 221/221 检查覆盖输入、样本对齐、环境变量、共线性、置换、multiplicity、VIF、conditioning/prevalence sensitivity 和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/25-environment-variance/`

## 固定环境

- R 4.4.1
- Quarto 1.9.38
- `env/renv.lock`：147 个 R/CRAN/Bioconductor 包；包括 knitr 1.51、rmarkdown 2.31、tidyverse 2.0.0、vegan 2.6-6.1、phyloseq 1.48.0、microeco 2.0.0、decontam 1.24.0、svglite 2.1.3 和 iNEXT 3.0.2
- `env/qiime2.yml`：QIIME 2 2026.4 官方 Linux 环境文件；本地验证记录见 `env/qiime2-validation.txt`
- `env/qiime2-itsxpress-overlay.yml`：在上述环境固定 ITSxpress 2.1.4 与 Biopython 1.87
- `env/gemelli.yml` 与 `env/gemelli-requirements.txt`：独立固定 Python 3.12.13、gemelli 0.0.13 及其完整依赖，用于 Robust Aitchison/rclr-RPCA
- `env/raw-qc.yml`：独立锁定 Python 3.14.6、OpenJDK 25.0.2、FastQC 0.12.1 与 MultiQC 1.33

恢复 R 环境：

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(lockfile = "env/renv.lock", library = .libPaths()[1], prompt = FALSE)'
```

创建并验证 QIIME 2 环境：

```bash
mamba env create -n microbiome-qiime2-2026.4 \
  -f env/qiime2.yml \
  --channel-priority flexible
conda run -n microbiome-qiime2-2026.4 qiime info

mamba env create -f env/gemelli.yml
conda run -n microbiome-gemelli-0.0.13 python -m pip check
```

## 验证入口

只检查 manifest 与写作契约：

```bash
python ../skills/best-practice-tutorial-style/scripts/render_tutorial.py \
  tutorial.yaml --variant github --check-only

python scripts/validate_series.py \
  --project-root . \
  --manifest tutorial.yaml \
  --output scaffold_validation.json

python3 scripts/validate_wsl2_conda.py \
  --project-root . \
  --output-dir results/07-wsl2-conda \
  --figure-dir figures

python3 scripts/validate_qiime2_install.py \
  --project-root . \
  --output-dir results/08-qiime2-install \
  --figure-dir figures

Rscript --vanilla scripts/validate_r_ecosystem.R \
  --project-root . \
  --output-dir results/09-r-ecosystem \
  --figure-dir figures

python3 scripts/validate_import_provenance_qc.py \
  --project-root . \
  --output-dir results/10-import-provenance-qc \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4 \
  --qc-env microbiome-fastqc-multiqc

python3 scripts/validate_primer_trimming.py \
  --project-root . \
  --output-dir results/11-primer-trimming \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4

python3 scripts/validate_dada2_asv.py \
  --project-root . \
  --input-artifact results/11-primer-trimming/nfcore-v4-trimmed.qza \
  --input-summary results/11-primer-trimming/nfcore-v4-trimmed-summary.qzv \
  --output-dir results/12-dada2-asv \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4

python3 scripts/validate_its_18s.py \
  --project-root . \
  --output-dir results/13-its-18s \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4

python3 scripts/validate_taxonomy_databases.py \
  --project-root . \
  --output-dir results/14-taxonomy-databases \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4

python3 scripts/validate_region_classifier.py \
  --project-root . \
  --output-dir results/15-train-classifier \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4

python3 scripts/validate_sepp_phylogeny.py \
  --project-root . \
  --reference-database downloads/sepp-refs-gg-13-8.qza \
  --output-dir results/16-sepp-phylogeny \
  --figure-dir figures \
  --qiime-env microbiome-qiime2-2026.4 \
  --threads 8

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_phyloseq_import.R \
  --project-root . \
  --input-dir data/small/phyloseq-import \
  --output-dir results/17-phyloseq-import \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_publication_graphics.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/18-publication-graphics \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_alpha_diversity.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/19-alpha-diversity \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_alpha_group_tests.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/20-alpha-group-tests \
  --figure-dir figures

GEMELLI_PYTHON="$(conda run -n microbiome-gemelli-0.0.13 \
  python -c 'import sys; print(sys.executable)')"
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_beta_distances.R \
  --project-root . \
  --input-dir data/small/beta-distances \
  --output-dir results/21-beta-distances \
  --figure-dir figures \
  --gemelli-python "$GEMELLI_PYTHON"

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_permanova_dispersion.R \
  --project-root . \
  --input-dir data/small/permanova-dispersion \
  --output-dir results/23-permanova-dispersion \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_environment_variance.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/25-environment-variance \
  --figure-dir figures
```

执行完整 Pilot QA：

```bash
python ../skills/tutorial-execution-qa/scripts/run_tutorial_qa.py \
  tutorial.yaml \
  --workspace .tutorial_runs/pilot \
  --output qa_report.json
```

直接渲染网站：

```bash
quarto render
```

Docker 构建的是 R + Quarto 网站环境；QIIME 2 上游环境由 `env/qiime2.yml` 单独管理：

```bash
docker build -t microbiome-best-practices:pilot .
docker run --rm -v "$PWD/_site:/project/_site" microbiome-best-practices:pilot
```

## 主要入口

- `tutorial.yaml`：内容、执行与发布契约
- `_quarto.yml`：55 篇 Quarto Book 目录
- `chapters/03-study-design.qmd`：实验单位、混杂审计和效力敏感性样板
- `chapters/04-contamination-controls.qmd`：blank、浓度、批次与 decontam 污染审计
- `chapters/05-batch-effects.qmd`：设计可识别性、群落/feature 批次诊断与敏感性分析
- `chapters/06-data-and-fastq.qmd`：paired-end 原理、FASTQ/Phred、文件同步与三件套 provenance
- `chapters/07-wsl2-conda.qmd`：Windows/WSL2 分层、Miniforge/mamba、环境锁与真实 Linux-side 审计
- `chapters/08-qiime2-install.qmd`：固定环境、真实 EMP import、Artifact 验收、导出回验与报错分层
- `chapters/09-r-ecosystem.qmd`：R/Bioconductor 版本配对、renv 锁定、三件套导入和双对象口径验收
- `chapters/10-import-provenance-qc.qmd`：multiplexed EMP 导入、barcode 拆样、provenance、FastQC/MultiQC 与质控解释
- `chapters/11-primer-trimming.qmd`：引物存在性探针、锚定 q2-cutadapt、成对保留与引物区分支
- `chapters/12-dada2-asv.qmd`：DADA2 截断、expected errors、overlap、逐阶段损失和 ASV 输出审计
- `chapters/13-its-18s.qmd`：ITS/18S 分支、ITSxpress、变长 ASV、UNITE 与 PR2 参考合同
- `chapters/14-taxonomy-databases.qmd`：SILVA/GTDB/Greengenes2 选择、区域参考、统一分类与命名敏感性审计
- `chapters/15-train-classifier.qmd`：区域 reference、known-lineage 留出、Naive Bayes 训练、confidence 敏感性与生产模型
- `chapters/16-sepp-phylogeny.qmd`：GG13_8 reference placement、rooted tree、jplace 支持度、枝长与 feature retention 审计
- `chapters/17-phyloseq-import.qmd`：QIIME 2 标准导出、方向与 ID gate、完整树裁剪、五组件 phyloseq 和 RDS 往返审计
- `chapters/18-publication-graphics.qmd`：色盲友好编码、字体、终稿尺寸与 PDF/SVG/TIFF/PNG 出版导出合同
- `chapters/19-alpha-diversity.qmd`：Alpha 指数、逐样本稀释、等 reads/coverage 标准化、Hill numbers 与四格式出版图
- `chapters/20-alpha-group-tests.qmd`：设计门禁、Kruskal–Wallis/Holm、计划成对比较、Cliff's delta、bootstrap 区间与标准化敏感性
- `chapters/21-beta-distances.qmd`：Bray/Jaccard/UniFrac/Aitchison/Robust Aitchison 的输入、标准化、根树、零值与 PCoA 门禁
- `chapters/22-ordination-unconstrained.qmd`：非约束排序样板
- `chapters/23-permanova-dispersion.qmd`：技术重复聚合、split-plot 交换性、受限置换 PERMANOVA、PERMDISP、计划比较与敏感性分析
- `chapters/24-ordination-constrained-cap.qmd`：约束排序样板
- `R/theme_pub.R`：整仓库共享作图函数
- `data/small/README.md`：Pilot 数据来源与再生说明
- `qa_report.json`：本地发布门禁报告
- `scripts/build_wechat_review_bundle.py`：从通过 QA 的 `_site` 生成 25 篇微信 HTML、优化正文图、确定性封面与草稿 payload
- `rendered/wechat_review_01_25/report.json`：本地生成的微信审阅包结构、尺寸、图片和内容检查报告（默认不进 Git）
- `rendered/wechat_review_01_25/live_report.json`：本地生成的 25 篇官方草稿及封面素材 ID 操作记录；明确记录未发布、未群发（默认不进 Git）

## 许可证

教程文字与原创图使用 CC BY 4.0；本仓库原创代码使用 MIT。第三方数据与软件保持各自许可证。
