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
- 第 01–55 篇均已完成；第 51–55 篇依次覆盖纵向混合模型与 volatility、piecewise SEM、因果中介、孟德尔随机化，以及因果证据与措辞校准。第 51–54 篇从真实数据独立实跑，第 55 篇使用七项一手研究的可追溯证据台账。
- 第 01–29 篇的上一轮 Pilot QA 已通过：固定真实数据及参考资产、55 篇目录、二十九篇 HTML、337 个 PDF/SVG/PNG/TIFF 图形文件和 473 份结果/审计文件均通过，25/25 个 workflow steps 与 809/809 个发布断言全部成功；run key 为 `09613ff88b76b076`。在当前 418-package 锁下，第 30–35 篇的 359/359 项、第 36–40 篇的 406/406 项、第 41–45 篇的 445/445 项、第 46–50 篇的 435/435 项和第 51–55 篇的 413/413 项发布检查均通过；01–55 合并发布门禁为 73/73，run key 为 `2cf0e712b936e213`。
- 经明确授权，第 01–25 篇已从同一份已通过 QA 的 Quarto 输出生成公众号审阅稿：102 张正文图和 25 张封面均已上传，25/25 篇草稿创建成功，草稿总数由 53 增至 78；未调用发布或群发接口。本地可追溯记录见 `rendered/wechat_review_01_25/report.json` 与 `rendered/wechat_review_01_25/live_report.json`（生成目录默认不进 Git）。
- 第 01–55 篇已生成最新本地公众号审阅包，共 55 篇、217 张正文图和 55 张确定性封面；其中第 51–55 篇含 20 张原创图。经明确授权，保留既有第 01–25 篇草稿，并新增第 26–55 篇：115 张正文图和 30 张封面均已上传，30/30 篇草稿创建成功，草稿总数由 78 增至 108；最新列表严格为第 55→26 篇，未调用发布或群发接口。生成记录见 `rendered/wechat_review_01_55/report.json`，在线草稿审计见 `rendered/wechat_review_01_55/live_report.json`（两者默认不进 Git）。
- GitHub 审阅仓库为 `petemeng/microbiome-best-practices`，第 01–55 篇已更新到 Draft PR [#1](https://github.com/petemeng/microbiome-best-practices/pull/1)；该 PR 尚未合并，GitHub Pages 生产部署仍关闭。
- QIIME 2 官方 2026.4 环境、固定 ITSxpress overlay、gemelli 0.0.13、PICRUSt2 2.6.3、原生 MMvec 1.0.5/MOFA2 0.7.4 与 SourceTracker 2.0.1 独立环境均已在本机验证；第 07–17 篇相应环境与上游验收均已通过，第 18–55 篇的出版图形、生态统计、组成数据、差异丰度、定量、网络、功能预测、机器学习、时间结局、跨组学、来源追踪、纵向与因果推断合同均设有独立门禁。

## 已锁定的真实数据

第 01、02、03、09、18、19、22、24、25、26、27、28、30、31、32、33、35–40、42 篇使用 CRAN 归档中的 `microeco 2.0.0`：

- 源码包：`microeco_2.0.0.tar.gz`
- SHA-256：`454a3b71ceeea86bdd54f475a12b5bac9c3b124f663b8dc6e560827fca89ebfd`
- 包许可证：GPL-3
- 表格规模：13,628 个 OTU、90 个 16S 样本、taxonomy、sample metadata 和 200 行环境观测
- 第 38 篇根树：14,096 tips、14,095 internal nodes，包含全部 13,628 个 count-table OTU；固定文件为 `data/small/rooted-tree.nwk.gz`
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

第 26 篇从同一 checksum-locked 湿地三件套独立重算群落组成：

- 门水平 Top 10 由 90 个逐样本闭合后的等权均值确定，平均覆盖 94.80% 的全部 reads；Unknown phylum 与 Other known phyla 分开保存
- Unknown genus 与 species 的等样本均值分别为 63.11% 和 99.63%，因此属水平图保留全部 reads 分母并把 Unknown genus 显式展示
- equal-sample 与 pooled-read 的门水平 Top 10 名单相同，但 coastal wetland 的 weighting total variation 为 2.12%；known-only 重归一化可使单个组–属均值最多增加 7.273 个百分点
- 195/195 检查覆盖输入哈希、七级 taxonomy、closure、Top-N、样本权重、denominator sensitivity 和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/26-community-composition/`

第 27 篇从同一 checksum-locked 湿地三件套独立重算 Phylum、Class、Order、Family、Genus 五级组成：

- 以全部 reads 为分母的等样本平均已知注释比例从 Phylum 的 99.65% 逐级降至 Genus 的 36.89%；最大分辨率损失发生在 Family→Genus，为 31.07 个百分点
- 324 个 feature（17 个已知 genus）存在 Family unknown、Genus known 的缺失父级，占等样本平均 2.52%；流程保留缺口，不用子级名称反向补造父级
- 2 个 Family 原始名称在不同父级下发生碰撞，另有 13 个名称跨 rank 重复；所有展示标签均带 rank，聚合键使用完整 lineage
- 各 rank 的 Top 10 coverage 从 Phylum 的 94.80% 降至 Genus 的 9.52%，因此不同 rank 面板不能把相同 Top-N 当成相同信息覆盖率
- 209/209 检查覆盖输入哈希、五级 closure、注释分辨率、缺失父级、完整 lineage、同名碰撞、Top-N、权重/分母敏感性和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/27-multirank-composition/`

第 28 篇从同一 checksum-locked 湿地三件套独立定义 feature-level 核心菌群与稀有生物圈：

- 主分析以每个样本全部 reads 为分母，固定检测阈值为相对丰度至少 0.01%、prevalence 至少 80%，并按样本等权；得到 97 个全局核心 feature，承载 17.04% 的等样本平均丰度，其中 34 个全局平均丰度仍低于 0.1%
- IW、CW、TW 的组内 core 分别为 166、160、278 个；三组共同 core 为 59 个，三组并集为 425 个，避免把全局 core 与任一组内 core 混为一谈
- detection × prevalence 网格给出 16 组 core 计数；主阈值在当前 10,364–37,374 reads 的文库中相当于 2–4 条 reads。固定 seed 20260728 的 10,000-read 稀释 one-read core 与主定义 Jaccard 为 0.581，因此深度审计只作敏感性证据
- 互斥状态包含 97 个 primary core、200 个 conditionally rare、1,155 个 sampling-limited、8,407 个 persistently rare 和 3,769 个 intermediate feature；样本内低于 0.1% 的 local rare mass 平均为 44.15%，但不被解释为休眠、激活或结构性缺失
- 320/320 检查覆盖输入哈希、occupancy–abundance、组内交集、阈值网格、原始深度与稀释敏感性、状态守恒、local rare mass 和四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/28-core-rare-biosphere/`

第 29 篇使用 Bioconductor `DirichletMultinomial 1.46.0` 随包分发的 Holmes/Turnbaugh Twins 真实属级计数，并导出独立三件套：

- `Twins.csv` 与 `TwinStudy.t` 被确定性整理为 130 个 source genus categories × 278 个样本、570,851 reads；metadata 含 154 个 suffix-derived participant keys，其中 124 个有两份来源样本
- 主模型直接使用未抽平、未加伪计数、未转相对丰度的非负整数 counts；seed `20260729`、候选 `k=1–7`，以 fresh-fit Laplace minimum 选择 `k=4`，AIC/BIC 仅作并列诊断
- 四个 canonical states 按 fitted *Bacteroides*、*Prevotella* 比例确定性命名为 DMM-A 至 DMM-D，样本量为 85、86、46、61；8 个样本的最大 posterior 低于 0.80
- 八个初始化均在“一模型一新 R 进程”中执行，跨 fit profile matching 后 8/8 回到同一 partition，最小 adjusted Rand index 为 1；5% prevalence、测序深度、确定性 1,000-read 抽平与单人单样本分支均重新拟合并保留完整 ledger
- 124 对重复来源样本中 85 对状态相同、39 对不同，agreement 为 0.685；BMIClass 只作拟合后描述，不执行疾病关联、因果或自然类别检验
- 数据研究：[Turnbaugh et al., *Nature*, 2009](https://doi.org/10.1038/nature07540)；DMM 方法与原分析：[Holmes et al., *PLOS ONE*, 2012](https://doi.org/10.1371/journal.pone.0030126)
- 222/222 检查覆盖输入、raw-count DMM、三种模型准则、初始化、posterior、profile label matching、五类敏感性、重复样本与四张图的 PDF/SVG/600 ppi PNG/LZW-TIFF 文件；明细在 `results/29-community-typing-dmm/`

第 30–33 篇继续从真实湿地三件套独立开始：

- 第 30 篇用“只增加一个属”的反事实直接展示 closure，并并列保存 raw reads、relative abundance、CLR 与三种 ALR reference frame；三张原创图均有 PDF/SVG/600 ppi PNG/LZW-TIFF。
- 第 31 篇把 CW 对 IW 的样本、80 行 feature universe、过滤、方向和 BH 门槛锁死，只改变模型；ANCOM-BC2、ALDEx2、MaAsLin2、LEfSe+BH、corncob 分别命中 44、47、52、51、62 个可报告属，union 不解释为 consensus。
- 第 32 篇分别重聚合 Family、Genus、Species；已知 reads 比例为 67.64%、36.96%、0.35%。Family 与 Genus 通过报告门禁；Species 明确失败并保留为注释审计，不输出主 biomarker。
- 第 33 篇让 volcano、审计后的 LEfSe cladogram、Manhattan 与 forest plot 共用同一 ANCOM-BC2 master table；79 个可报告属保持同一 contrast、方向和 FDR 定义。

第 34 篇使用 Bioconductor `reconsi 1.16.0` 分发的 Vandeputte 定量微生物组数据：

- `scripts/prepare_vandeputte_absolute_data.R` 导出 234 个属 × 135 个样本、4,080,996 reads 的三件套和逐样本 flow-cytometry `cell-load.tsv`；来源摘要 SHA-256 为 `8ac35a7119d1772d5d1c92f30c051f2fb3bdb2b5d3dca3f83c885528352d9d94`。
- 主比较只保留同一 Disease cohort 的 29 个 Crohn's disease 与 66 个 healthy 样本；另 40 个 Study healthy 样本不跨 cohort 混并。
- 相对丰度与 quantitative microbiome profiling 使用相同 feature screen，并以 seed `20260734` 的 2,000 次 bootstrap 给出区间；三张原创图均导出四种格式。
- 数据研究：[Vandeputte et al., *Nature*, 2017](https://doi.org/10.1038/nature24460)。

第 35 篇从湿地三件套构建属级 association network：

- 先用与分组无关的 prevalence `>=20%`、total reads `>=100` 规则保留 40 个属，再对 780 个 taxon pairs 推断；seed 固定为 `20260735`。
- SparCC 使用 200 次 bootstrap/permutation、add-one empirical p-value 与 BH，得到 264 条 `|r|>=0.30` 且 q<0.05 的边；SPIEC-EASI MB 使用 20 个 lambda 与 50 次 StARS subsampling，得到 30 条边，instability 为 0.0372。
- 30 条 SPIEC-EASI 边全部属于 SparCC edge；按 wetland group 中心化 CLR 后，156/264 条 SparCC edge 保持同方向且 `|r|>=0.20`。这些边只解释为统计关联，不写成互作或因果。
- 三张原创图均导出 PDF/SVG/600 ppi PNG/LZW-TIFF；SPIEC-EASI 固定到官方仓库提交 `c463727a51d0df34db0c670d3b170195bb3d4eba`。

第 36–40 篇继续从同一真实湿地三件套独立开始：

- 第 36 篇固定同一组 40 个命名属；SparCC 网络有 266 条边、4 个模块和 3 个 connector，SPIEC-EASI 网络有 30 条边。三组 density-matched 网络平均 edge Jaccard 为 0.0528，499 次标签置换 lower-tail P=0.002；鲁棒性只解释为推断图属性。
- 第 37 篇对 422 个命名属执行 CLR-WGCNA；power 7 首次同时通过 signed scale-free R² ≥0.80 与 mean connectivity ≥5。得到 4 个非灰模块、115 个 grey 属、32 项 module–trait tests 中 23 项 BH q<0.05，并以 pseudocount sensitivity（ARI 0.715）审计 19 个 candidate hubs。
- 第 38 篇从同一包的 rooted branch-length tree 开始；主过滤保留 585 个 OTU 与 47.4% reads，用 999 次 null randomization 计算 βNTI/RCbray，并显式固定 parallel worker RNG streams。三个预设 short-range niche-signal tests 中 2 个通过 BH 门禁，另用 220-OTU/199-randomization 分支审计过滤敏感性。
- 第 39 篇对 13,628 个 OTU 拟合 Sloan NCM：pooled migration parameter m=0.1029、abundance–occupancy R²=0.7115；4,531/1,086/8,011 个 OTU 分别落在 envelope 上方/下方/内部，并报告 pool、detection 和 10,000-read rarefaction sensitivity。
- 第 40 篇在 3,208 个过滤 OTU 上计算 standardized Levins sample-use breadth；90 个逐样本 RAD fits 中 Zipf–Mandelbrot 与 Zipf 分别胜出 75 和 15 次。CW/IW/TW 的 within-region distance-decay slopes 均为负，但 TW Mantel P=0.287，因而不把 pooled 空间格局写成单一 dispersal 机制。
- 二十张原创图均导出 PDF/SVG/600 ppi PNG/LZW-TIFF；结果 ledger 位于 `results/36-network-robustness/` 至 `results/40-niche-distance-decay/`。

第 41–45 篇使用三套新增真实数据与两套已锁定三件套：

- 第 41 篇把官方 chemerin 37-ASV × 24-sample 数据真实跑入 PICRUSt2 2.6.3；37/37 个 ASV 通过 NSTI 2.0 门槛，maximum NSTI 为 0.2175，sample-weighted median NSTI 为 0.0423。225 个 MetaCyc pathways 的 facility-adjusted genotype tests 完整做 BH；没有 q<0.05 pathway，因此正文不虚构显著通路。
- 第 42 篇从 microeco 湿地三件套重跑 FAPROTAX 1.2.12；3,992/13,628 个 features 被映射，覆盖 34.75% reads。genus-only coverage 仅 0.478%，作为 taxonomy-resolution sensitivity 的实证边界；93 个 functions 在同一 FDR family 中检验并报告 epsilon-squared。
- 第 43 篇在 Vandeputte Disease cohort 上执行 3 × 5 nested CV；classification sample-level AUC 为 0.9948（95% CI 0.9843–1.0000）、Brier 0.0361、50-run label-permutation p=0.0196。microbial-load regression 的 nested-CV R² 为 0.4358，而 apparent training R² 为 0.9160，直接展示过拟合差距。
- 第 44 篇锁定 MicrobiomeHD Xiang/Zhao/Zackular 三个 CRC 16S cohorts；Peptostreptococcus 的 REML pooled CLR effect 为 0.6007、BH p=0.000118。leave-one-study-out AUC 分别为 0.588、0.436、0.578，macro mean 仅 0.534，因此明确报告 transportability failure，而不只展示内部高分。
- 第 45 篇锁定 MiSurv 173 只 NOD mice 的 baseline 16S 与 T1D onset/censoring；118 events、55 censored，每只动物从实际 sampling week 6–8 计算随访。逐属 adjusted Cox 没有 BH-significant genus；training/test 为 120/53，locked-test C-index 0.538，10/15/20-week AUC 为 0.489/0.611/0.603。
- 二十张原创图均导出 PDF/SVG/600 ppi PNG/LZW-TIFF；数据来源和 checksum 见 `data/small/picrust2-chemerin/`、`cross-cohort-crc/`、`survival-t1d/`，结果 ledger 位于 `results/41-picrust2/` 至 `results/45-survival-analysis/`。

第 46–48 篇使用 Franzosa IBD 队列的配对微生物组与代谢组数据：

- 数据由 Muller 等整理的公开仓库固定到提交 `89a519d8c832008fbc6e650453e83e2f04858d02`；220 个配对样本包含 56 个 Control、88 个 CD 和 76 个 UC，整理为 16S 三件套、277 个高置信命名代谢物及其注释表。原始文件、过滤规则和 SHA-256 见 `data/small/paired-ibd-multiomics/source-summary.json`。
- 第 46 篇在相同样本上完成 Procrustes、Mantel、逐对 Spearman 与 HAllA；组别残差化后，25 × 40 个特征中有 350 个 FDR 显著 pair，归入 70 个 HAllA blocks。全局 Procrustes/Mantel 统计量分别为 0.5288/0.5543，受限置换 q 值均为 0.001。
- 第 47 篇使用原生 mixOmics 6.26.0，在固定的 153/67 train/test split 上拟合 sPLS 与 DIABLO；test balanced accuracy 为 0.6758，macro AUC 为 0.7984，26 个特征达到至少 80% 的重复选择率。
- 第 48 篇在隔离环境中运行原生 MMvec 1.0.5（TensorFlow 1.15）和 mofapy2 0.7.4。MMvec 使用 153 个训练样本与 67 个 holdout 样本，100 epochs 的最佳边界 epoch 为 100；MOFA 保留 4 个 factors，两组学总方差解释率分别为 34.83% 和 32.59%。
- 数据研究：[Franzosa et al., *Nature Microbiology*, 2019](https://doi.org/10.1038/s41564-018-0306-4)；整理资源：[Muller et al., *Scientific Data*, 2022](https://doi.org/10.1038/s41597-022-01644-2)。

第 49 篇使用 Duran 等拟南芥根部多界微生物组数据：

- 官方补充仓库固定到提交 `6db5e85cc5d442fd95fcdcb7250b72fa9e2ff900`；36 个严格配对样本覆盖 3 个 soil × 3 个 compartment × 4 个 biological replicates，包含 772 个 bacterial、1,063 个 fungal 和 219 个 oomycete features。
- 调整 soil 与 compartment 后，1,200 个跨界 pair 中仅 2 个通过 BH；491 个 raw association 在调整后改变方向。三组 Procrustes 比较中 2 组通过调整后的 FDR 门槛，30 个 bootstrap candidates 中有 2 个保留稳定证据。
- 数据研究：[Duran et al., *Cell*, 2018](https://doi.org/10.1016/j.cell.2018.10.020)；来源、配对键和 SHA-256 见 `data/small/multi-kingdom-duran/source-summary.json`。

第 50 篇使用 FEAST 官方来源追踪示例：

- 官方仓库固定到提交 `2f8f3df8051e0e08341f597a9f4693bfb76b3bf6`；数据含 1 个 sink、9 个 candidate-source samples、4 个 source classes 和 1,839 个匿名 features，已整理为三件套并保留 source/sink 标签。
- FEAST 与 SourceTracker 都把 Infant gut 识别为主来源，但 unknown 分别为 18.54% 和 0.44%；两种方法的 leave-one-source-out known-class accuracy 均为 0.50。缺失来源敏感性显示 unknown 不保证单调上升，因此正文不把单次比例当作传播证明。
- 数据研究与方法：[Shenhav et al., *Nature Methods*, 2019](https://doi.org/10.1038/s41592-019-0431-x)；来源和 SHA-256 见 `data/small/source-tracking-feast/source-summary.json`。

第 46–50 篇共二十张原创主图，均导出 PDF/SVG/600 ppi PNG/LZW-TIFF；435/435 项独立发布检查通过，结果 ledger 位于 `results/46-integration/` 至 `results/50-source-tracking/`。

第 51 篇使用 O'Keefe 等人的 DietSwap 真实纵向数据：

- Bioconductor `microbiome 1.26.0` 的 `dietswap` 对象被整理为 130 features × 222 observations，来自 38 位受试者和 6 个 timepoints；这是 HITChip 16S 系统发育芯片数据，适合演示重复测量模型，但其 assay-specific preprocessing 不等同扩增子测序。
- random-intercept 与 AR(1) 模型的 AIC 分别为 326.2578 和 326.2104；AR(1) `phi=0.1383`，likelihood-ratio `P=0.1525`。相邻时点共有 180 对，cohort-by-time 的 T4 interaction estimate 为 0.5288、`P=0.0112`。
- 数据研究：[O'Keefe et al., *Nature Communications*, 2015](https://doi.org/10.1038/ncomms7342)；原始记录：[Dryad](https://doi.org/10.5061/dryad.1mn1n)；来源与 checksum 见 `data/small/longitudinal-dietswap/source-summary.json`。

第 52 篇从真实湿地三件套与实测环境表建立预设 SEM：

- 90 个样本同时进入 piecewise SEM 与 PLS-PM；完整预设 DAG 的 Fisher's C 为 4.514、`P=0.105`，删去直接环境路径的竞争模型 Fisher's C 为 51.796、`P<0.001`。
- Shannon→TOC 的标准化路径为 -0.0123、`P=0.8647`，PLS-PM goodness-of-fit 为 0.4597；正文保留空结果并把 directed separation、竞争图和残差诊断作为主审计。

第 53 篇使用 Franzosa IBD 配对微生物组—代谢组数据演示中介分析：

- 预先指定 steroid exposure→*Faecalibacterium* CLR→胆汁酸 outcome 的 108 个完整样本中，28 个 exposed、80 个 unexposed。
- 2,000 次模拟得到 ACME 0.0102（95% CI -0.1046–0.1370，`P=0.846`）；横断面时间顺序不足，因此结果只作为观测中介示范，不写成机制证明。

第 54 篇使用 MiBioGen 与公开 IBD GWAS 汇总统计：

- 暴露为 `genus.Bifidobacterium.id.436`，结局为 GCST004131（25,042 cases、34,915 controls）；6 个跨染色体候选工具经 harmonisation 删除回文 SNP `rs10841473` 后保留 5 个，minimum F statistic 为 20.7695。
- IVW beta 为 -0.0167、`P=0.8756`、OR 0.9835；heterogeneity `P=0.1063`、MR-Egger intercept `P=0.8254`、MR-PRESSO global `P=0.183`。当前冻结数据不支持 Steiger 与 colocalisation，正文明确标记不可得，而非补造诊断。
- 暴露研究：[Kurilshikov et al., *Nature Genetics*, 2021](https://doi.org/10.1038/s41588-020-00763-1)；结局研究：[de Lange et al., *Nature Genetics*, 2017](https://doi.org/10.1038/ng.3760)；来源与 API ledger 见 `data/small/mr-mibiogen/`。

第 55 篇使用七项一手研究构建从横断面关联到机制链与三角验证的阅读框架；该篇不执行统计模型，也不把证据阶梯当自动评分。结构化 DOI、可支持结论与残余威胁见 `data/small/causal-evidence/evidence-cases.tsv`。

第 51–55 篇共二十张原创主图，均导出 PDF/SVG/600 ppi PNG/LZW-TIFF；413/413 项独立发布检查通过，结果 ledger 位于 `results/51-longitudinal-analysis/` 至 `results/55-causal-evidence/`。

## 单篇复现入口

计算章节先给出真实数据下载与读取，再进入分析。第 02–06、17–40、42–54 篇的 `examples/` 脚本与网页版的完整 R 工作流逐块保持一致；装好所需依赖后，在新建文件夹中按顺序运行，不需要先生成其他章节的中间对象。

第 17、18、21、46–50 篇还会自动下载公开、固定版本的计算脚本；涉及 Python 的章节在正文中列出独立环境文件和解释器选择。第 07–16、41 篇保留完整命令与一次性执行记录，不能把环境检查或重绘既有输出等同于重新运行重型上游。第 01、55 篇以数据结构、真实图表和证据判读为主，不强制运行统计代码。

微信保留数据入口、关键分析代码和完整脚本链接，省略通用绘图函数定义；完整 R 工作流仍逐块来自网页版。维护时先运行 `scripts/validate_reader_series.py` 的空目录检查，再验证微信中的分析代码未被改写或打乱，省略范围仅限通用绘图脚手架。方法学与数值检查仍使用各章原有验证程序，空目录检查不替代它们。

## 固定环境

- R 4.4.1
- Quarto 1.9.38
- `env/renv.lock`：418 个 R/CRAN/Bioconductor/GitHub 包；新增并固定 piecewiseSEM 2.3.0.1、plspm 0.6.0、mediation 4.5.1、emmeans 1.10.2、TwoSampleMR 0.6.6、MRPRESSO 1.0 及其传递依赖，同时保留 mixOmics、FEAST、ggpicrust2、ranger、pROC、metafor、timeROC、WGCNA、iCAMP 与既有差异丰度/网络依赖
- `env/qiime2.yml`：QIIME 2 2026.4 官方 Linux 环境文件；本地验证记录见 `env/qiime2-validation.txt`
- `env/qiime2-itsxpress-overlay.yml`：在上述环境固定 ITSxpress 2.1.4 与 Biopython 1.87
- `env/gemelli.yml` 与 `env/gemelli-requirements.txt`：独立固定 Python 3.12.13、gemelli 0.0.13 及其完整依赖，用于 Robust Aitchison/rclr-RPCA
- `env/picrust2.yml`：独立固定 Python 3.12.13 与 PICRUSt2 2.6.3，用于第 41 篇真实功能预测
- `env/mmvec-native.yml`：隔离固定 Python 3.7.16、TensorFlow 1.15.0 与 MMvec 1.0.5，用于第 48 篇原生条件共现模型
- `env/multiomics.yml`：固定 Python 3.10.19、HAllA 0.8.40 与 mofapy2 0.7.4，用于第 46、48 篇原生跨组学分析
- `env/sourcetracker2.yml`：固定 Python 3.10.20 与 SourceTracker 2.0.1，用于第 50 篇来源追踪对照
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

mamba env create -f env/picrust2.yml
conda run -n picrust2-2.6.3 picrust2_pipeline.py --version

CONDA_CHANNEL_PRIORITY=flexible mamba env create -f env/mmvec-native.yml
conda run -n mmvec-native python -c 'import mmvec, tensorflow as tf; print(mmvec.__version__, tf.__version__)'

mamba env create -f env/multiomics.yml
R_LIBS_USER="$PWD/.r-lib" conda run -n multiomics-native python -c 'from importlib.metadata import version; from halla import HAllA; from mofapy2.run.entry_point import entry_point; print(version("halla"), version("mofapy2"))'

mamba env create -f env/sourcetracker2.yml
conda run -n sourcetracker2 sourcetracker2 --help
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

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_community_composition.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/26-community-composition \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_multirank_composition.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/27-multirank-composition \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_core_rare_biosphere.R \
  --project-root . \
  --input-dir data/small \
  --output-dir results/28-core-rare-biosphere \
  --figure-dir figures

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla \
  scripts/validate_community_typing_dmm.R \
  --project-root . \
  --input-dir data/small/community-typing-dmm \
  --output-dir results/29-community-typing-dmm \
  --figure-dir figures
```

验证第 30–35 篇正文、HTML、图形、数据和环境锁：

```bash
python3 scripts/validate_articles_30_35.py \
  --project-root . \
  --site-root . \
  --output results/articles-30-35-validation.json
```

验证第 36–40 篇正文、HTML、20 张主图四格式、根树、结果和环境锁：

```bash
python3 scripts/validate_articles_36_40.py \
  --project-root . \
  --site-root . \
  --output results/articles-36-40-validation.json
```

验证第 41–45 篇正文、HTML、20 张主图四格式、三套真实数据、结果和环境锁：

```bash
python3 scripts/validate_articles_41_45.py \
  --project-root . \
  --site-root . \
  --output results/articles-41-45-validation.json
```

验证第 46–50 篇正文、HTML、20 张主图四格式、三套真实数据、原生环境与结果：

```bash
python3 scripts/validate_articles_46_50.py \
  --project-root . \
  --site-root . \
  --output results/articles-46-50-validation.json
```

验证第 51–55 篇正文、HTML、20 张主图四格式、纵向/SEM/中介/MR 数据与因果证据台账：

```bash
python3 scripts/validate_articles_51_55.py \
  --project-root . \
  --site-root . \
  --output results/articles-51-55-validation.json

python3 scripts/validate_release_01_55.py \
  --project-root . \
  --output results/articles-01-55-release-validation.json
```

执行完整 Pilot QA：

```bash
python ../skills/tutorial-execution-qa/scripts/run_tutorial_qa.py \
  tutorial.yaml \
  --workspace .tutorial_runs/pilot \
  --output qa_report.json
```

从已通过 QA 的站点生成 01–55 篇本地公众号审阅包（不调用微信接口）：

```bash
python3 scripts/build_wechat_review_bundle.py \
  --site-dir _site \
  --qa-report results/articles-01-55-release-validation.json \
  --formal-count 55 \
  --output-dir rendered/wechat_review_01_55
```

第 26 篇的可下载脚本位于 `examples/26-community-composition.R`。正文中的 `reader-*` / `fig-*` R 代码、该脚本和微信正文必须保持相同顺序和内容；生成微信包时会自动检查，必要的数据加载代码被移除会直接报错。

发布前还需从实际微信 payload 提取代码，在**不存在的新目录**和 `Rscript --vanilla` 会话中运行。先安装正文声明的 R 包，再执行：

```bash
python3 scripts/validate_reader_composition.py \
  --project-root . \
  --draft-json rendered/wechat_review_01_55/26/draft.json \
  --work-dir /tmp/ch26-reader-clean-run \
  --output /tmp/ch26-reader-validation.json
```

此检查不复制项目数据：代码自行下载三张表，核对门/属相对丰度与独立参考结果，检查四张图的四种导出格式、600 ppi 栅格尺寸、英文标签及其与网站导出 PNG 的一致性。它验证的是空白工作目录中的复现，不是未安装任何软件的操作系统；也不替代独立的统计结果检查。

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
- `chapters/25-environment-variance.qmd`：envfit、Mantel、VPA、共线性、置换与解释边界
- `chapters/26-community-composition.qmd`：堆叠、气泡、冲积、热图、Top-N 与分母/样本权重审计
- `chapters/27-multirank-composition.qmd`：Phylum–Genus 分辨率、缺失父级、lineage、同名碰撞与跨 rank Top-N 审计
- `chapters/28-core-rare-biosphere.qmd`：detection/prevalence 核心定义、组内交集、深度/阈值敏感性与互斥稀有状态
- `chapters/29-community-typing-dmm.qmd`：raw-count DMM、component 数、posterior uncertainty、初始化/预处理敏感性与重复样本稳定性审计
- `chapters/30-compositional-data.qmd`：closure、raw/relative/CLR 与 reference-frame 反事实审计
- `chapters/31-da-methods.qmd`：ANCOM-BC2/ALDEx2/MaAsLin2/LEfSe/corncob 同设计横评
- `chapters/32-multirank-da.qmd`：Family/Genus/Species 独立聚合、检验族与注释报告门禁
- `chapters/33-da-visualization.qmd`：同一主结果表驱动 volcano、cladogram、Manhattan 与 forest plot
- `chapters/34-absolute-quantification.qmd`：Vandeputte flow-cytometry QMP、相对/定量尺度与 cohort 门禁
- `chapters/35-cooccurrence-networks.qmd`：SparCC、SPIEC-EASI、StARS、方法交集与组驱动边审计
- `chapters/36-network-robustness.qmd`：Zi–Pi、candidate keystone、节点移除鲁棒性与 density-matched 组间 rewiring
- `chapters/37-microbial-wgcna.qmd`：CLR-WGCNA、soft-threshold 双门禁、模块–性状关联与 hub sensitivity
- `chapters/38-community-assembly-bnti.qmd`：根树合同、phylogenetic-signal gate、βNTI/RCbray 与 deterministic parallel nulls
- `chapters/39-neutral-community-model.qmd`：Sloan abundance–occupancy、migration parameter、predictive envelope 与 detection/depth sensitivity
- `chapters/40-niche-distance-decay.qmd`：Levins sample-use breadth、逐样本 RAD 与分层 geographic distance-decay
- `chapters/41-picrust2.qmd`：PICRUSt2 2.6.3 实跑、NSTI/coverage 审计与 ggpicrust2 MetaCyc 图
- `chapters/42-functional-guilds.qmd`：FAPROTAX coverage、全族 BH/效应量与 taxonomy-resolution sensitivity
- `chapters/43-random-forest.qmd`：classification/regression nested CV、calibration、label permutation 与稳定 importance
- `chapters/44-cross-cohort-validation.qmd`：三队列 REML meta、异质性审计与 leave-one-study-out 外部验证
- `chapters/45-survival-analysis.qmd`：逐样本 time origin、adjusted Cox/PH audit、ridge Cox 与 locked-test time-dependent ROC
- `chapters/46-procrustes-mantel-halla.qmd`：Procrustes/Mantel 全局一致性、组别调整与 HAllA 层级关联
- `chapters/47-spls-diablo.qmd`：原生 mixOmics sPLS/DIABLO、锁定测试集、调参与选择稳定性
- `chapters/48-mmvec-mofa.qmd`：原生 MMvec 条件共现、MOFA2 潜因子与边界/holdout 审计
- `chapters/49-multi-kingdom.qmd`：16S＋ITS/18S 配对设计、跨界 Procrustes 与协变量调整
- `chapters/50-source-tracking.qmd`：FEAST/SourceTracker、unknown 来源、leave-one-source-out 与缺失来源敏感性
- `chapters/51-longitudinal-analysis.qmd`：subject-level 混合模型、AR(1)、群落轨迹与相邻时点 volatility
- `chapters/52-structural-equation-model.qmd`：预设 DAG、piecewise SEM、directed separation、竞争图与 PLS-PM 审计
- `chapters/53-mediation-analysis.qmd`：ACME/ADE、bootstrap/simulation 区间与未测 mediator–outcome 混杂敏感性
- `chapters/54-mendelian-randomization.qmd`：MiBioGen/GWAS harmonisation、IVW/weighted median/MR-Egger、leave-one-out 与 pleiotropy 审计
- `chapters/55-causal-evidence.qmd`：七级证据阅读框架、残余威胁矩阵、可辩护措辞与三角验证
- `scripts/prepare_vandeputte_absolute_data.R`：从 reconsi 1.16.0 真实数据重建第 34 篇三件套与 cell-load 表
- `scripts/validate_articles_30_35.py`：第 30–35 篇正文、HTML、19 张主图四格式、真实数据与环境锁的独立发布门禁
- `results/articles-30-35-validation.json`：第 30–35 篇 359 项发布检查结果
- `scripts/validate_articles_36_40.py`：第 36–40 篇正文、HTML、20 张主图四格式、真实根树、环境锁与关键结果的独立发布门禁
- `results/articles-36-40-validation.json`：第 36–40 篇发布检查结果
- `scripts/validate_articles_41_45.py`：第 41–45 篇正文、HTML、20 张主图四格式、三套真实数据、PICRUSt2/R 环境与关键结果的独立发布门禁
- `results/articles-41-45-validation.json`：第 41–45 篇 445 项发布检查结果
- `scripts/validate_articles_46_50.py`：第 46–50 篇正文、HTML、20 张主图四格式、三套真实数据、原生环境与关键结果的独立发布门禁
- `results/articles-46-50-validation.json`：第 46–50 篇 435 项发布检查结果
- `scripts/validate_articles_51_55.py`：第 51–55 篇正文、HTML、20 张主图四格式、四类分析输入、418 包环境与关键结果的独立发布门禁
- `results/articles-51-55-validation.json`：第 51–55 篇 413 项发布检查结果
- `scripts/validate_release_01_55.py`：组合既有 Pilot QA、系列合同和第 30–55 篇专项门禁
- `results/articles-01-55-release-validation.json`：01–55 合并发布门禁与 run key
- `scripts/validate_multirank_composition.R`：第 27 篇 209 项确定性数据、层级、敏感性与图形验收
- `results/27-multirank-composition/`：第 27 篇完整 lineage 映射、五级组成、审计表、摘要与日志
- `scripts/validate_core_rare_biosphere.R`：第 28 篇 320 项确定性定义、阈值、深度、状态与图形验收
- `results/28-core-rare-biosphere/`：第 28 篇 occupancy、组内 core、敏感性、状态、摘要与日志
- `scripts/validate_community_typing_dmm.R`：第 29 篇 222 项输入、模型选择、posterior、稳定性、敏感性与图形验收
- `results/29-community-typing-dmm/`：第 29 篇完整模型、component profile、样本 posterior、稳定性 ledger、摘要与日志
- `R/theme_pub.R`：整仓库共享作图函数
- `data/small/README.md`：Pilot 数据来源与再生说明
- `qa_report.json`：本地发布门禁报告
- `scripts/build_wechat_review_bundle.py`：从通过 QA 的 `_site` 默认生成 55 篇微信 HTML、优化正文图、按术语安全换行的确定性封面与草稿 payload
- `rendered/wechat_review_01_55/report.json`：本地生成的 55 篇微信审阅包结构、尺寸、217 张正文图和内容检查报告（默认不进 Git）
- `rendered/wechat_review_01_25/live_report.json`：此前已授权创建的 25 篇官方草稿及封面素材 ID 操作记录；明确记录未发布、未群发（默认不进 Git）

## 许可证

教程文字与原创图使用 CC BY 4.0；本仓库原创代码使用 MIT。第三方数据与软件保持各自许可证。
