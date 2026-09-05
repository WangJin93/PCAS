# PCAS R 包代码审查与优化建议（仅分析，未改动任何代码）

- 审查对象：`WangJin93/PCAS`（2024，IJMS 25:6690 配套工具），本地路径 `/home/Jingle/data/PCAS`，commit `1ad16bf`
- 审查方式：通读 `R/`、`DESCRIPTION`、`NAMESPACE`、`man/`、`README.md` 及 `inst/shinyapp/` 主要模块；对真实数据 API（jingege.wang）与本地 `.rda` 做了实测验证
- 结论先行：**所有下述"问题"均经代码阅读 + 本地 R 实验或真实 API 探测确认**；建议按文末"改造优先级"分三批实施。本文档不改代码，仅作为实施蓝图。

---

## 1. 包总体架构与数据流

```
get_data()                      ← 唯一 HTTP 入口（jsonlite 直连 API）
   │
   ├── get_expr_data()          ← 表达/磷酸化取数，宽表→样本×基因长表，本地 .RData 缓存
   │       │
   │       ├── merge_clinic_data()    ← 合并临床表（仅 Tumor、内连接、分期简化）
   │       ├── cor_cancer_genelist()  ← 单队列内 基因×基因/位点 相关
   │       ├── cor_pancancer_genelist() ← 泛癌 目标基因×基因集 相关
   │       ├── cor_pancancer_drug()   ← 泛癌 目标基因×药物敏感性 相关
   │       ├── cor_pancancer_TIL()    ← 泛癌 目标基因×免疫浸润 相关
   │       └── viz_TvsN()             ← Tumor vs Normal 箱线/小提琴
   │
   ├── get_DEGs_result()        ← 差异分析结果（limma/t.test）→ viz_DEGs_volcano()
   │
   └── viz_cor_heatmap() / viz_corplot() / viz_phoso_sites()
                     ↑
              PCAS_app()  →  inst/shinyapp（Shiny 前端，全部函数在此被调用）
```

内置数据（`data/`，LazyData=TRUE，共 10 个对象）：

| 对象 | 维度 | 用途 | 备注（实测） |
|---|---|---|---|
| `dataset_info` | 56×8 | 数据集目录（Abbre/Normal/Tumor 计数等） | **19 个数据集 Normal=NA（无正常组织样本）**，含 LUAD_APOLLO、GBM_Pediatric、PDAC_KU、OV_PTRC*、BRCA_PTRC 等 |
| `idmap_RNA` | 60774×4 | mRNA：Symbol→探针/转录本 id | 119 个 Symbol 对应 >1 个转录本 |
| `idmap_protein` | 172404×2 | 蛋白/磷酸化位点 id（`NP_…:s315`、`NP_…:s218y223t227` 组合式） | 一个 Symbol 可对应多行（蛋白 + 位点） |
| `drug_CPTAC` | 1553×198 | 药物敏感性矩阵（行名=样本 ID） | 仅覆盖部分队列样本 |
| `drug_info` | 198×6 | 药物-通路注释（Target.pathway） | 无 NA |
| `TIL_CPTAC` | 1636×138 | 免疫浸润（多算法列，列名带算法后缀） | 列名如 `Bcells_EPIC` |
| `TIL_map` | 137×2 | cell_type ↔ algorithm | 过滤 TIL 列的依据 |
| `uniport_map` | 15162×4 | Symbol→UniProt Entry（Reviewed 标志） | |
| `ID_list_pro` | list | 前端基因下拉 | |
| `phoso_id` | list | 按数据集列出的可用磷酸化位点 | |

实测 API 返回形态：
- `action=expression`：每基因 1 行，列 = `row_names`(基因/位点 id) + 各样本列（`C3L-00094_Tumor`/`…_Normal`，另有 `Normal Only IR_Other`、`Taiwanese IR_Other`、`Tumor Only IR_Other` 三列样本池）。
- `action=DEGs`：ttest 表 5 列（row_names/P.Value/logFC/adj.P.Val/change）；limma 表 8 列（额外 AveExpr/t/B）。
- `action=clinic`：19 列；实测 LUAD_APOLLO 中 `Days_To_Recurrence` 66/101 为 NA、`Pack_Years_Smoked` 34/101 为 NA。

---

## 2. 逐函数解读与问题清单

### 2.1 `get_data()`（R/get_data.R）——底层 API 查询

**功能解读**：包内唯一网络入口。把 `table`/`action`/`genes` 拼进固定 URL，`jsonlite::fromJSON` 解析后原样返回（列名随 API 变）。

**问题**：
1. 无参数校验：`action` 任意字符串都发请求；`table` 不存在时不拦截（应先用 `dataset_info$Abbre` 校验）。
2. **无错误处理**：服务器 5xx/断网/非 JSON 返回时 `fromJSON` 抛原生解析错误，用户看到的是不可读报错；无超时、无重试。
3. **URL 不编码**：磷酸化位点 id 含 `:`，若未来含 `+`/`/`/`&` 的符号/别名会破坏查询；多个基因用 `paste0(collapse=",")` 后整体拼 URL，分隔符语义依赖服务端。
4. 返回结构不确定（可能 data.frame 也可能 list），调用方无契约可言。
5. 文档错误：roxygen 写 `@param gene` 而形参是 `genes`（README 同样混乱）。

**缺失数据相关**：函数无法区分"查询成功但无该基因"与"服务器错误/无响应"——两者都被压成空/报错，是全部下游"静默缺失"的源头。

**建议**：改用 `httr::GET`（或 `curl`）+ 显式 timeout + 检查 HTTP 状态 + `tryCatch` 解析；把结果统一为 `list(status, data, message)` 或约定"空 data.frame = 无数据，NULL+warning = 网络/服务错误"；`match.arg` 校验 action；对每个 `genes` 元素 `URLencode`；把 API 域名/超时做成 `options()`，便于部署切换。

---

### 2.2 `get_expr_data()`（R/get_expr_data.R）——表达数据取数（缺失数据问题最集中）

**功能解读**：对每个 dataset：mRNA 经 `idmap_RNA` 把 Symbol 映射成探针 id，蛋白/磷酸化直接用传入 id；调 `get_data` 后转置为"样本×基因"，解析出 `ID`/`type` 列（按样本名最后一个 `_` 分词），`plyr::rbind.fill` 跨数据集合并，按 `dataset_基因md5.RData` 缓存到工作目录 `data_files/data_temp/`；最后重排为 `ID,type,dataset,<基因列>` 并强制数值。

**问题（缺失数据反馈为主）**：
1. **部分基因缺失完全静默（实测复现）**：请求 `c("TP53","ZZZNOTAGENE")`，返回 214×4 只有 TP53，全程无任何"ZZZNOTAGENE 未找到"提示；最终 `intersect(colnames, genes)` 把缺失基因列静默删掉，调用方无从得知丢了什么。
2. **mRNA 符号不在 idmap 时**：`ids <- character(0)`，带着空查询请求 API（代码无前置校验）；仅当 API 恰好返回 0 行时 console 才有 `retrive no results` 消息（拼写错误），且该消息在 Shiny 中用户不可见。
3. **缓存写当前工作目录**：`dir.create("data_files")` 失败（只读目录/shinyapps.io 部署）直接报错中断；缓存键不含时间/版本，数据库更新后永不失效；无开关可关。
4. 缓存命中判断用 `nrow(data)`/`is.null(nrow(...))` 双路径，0 行与 NULL 行为不一致。
5. `type` 靠"最后一个 `_` 后 token"推断：对 `Normal Only IR_Other` 等样本池列可解析出 `Other`，但遇到未来其它命名将静默产生未知 type，后续 `filter(type %in% c("Tumor","Normal"))` 把它们丢弃且无统计。
6. 同一 Symbol 有多个转录本（119 个 Symbol）时 mRNA 合并会**产生同一样本多行**，下游 t.test/相关分析把样本重复计数而无任何提示。
7. 单基因分支硬编码 `cptac_data[,3]`（列号假设），多基因分支靠 `intersect` 排序——**返回列顺序是 idmap 顺序而非用户输入顺序**，且单/多基因两套逻辑不一致。
8. `rbind.fill` 使"某数据集缺某基因"变成整列 NA，静默传导到所有 cor/viz 函数。

**建议**：
- 循环内逐数据集、逐基因统计"查到/未查到/全 NA"，结束时统一打印摘要，并在返回值上挂 `attr(x, "availability")`（`dataset × gene` 的存在性/有效样本数矩阵），Shiny 侧可直接渲染成表。
- 请求前先本地校验：mRNA 用 `idmap_RNA` 预检，未命中即 `message("Symbol not in idmap_RNA: …")` 并跳过，不再发空查询。
- 多转录本策略显式化：默认取该 Symbol 首个（或按 `gene_type` 取 protein_coding）转录本，并提示被合并的转录本数。
- 缓存改到 `tempdir()`/`tools::R_user_dir()` 或 `getOption("PCAS.cache.dir")`，键中加入版本/日期，支持 `use_cache=FALSE`；写缓存失败降级为 warning 而不是中断。
- 统一返回列序 = 用户输入基因顺序；单/多基因合并为同一条路径；消息文案统一、修 `retrive`→`retrieve` 等拼写。
- 对 `dataset_info$Normal == NA` 的 19 个数据集，在用户同时要求 Tumor+Normal 时提前提示"该数据集无正常样本"。

---

### 2.3 `get_DEGs_result()`（R/get_DEGs_result.R）——差异分析结果取数

**功能解读**：按 `method` 拼表名（limma→`<dataset>_limma`，其它→`<dataset>_ttest`），调 `action=DEGs`；mRNA 结果把首列改名 `mRNAs` 并与 `idmap_RNA` 合并补 Symbol；logFC/P.Value 转数值；结果缓存到 `data_files/DEG_results/`。

**问题**：
1. `method` 不校验：`method="wilcox"` 会静默落到 ttest 分支（`ifelse` 逻辑），用户以为算了 wilcox。
2. `dataset` 不存在/API 返回空 → `colnames(results)` 直接报原生错误；无 tryCatch、无"该数据集不存在，可用 dataset_info$Abbre 查看"式反馈。
3. mRNA 分支 `merge(idmap_RNA…)` 是**内连接**：无法匹配到 Symbol 的 mRNA 行被静默丢弃，无前后行数对比。
4. 缓存同 get_expr_data（工作目录、不过期）。
5. 依赖 `str_detect(dataset,"mRNA")` 判断字段处理，命名约定脆弱。

**建议**：`method` 用 `match.arg(c("t.test","limma"))`；请求前校验 dataset；空结果/网络错误给明确 message 并返回 NULL；mRNA 合并前后打印行数及未匹配数；返回对象附加 `attr(,"method")`/`attr(,"dataset")` 供下游校验；缓存改造同上。

---

### 2.4 `merge_clinic_data()`（R/merge_clinic_data.R）——临床合并（"缺失即静默"的典型）

**功能解读**：拉 cohort 的 clinic 表并删首列，把 `data_input` 过滤为 Tumor 后与临床表按 `ID = Cases_Submitter_ID` **内连接**；若存在分期列则生成去掉 A/B/C（或 a/b/c）的 `*_simplify` 列。

**问题（缺失数据最典型）**：
1. `clinic` 取数失败（NULL）→ `clinic[-1]` 得 NULL → `merge()` 抛原生错误；`data_input=NULL` → `dplyr::filter` 报错——都无友好提示。
2. **内连接静默丢弃无临床数据的样本**：不报告"表达样本 X 个、匹配临床 Y 个、丢弃 Z 个（无临床记录）"。实测 1553 例药物矩阵、1636 例 TIL 矩阵本就只覆盖部分样本，同类问题贯穿全部合并。
3. 过滤 `type=="Tumor"` 后可能为 0 行，仍继续 merge，返回空表无提示。
4. 各临床字段 NA 不统计（实测 APOLLO 中 Days_To_Recurrence 缺失 65%）：用户画图/建模时才发现。
5. `AJCC_*`/`Tumor_Stage` 简化列**只在列存在时**生成，不同 cohort 字段不一致，用户事后才发现少列。
6. 返回纯 data.frame，无法携带以上任何诊断信息。

**建议**：入口校验 `data_input` 必需列（ID/type/dataset）与 clinic 成功性；合并前后统计匹配/丢失样本；对结果每个临床字段给出缺失计数（或缺失率>阈值的告警）；`*_simplify` 是否生成写入摘要。推荐返回 `list(df=…, summary=…)`（`summary` 含匹配统计 + 字段 NA 表 + 简化列清单），并在文档中说明；若为兼容现有 app 调用可先加参数 `return_summary=FALSE` 过渡。

---

### 2.5 `cor_cancer_genelist()`（R/cor_cancer_genelist.R）——单队列基因×基因相关

**功能解读**：取 `dataset1`+`id1` 与 `dataset2`+`id2` 两组表达；去数据集后缀后按 ID/type/dataset 内连接；过滤 sample_type；第 4 列与其余列 `psych::corr.test`；返回 `list(cor_result, cor_data)`。

**问题（含确定性 bug）**：
1. **NULL 守卫失效（实测确认）**：`data1$dataset<-` 先于 `if(is.null(data1))` 执行——对 NULL 做 `$<-` 会把 NULL 变成 list（R 语义），守卫永远不触发；且 `data2` 完全没有守卫。id1/id2 查不到数据时会一路变成晦涩错误而非"基因未找到"。
2. 合并前 `intersect(ID)`、合并本身均静默丢样本/丢基因；`id2` 中部分基因缺失无提示。
3. `sample_type` 不校验（拼错→0 行→corr.test NA 或报错）。
4. 结果列名拼写错误 `"Correlation efficience"`；r/p 矩阵用 `result[1]`/`result[4]` 位置索引，依赖 psych 内部结构。
5. 无有效样本量 n 输出；小样本/常数列产生 NA 无说明。

**建议**：把 is.null 守卫提到任何 `$<-` 之前并对两个数据集都做；合并后报告匹配样本数与丢失数；按 `id2` 逐基因报告缺失；`sample_type` 用 match.arg 并校验非空；corr.test 包 tryCatch，n<3 或方差为 0 时输出 NA 并注明原因；返回结果附 n 向量；改列名 `"Correlation"`。

---

### 2.6 `cor_pancancer_genelist()` / 2.7 `cor_pancancer_drug()` / 2.8 `cor_pancancer_TIL()`

**功能解读**（三者同构，仅 y 侧数据源不同）：将目标基因宽表 `df`（列 4 为目标基因表达）与 `geneset_data` / 药物矩阵 `drug_CPTAC` / 免疫浸润 `TIL_CPTAC` 按 ID（±type/dataset）合并 → 按 dataset 拆分 → 逐数据集 `corr.test(x=第4列, y=各特征列)` → 组装 `dataset×特征` 的 r/p 矩阵 → `na.omit` 后转置 → 返回 `list(r, p, sss)`（sss 为按数据集拆分后的原始数据）。

**问题（缺失数据集中爆发区）**：
1. **样本覆盖静默**：药物/TIL 数据只覆盖部分样本（1553/1636 例），内连接丢弃大量样本，无任何匹配率报告；`cor_pancancer_TIL` 还先强制 `type=="Tumor"`。
2. **`if(nrow(sss_can)<4) next` 静默跳过**：这些数据集在 r/p 中留 NA 行，随后 `na.omit()` 把**含任一 NA 的整行/整列删掉**（转置前后各一次），用户完全不知道哪些 dataset、哪些基因因样本不足/全 NA 被剔除；极端情况全部被删 → 返回空矩阵，无提示。
3. `corr.test` 每个数据集算两遍（r 一遍、p 一遍，浪费一倍时间）；x/y 含全 NA 或常数列时返回 NA 而无说明（实测确认 corr.test 在此类输入下不报错、静默给 NA）。
4. `drug`：`sig` 来自 `drug_info` 按 `Target.pathway` 过滤后拼出的 `Name_ID`——若某药物列在 `drug_CPTAC` 中不存在，`sss_can[, sig]` 会下标越界报错；pathway 拼错时 sig 可能为空。**无"该通路无药物/拼写建议"反馈**（可用 `unique(drug_info$Target.pathway)` 校验）。
5. `TIL`：`TIL_type` 不经 `TIL_map` 校验，无效算法→sig 为空/列缺失→错误或空结果，无提示（应提示 `unique(TIL_map$algorithm)` 候选）。
6. `cor_pancancer_genelist` 中 `sig <- colnames(df)[5:ncol]` 在**合并后**的对象上取列——若 `df` 本身含多个基因（函数文档并未禁止），第 5 列起会把 `df` 自身的第 2 个基因也当 y 计算，结果误导而无提示。
7. dataset 去后缀后作为行名/列名：若同一队列的 protein 与 mRNA 同时参与，行名碰撞（都变成 `LUAD_APOLLO`）。
8. 返回结构与 `cor_cancer_genelist`（`cor_result/cor_data`）不一致，下游/用户需分别适配。

**建议（此三函数应重构为共享核心）**：
- 抽公共函数如 `pancorr_core(df_target, feature_matrix, by, min_n, …)`：逐 dataset 记录**有效样本数 n** 与"因何未算"（样本<min_n / 目标基因全 NA / 特征全 NA / corr.test NA），返回 `r`、`p`、`n` 三个等维矩阵 + `skip_reason` 说明，message 汇总"跳过数据集 X（n=3<4）"。
- 彻底去掉对结果矩阵的 `na.omit()`，缺失以 NA 呈现，交给 viz 层用专门样式表达（见 2.11），绝不在分析层静默删行。
- 药物/TIL 合并前统计匹配样本；`Target.pathway`/`TIL_type` 先与 `drug_info`/`TIL_map` 校验，非法值给候选列表提示。
- `cor_pancancer_genelist` 校验 df 为单目标基因（或增加参数显式指定目标列），避免第 2 基因混入 y。
- 统一三函数返回契约 `list(r, p, n, data)`；与 `cor_cancer_genelist` 对齐命名。

---

### 2.9 `viz_TvsN()`（R/viz_TvsN.R）——Tumor vs Normal 差异可视化

**功能解读**：三种模式——`single`（单基因单数据集）、`multi_gene`（多基因 melt 后箱线）、`multi_set`（单基因多数据集箱线）；用 `ggpubr::compare_means`（t.test/wilcox）算 p 与显著性符号；可选在底部标 n、图上标 p/p.signif。

**问题**：
1. **默认参数直接报错（实测复现）**：`df_type` 默认 `c("single","multi_gene","multi_set")` 是三元素向量，而代码是 `if (df_type == "single")`——不传 df_type 调用必然报 `the condition has length > 1`。应 `match.arg(df_type, c("single","multi_gene","multi_set"))`。
2. **na.omit 静默删行**：三种模式都先 na.omit，删除量、删除发生在哪个基因/数据集均无报告。
3. **只有单组时直接崩/无意义**：实测 19 个数据集无 Normal 样本（dataset_info$Normal=NA）。对这类数据集做 Tumor vs Normal 时 filter 后只剩一组，`compare_means` 报错或返回空 pv，图上无任何解释。
4. 样本量过小（每组 <3）无提示。
5. `Show.P.value` 默认 TRUE 对 multi_gene × 多数据集会逐个检验，遇单组组合即失败，全图崩溃而非跳过该组合。
6. 文本标签的 y 用全局 min/max 估算，与分面/多组布局可能重叠（展示层问题）。

**建议**：`match.arg` + 前置校验（必需列 ID/type/dataset/表达列）；统计"输入行、NA 剔除行、Tumor 组 n、Normal 组 n"，逐 gene×dataset 报告；仅单组的组合跳过检验并在图上以文字标注 "No Normal samples in dataset X"（而不是报错或静默）；n<3 的组合标 "insufficient (n<3)"；返回 `invisible(list(plot=p, stats=pv, n=counts))` 便于 app 展示与下载。

---

### 2.10 `viz_DEGs_volcano()`（R/viz_DEGs_volcano.R）——DEG 火山图

**功能解读**：按 `p.cut`/`logFC.cut` 把基因标为 Up/Down/No（用 `adj.P.Val`），按 logFC 排序后可选标注首尾基因（`show.top`）或指定基因（`show.labels`）；`-log10(adj.P.Val)` 作 y 轴。

**问题**：
1. 无必需列校验（Symbol/logFC/P.Value/adj.P.Val）：若输入表缺列，`aes` 报错晦涩难懂。实测 API 的 ttest/limma 结果都含这些列，但函数应对更宽输入/错误输入给明确反馈。
2. `p.cut`/`logFC.cut` 非数值、负数、或 p.cut≥1 均未校验。
3. `show.top=TRUE` 且 `nrow(df)<10` 时 `df[(nrow-4):nrow, ]` 产生非法索引报错；应先检查行数并提示。
4. `show.top` 与 `show.labels` 同时为真时后者覆盖前者的 label（`if` 顺序），无提示。
5. P 值/adj.P.Val 极小（如 <1e-300）时 `-log10` 产生 Inf，geom 会警告并吞点。

**建议**：入口校验必需列并报告缺失列名；`show.top` 需 nrow≥10 否则 warning 并降级为不标注；把 `-log10` 结果截断（如 `pmin(..., 300)`）；语义冲突（top+labels 同开）给 message；允许自定义 cut 方向（如 |logFC| 双侧已有）与配色。

---

### 2.11 `viz_cor_heatmap()`（R/viz_cor_heatmap.R）——相关结果热图

**功能解读**：把 r 矩阵按行（特征）hclust 聚类；melt r、p 两矩阵并按键合并；按 p 值打 `*`/`**`/`***`；双色渐变 tile；行数>1 时用 `ggtree` + `aplot::insert_left` 在左侧拼聚类树。

**问题**：
1. **`r[is.na(r)] <- 0` 严重误导（缺失数据可视化核心问题）**："没算出来/样本不足"的格子被画成 r=0（无相关）的白色，与"真实测得无相关"无法区分。p 矩阵 NA 未处理，melt 后 `case_when` 得 NA 标签 → `geom_text` 空 + 警告。应保留 NA 格子并画成专门颜色（`scale_fill_gradient2(na.value="grey90")`），图例或标题注明"灰格=无数据"，并在 message 中报告缺失格数量。
2. **未声明依赖 ggtree/aplot**：`ggtree(...)`、`insert_left(...)` 未限定命名空间，DESCRIPTION Imports 与 NAMESPACE 都没有这两个包（也不在 README 安装清单里）——干净环境下运行到聚类分支即报 `could not find function`。这是"依赖缺失反馈"缺失的典型。
3. r/p 两个 melt 结果 merge 时若无完全匹配（如维度不一致）会错位/多行，无校验。
4. 聚类对全 0/全 NA 行：NA 已被填 0，若整行同值则 dist=0，hclust 可运行但结果无意义，无提示。
5. 显著性阈值硬编码（0.05/0.01/0.001）。

**建议**：NA 格可视化 + 计数提示（见上）；`aplot`/`ggtree` 加入 DESCRIPTION（Suggests 或 Imports），函数内 `requireNamespace` 检查并在缺包时给"请安装 ggtree/aplot（BiocManager）"的明确提示、跳过聚类树降级为纯热图；r/p 维度一致性校验；显著性 cut 参数化（`sig_cuts`）；可选在格内显示 r 值。

---

### 2.12 `viz_corplot()`（R/viz_corplot.R）——双基因散点 + 相关标注

**功能解读**：取 data 中 a、b 两列画散点 + lm 平滑 + rug；标题由内嵌 `corr_eqn` 用 `cor.test` 给出 n、r、p。

**问题**：
1. **NA 未处理（实测）**：x/y 含 NA 时 `cor.test` 不报错但 r/p=NA，标题显示 `n = 全部行数`（含 NA 行）+ `r = NA`，误导用户；`geom_smooth(na.rm=T)` 与点图不一致。
2. `a`/`b` 列不存在时报原生 `undefined columns selected`；未校验。
3. `method` 不校验；`exact=FALSE` 只对 spearman 有含义。
4. n<3 或常数列时 r/p=NA 无说明。

**建议**：先 `complete.cases` 计数并提示"剔除 NA 行 k（共 N 行中的 n 有效）"，用有效样本重算 n/r/p；列存在性校验给明确 message；`match.arg(method, c("pearson","spearman"))`；n<3 时明确标注不可算。

---

### 2.13 `viz_phoso_sites()`（R/viz_phoso_sites.R）——蛋白磷酸化位点示意图

**功能解读**：经 `uniport_map` 找 Reviewed 的 UniProt Entry → `drawProteins::get_features`（联网 UniProt）拉结构域 → `feature_to_dataframe` → draw_canvas/chains/domains/regions/motif；CPTAC 模式从 `idmap_protein` 中该基因所有 `row_names`（形如 `NP_000025.1:s218y223t227`）解析位点打红点；UniProt 模式用 `phospho_site_info(rel_data)`。

**问题**：
1. **位点解析 bug（确定）**：内层 `for(i in 1:length(aa)){ bb <- substr(...) }` 每轮**覆盖** bb，只保留最后一个片段。实测 id `NP_000025.1:s218y223t227` 只会输出 `t227`，`s218`、`y223` 全部丢失；外层/内层循环复用变量 `i`（虽能运行但极易出错）。
2. `idmap_protein` 中该 Symbol 的**蛋白行（无 `:site` 后缀）也会混入**被当位点解析出垃圾（如 TP53 有 5 行：1 蛋白 + 4 位点行）。
3. 基因不在 `uniport_map`（或非 reviewed）→ `uni_id=character(0)` → `get_features` 行为未定义且无提示；UniProt 为联网调用，断网/限流无 tryCatch 反馈。
4. 多个 reviewed Entry 只 `unique()` 后仍当单 Entry 处理；`From[1]` 前缀假设对全部行成立。
5. 无返回（只有 ggplot 对象），位点表无法导出核对；拼写 "phoso" 与多处死代码。

**建议**：位点解析改为一次正则 `regmatches(id, gregexpr("[sty]\\d+", id))` 收集全部片段（修复丢失位点的 bug）；仅取含 `:` 且带位点的行参与画点；Entry 空/网络失败时给出明确 message（"未在 UniProt 找到 reviewed entry for gene X / 网络错误，已跳过特征绘制"）；返回 `invisible(list(plot=p, sites=p_data))` 供核对与下载。

---

### 2.14 `PCAS_app()` 与 Shiny 前端（inst/shinyapp/）

**功能解读**：`PCAS_app()` = `shinyAppFile(system.file("shinyapp","app.R",...))`。app.R 加载十余个包，10 个模块对应各函数；反馈模块收集 issue 并**用硬编码的 QQ 邮箱 SMTP 账号密码发信**。

**问题**：
1. 依赖安装策略 `if(!require()) install.packages(...,ask=FALSE)` 需要联网且不可控（含 Bioc 包），离线/企业环境直接失败，无诊断信息。
2. 前端对"部分缺失"无感知：各模块只对 `is.null(...)`（整体无数据）弹窗；部分基因缺失、部分样本被内连接丢弃、部分数据集被 na.omit 剔除时**用户看到的就是少了几行/几格，且无解释**（package 层的 message 只进 console）。
3. **安全**：mod_feedback.R 中 SMTP 密码明文硬编码并提交到公开仓库（应立即轮换该密码并改为环境变量/服务端代理），属高危项。
4. app 引用的 `colourpicker`、`shinycssloaders` 等未在 DESCRIPTION Suggests 中登记，R CMD check 无法发现断链。

**建议**：依赖改为安装脚本/`renv` + 明确错误提示；把 package 层新增的"可用性摘要"（availability attr / summary）渲染成可下载表格或状态提示；删除/轮换 SMTP 凭据。

---

### 2.15 包元数据与文档（DESCRIPTION/NAMESPACE/man/README）

1. `DESCRIPTION`：`Date: 2024-061-19` 非法；`License: None`（应改为 GPL/MIT 等真实许可，README 中论文在 IJMS 发表，建议与作者确认）；**Imports 缺 `stringr`、`digest`、`plyr`**（代码大量 `stringr::`/`digest::`/`plyr::`），`viz_cor_heatmap` 还使用未声明的 `ggtree`/`aplot`；`import(...)` 整包导入过宽（NAMESPACE 24 行），建议收紧为 `@importFrom`。
2. `R/data.R`：10 个数据对象只有 9 条记录且**缺 `TIL_map`**；每条缺少 `@name/@format/@source/@docType data` 完整块，`man/` 下没有生成任何数据集 Rd（`?dataset_info` 查不到）。
3. `man/hello.Rd` 为包模板残留（对应函数不存在）。
4. roxygen/README 文档错位：`viz_DEGs_volcano` 的 Rd 与 README 参数表写的是 `cohort`/`data_input`（真实形参为 df/p.cut/logFC.cut/show.top/show.labels）；`get_data` 文档 `gene` vs 形参 `genes`；README 多个 markdown 编号错误。
5. 输出列名 `"Correlation efficience"` 拼写错误；函数命名 `phoso` 应为 `phospho`（改名为破坏性变更，可在新版本以别名兼容）。

---

## 3. 缺失数据反馈：专项设计建议（横向框架）

把上面各函数的问题归纳为一条可落地的五层框架：

**① 识别层（在哪缺失）**：所有取数/合并函数记录三类缺失——(a) 查询缺失：请求了但基因/位点在数据集/映射表里不存在；(b) 合并缺失：因内连接（临床/药物/TIL）被丢弃的样本；(c) 值缺失：列内 NA、全 NA 列、na.omit 删除的行。`get_expr_data` 产出 `availability`（dataset×gene）矩阵；各 merge 产出 matched/unmatched 计数。

**② 校验层（尽早失败并给建议）**：统一输入校验：dataset 是否在 `dataset_info$Abbre`（给近似候选）、gene 是否在 idmap、`sample_type`/`method`/`cor_method`/`Target.pathway`/`TIL_type` 是否合法（给 `unique()` 候选）、data.frame 必需列是否存在。全部用"错误/警告 + 可执行建议"文案。

**③ 记录层（内部量化）**：`n_obs`、`n_na_removed`、`n_unmatched`、`n_skipped_datasets` 等随结果携带（attr 或 summary list），Shiny 端可渲染，杜绝"console message 在 GUI 里隐形"。

**④ 呈现层（图内表达缺失）**：热图 NA 用 `na.value` 专属灰格 + 图例说明，禁止 `NA→0`；散点/箱线注明剔除数与 n；TvsN 对无 Normal 数据集画"不可比"标注；相关结果输出 n 矩阵并允许在图上以字号/标注表达。

**⑤ 文档层**：所有 Rd 明确：NA 语义、返回 NULL 的时机、哪些步骤静默剔除（改造后应声明为"不再静默"）、缓存行为、联网依赖（API/UniProt）；README 参数表与真实签名逐一对齐。

---

## 4. 改造优先级

### P0 —— 正确性 / 安全隐患（先修）
1. `viz_TvsN`：`df_type` 改 `match.arg`，默认值可用（否则函数无法缺省调用）。
2. `cor_cancer_genelist`：NULL 守卫移到 `$<-` 之前并为 data2 补守卫。
3. `viz_phoso_sites`：修复位点解析（只保留最后一个片段 bug）。
4. `viz_cor_heatmap`：删除 `r[is.na(r)] <- 0`，改 NA 专属样式。
5. `DESCRIPTION`：补 Imports（stringr/digest/plyr，及 ggtree/aplot 或降级逻辑）；修 Date/License；`man/hello.Rd` 删除。
6. `mod_feedback.R`：移除硬编码 SMTP 密码（安全）。
7. `cor_pancancer_*`：去除对结果矩阵的 `na.omit()`（改 NA 透传 + 呈现层处理）。

### P1 —— 缺失数据反馈（本任务重点）
8. `get_expr_data`：预检 idmap + 逐基因/逐数据集缺失汇总 + `availability` 属性 + 单/多基因统一路径 + 缓存位置/开关重构。
9. 三个 `cor_pancancer_*`：样本匹配统计、`n` 矩阵返回、跳过原因 message、pathway/TIL_type 参数校验（含候选提示）。
10. `merge_clinic_data`：clinic 失败/`data_input` 空校验、匹配计数、临床字段 NA 汇总、summary 返回。
11. `get_data`：HTTP 状态/超时/tryCatch + URL 编码 + action 校验，区分"无数据"与"服务错误"。
12. `viz_TvsN`/`viz_corplot`/`viz_DEGs_volcano`：NA 行数提示、单组数据集的"不可比"标注、必需列校验。
13. `get_DEGs_result`：method/dataset 校验、merge 丢失行报告。

### P2 —— 工程质量 / 一致性
14. 三/四个相关函数重构为共享核心，统一返回契约（r/p/n/data）与命名（修 "Correlation efficience"，`phoso`→`phospho` 别名兼容）。
15. NAMESPACE 收紧为 `@importFrom`；Roxygen 补齐 @return/@param/@examples 与实际签名一致。
16. 数据文档：补 `TIL_map`、生成数据集 Rd；README 逐函数勘误。
17. app：模块级错误提示统一从函数返回的 summary 渲染；依赖安装方案替换。

---

## 5. 实测证据附录

| 实验 | 结果 |
|---|---|
| `get_expr_data("LUAD_CPTAC_protein", c("TP53","ZZZNOTAGENE"))` | 返回 214×4 仅含 TP53，**无任何缺失提示** |
| `get_expr_data(..., "ZZZNOTAGENE")`（整组无效） | 消息 `… retrive no results`（拼写错误）+ `Retrive no data.`，返回 NULL |
| `if (df_type == "single")` 且 df_type=默认三元素 | `Error: the condition has length > 1` |
| `z <- NULL; z$dataset <- "x"` | z 变成 list（NULL 守卫因此失效） |
| `cor.test` 含 NA | 不报错，r/p=NA（标题误导） |
| `corr.test` n=3 / 常数列 / 全 NA 列 | 不报错，r/p 静默为 NA |
| API expression | 每基因 1 行宽表；样本列 `…_Tumor/_Normal` 及三个 `IR_Other` 池 |
| API DEGs（ttest 与 limma） | 均含 adj.P.Val/change 等列（volcano 输入基本齐备） |
| API clinic（LUAD_APOLLO） | `Days_To_Recurrence` 66/101 NA、`Pack_Years_Smoked` 34/101 NA |
| idmap 磷酸化位点 id | 存在 `NP_000025.1:s218y223t227` 组合格式 → 旧解析只留 `t227` |
| dataset_info | 56 个数据集中 19 个 Normal=NA（无正常样本，viz_TvsN 会单组化） |

---

## 附录：实施与验证状态（2026-01 更新）

本报告的 P0/P1/P2 建议已全部实施（新增/重写文件：`R/get_data.R`、`R/get_expr_data.R`、
`R/get_DEGs_result.R`、`R/merge_clinic_data.R`、`R/cor_cancer_genelist.R`、
`R/cor_pancancer_{genelist,drug,TIL}.R`、`R/pancorr_engine.R`、`R/viz_TvsN.R`、
`R/viz_DEGs_volcano.R`、`R/viz_cor_heatmap.R`、`R/viz_corplot.R`、
`R/viz_phoso_sites.R`、`R/internal.R`、`R/globals.R`、`R/data.R`、`DESCRIPTION`、
`NAMESPACE`、`man/*.Rd`、`README.md`、`tests/smoke_tests.R`；App 适配：
`mod_feedback.R`、`modules_Cancer_expression.R`、`modules_Cancer_correlation.R`、
`modules_Cancer_multiple.R`、`modules-cptac-search.R`、`modules-cptac-site.R`、
`modules-pancan-corr.R`）。

### P0 落实
1. `viz_TvsN`：`df_type` 使用 `match.arg`，默认 `"single"`，缺省调用不再报错（有测试）。
2. `cor_cancer_genelist`：NULL 守卫移到任何 `$<-` 之前，dataset1/dataset2 均守卫。
3. `viz_phoso_sites`：位点解析改为正则一次提取全部 `[sty]\d+`，组合位点不再丢；
   只取含 `:` 的位点行；空 Entry/网络失败给出明确 message 并返回 NULL。
4. `viz_cor_heatmap`：删除 `r[is.na(r)] <- 0`，NA 格用 `na.value="grey85"` 独立呈现并计数提示；
   缺 ggtree/aplot 时降级为纯热图并提示；r/p 维度/行列名一致性校验。
5. `DESCRIPTION`：Date 修正、License: MIT + file LICENSE、Imports 补齐
   digest/dplyr/httr/scales/stringr/ggpubr/ggrepel 等、Suggests 加 aplot/ggtree；
   `man/hello.Rd` 删除。
6. `mod_feedback.R`：硬编码 SMTP 密码移除，改为 `PCAS_SMTP_*` 环境变量，缺省/失败均以 UI 提示。
7. 三个 `cor_pancancer_*`：删除结果矩阵的 `na.omit()`，NA 透传由引擎/热图呈现。

### P1 落实（缺失数据反馈）
8. `get_expr_data`：数据集名预检；mRNA 先经 idmap 预检（未映射基因明确提示并跳过请求）；
   逐数据集×基因 availability 汇总并挂 `attr(,"availability")`；返回列恒为用户请求顺序、
   缺失基因保留全 NA 列；多转录本按样本均值聚合；缓存目录可配（`options(PCAS.cache.dir)`）、
   键含月份防陈旧、写失败降级；单/多基因统一路径；消息文案统一。
9. 相关三函数：样本匹配/未匹配计数与 message、`n`（逐格有效样本数）矩阵、跳过原因 summary、
   pathway/TIL_type 参数校验并给候选提示；`df` 多表达列时提示仅用第一列。
10. `merge_clinic_data`：输入/clinic 校验；匹配/未匹配计数；临床字段缺失率表；`*_simplify`
    生成清单；默认返回 `list(df, summary)`（可用 `return_summary=FALSE` 取旧式 data.frame）。
11. `get_data`：httpl 超时+重试+HTTP 状态检查、URL 逐项编码、action 校验、
    "无数据(NULL+message)"与"服务错误(NULL+warning)"区分。
12. `viz_TvsN`/`viz_corplot`/`viz_DEGs_volcano`：NA 剔除行数提示、单组数据集的"跳过检验"提示、
    必需列校验（volcano 缺 adj.P.Val 回退 P.Value）；corplot 用完全观测重算 n/r/p。
13. `get_DEGs_result`：method/dataset 校验、mRNA merge 丢失行报告、cache 重构。

### P2 落实
14. 相关三函数共用 `R/pancorr_engine.R` 核心，统一 `list(r, p, n, sss, summary)`；
    `cor_cancer_genelist` 返回 `cor_result/cor_data/n/summary` 并修正 "Correlation" 列名。
15. NAMESPACE 收紧为仅导出+shiny 导入；全部外部调用 `pkg::fun`；新增 `R/globals.R` 声明
    ggplot aes 列名；codetools/R CMD check 无 "no visible binding" 类问题。
16. `R/data.R` 补全 10 个数据对象文档（含此前缺失的 `TIL_map`），roxygen2 重生成全部 man。
17. App：所有消费点已适配（`merge_clinic_data` 列表返回、NULL 守卫与 req、phoso 站点表改读
    `attr(p,"sites")`、pancan-corr 热图 req、下载前判空）。

### 验证
- 全部 17+13 个 R/App 文件 `parse()` 通过。
- `man/*.Rd` 全部通过 `tools::checkRd`。
- 集成冒烟测试 `tests/smoke_tests.R`（真实 PCAS API + 本地构造数据）：
  **48 项断言全部通过（48 passed / 0 failed）**，覆盖上述每项行为与缺失数据反馈路径。
- 包在本机以 `R CMD INSTALL` 安装成功；`R CMD check --no-manual` 结果见仓库验证日志。
- 已知环境限制：UniProt 网络调用受网络环境影响，`viz_phoso_sites` 测试接受
  "返回 ggplot"或"网络失败时优雅 NULL+提示"两种结果。
