# PCAS 更新说明 / Changelog

## 0.2.0 (2026-09-05)

本版本是一次系统性的代码加固发布：在保持全部 14 个导出函数接口可用（并已同步适配
Shiny App）的前提下，修复了多个正确性 bug，全面补齐"缺失数据"的识别、记录与用户
反馈，并完善了工程规范与文档。

This is a hardening release: it keeps every exported function usable (the bundled
Shiny app was adapted in sync) while fixing several correctness bugs, adding
comprehensive missing-data detection/feedback, and improving package hygiene
and documentation.

### 中文更新摘要 / Summary (中文)

1. **缺失数据反馈（核心）**
   - `get_expr_data()`：请求前先校验数据集名；mRNA 符号先经 `idmap_RNA` 预检，
     未映射的符号会提示且不再发出空查询；输出恒为用户请求的基因顺序，某数据集
     测不到的基因保留为全 NA 列（不再被静默删除）；结果附带
     `attr(result, "availability")`（dataset × gene 的 present/n_valid 明细）；
     同基因多转录本按样本取均值并提示；所有跳过/缺失均给出说明消息。
   - 三个泛癌相关函数 `cor_pancancer_genelist()/cor_pancancer_drug()/
     cor_pancancer_TIL()`：改为共享计算引擎；结果新增逐格有效样本数矩阵 `n` 与
     `summary`（跳过原因、缺失特征、样本匹配统计）；彻底移除对相关矩阵的
     `na.omit()` 静默整行删除；`Target.pathway`/`TIL_type` 等参数先校验并给出
     候选值提示。
   - `merge_clinic_data()`：报告表达样本与临床记录匹配/未匹配数量、临床字段缺失率、
     生成的分期简化列清单；默认返回 `list(df, summary)`（旧用法可用
     `return_summary = FALSE` 恢复）。
   - `get_DEGs_result()`：`dataset`/`method` 校验；mRNA 行无法映射到基因符号时报告
     丢弃行数。
   - 可视化函数均报告被剔除的 NA 行数；热图将"未计算"格子以灰色独立呈现
     （不再把 NA 当作相关系数 0）。

2. **正确性修复**
   - `viz_TvsN()` 的 `df_type` 改为 `match.arg` 校验（默认可直接调用，不再报
     "condition has length > 1"）；仅含单一样本类型的数据集跳过检验并提示。
   - `cor_cancer_genelist()`：NULL 守卫移到任何转换之前，两个数据集都检查；
     修正输出列名 `Correlation`，新增每特征样本数 `n`。
   - `viz_phoso_sites()`：修复组合位点（如 `NP_…:s218y223t227`）只画出最后一位点
     的解析 bug；空基因/无 reviewed Entry/UniProt 网络失败时给出明确提示并优雅返回
     NULL；站点表以 `attr(p, "sites")` 供导出。
   - `viz_cor_heatmap()`：维度/行列名校验；缺少 ggtree/aplot 时自动降级为纯热图。
   - `get_data()`：HTTP 超时/重试/状态检查、查询参数 URL 编码、action 校验；
     "无数据"（message）与"服务错误"（warning）区分。
   - 缓存：目录改为可配置（`options(PCAS.cache.dir)`），键含月份防止陈旧数据，
     只读文件系统上写入失败自动降级。

3. **工程与文档**
   - `DESCRIPTION`：License 修正为 `MIT + file LICENSE`、补全 Imports/Suggests、
     数据启用 xz 压缩；删除模板文件 `man/hello.Rd`。
   - 补齐 10 个内置数据对象的 Rd 文档（含此前缺失的 `TIL_map`）；README 参数表勘误。
   - Shiny App：适配 `merge_clinic_data()` 新返回结构并增加空结果提示；反馈模块不再
     硬编码 SMTP 密码，改由 `PCAS_SMTP_*` 环境变量提供；多处 `req()`/判空防崩溃。
   - 新增自动化冒烟测试 `tests/smoke_tests.R`（48 项断言）。
   - 安全：删除 `mod_feedback.R` 中明文邮箱授权码。

### Bug fixes (English)
- `viz_TvsN()`: `df_type` is now validated with `match.arg` and works with its
  default value; single-group datasets skip the significance test with a message.
- `cor_cancer_genelist()`: NULL guards run before any transformation on both
  datasets; corrected result column name (`Correlation`); added per-feature `n`.
- `viz_phoso_sites()`: combined phosphosite ids are now parsed into *all* sites
  (previously only the last one was drawn); graceful NULL + message on unknown
  gene / missing reviewed UniProt entry / network failure.
- `viz_cor_heatmap()`: no longer turns missing values into a coefficient of 0;
  dimension/name checks; falls back to a plain heat map without ggtree/aplot.
- `get_data()`: request timeouts/retries, HTTP status checks and per-argument
  URL encoding; distinguishes "no data" from "server/network error".

### New features & missing-data feedback (English)
- `get_expr_data()`: dataset-name pre-validation; mRNA symbols pre-checked against
  `idmap_RNA`; requested-order columns kept (all-NA when a dataset does not
  measure an identifier); `attr(result, "availability")` matrix; multi-probe
  symbols averaged per sample with a message; configurable cache
  (`options(PCAS.cache.dir)`) with monthly invalidation.
- Pan-cancer correlation functions: shared engine; new `n` matrix of
  pairwise-complete sample sizes and a `summary` element; `na.omit()` removed
  from result matrices; `Target.pathway` / `TIL_type` validated with suggestions.
- `merge_clinic_data()`: matching/missingness report; returns
  `list(df, summary)` by default (`return_summary = FALSE` keeps the old shape).
- `get_DEGs_result()`: dataset/method validation; reports mRNA rows dropped by
  the symbol mapping.
- Visualizations report rows removed for missing values; heat maps render
  "not computed" cells in grey instead of as r = 0.

### Documentation & internals (English)
- DESCRIPTION: `MIT + file LICENSE`, completed Imports/Suggests,
  `LazyDataCompression: xz`; removed the template `man/hello.Rd`.
- Full Rd documentation for all 10 bundled datasets (incl. the previously
  undocumented `TIL_map`); README argument tables corrected.
- Shiny app adapted to the new returns with empty-result guards; feedback module
  reads SMTP credentials from `PCAS_SMTP_*` environment variables (no more
  hard-coded password); several `req()` guards added.
- New smoke/integration test suite `tests/smoke_tests.R` (48 assertions).
- Validation: `R CMD check --no-manual` reports **Status: OK** (0 errors /
  0 warnings / 0 notes).
