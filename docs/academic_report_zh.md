# 学术中文版报告 - MIMIC-IV ICU 早期升压药治疗策略 RWE 项目

## 题目
MIMIC-IV 中成人 ICU 脓毒症编码住院的早期升压药启动与院内死亡风险：基于 target-trial-style 设计的真实世界证据分析

## 摘要
本项目使用 MIMIC-IV ICU 数据，评估成人首次 ICU 入室且有脓毒症相关编码的住院中，早期启动升压药与院内死亡之间的关联。研究将 ICU 入室时间定义为 time zero。治疗策略组定义为 ICU 入室后 0-6 小时内首次启动升压药；对照组定义为 ICU 入室后 0-6 小时内未启动升压药。需要特别说明的是，对照组不是“从未使用升压药”的患者，而是在 6 小时 landmark 之前没有启动升压药的患者，之后可能存在延迟使用。最终分析队列包括 9,502 个 ICU stay，其中早期升压药组 3,772 例，对照组 5,730 例。粗院内死亡率分别为 32.9% 和 25.2%。经过 propensity-score weighting 后，13 个数值型协变量的加权后绝对 SMD 均小于 0.1，最大绝对 SMD 从 0.288 降至 0.046。AIPW 估计的早期升压药启动相对于 6 小时内未启动的风险差为 1.17 个百分点，95% CI 为 -0.52 至 2.87。该结果不能解释为升压药有害，而应理解为早期升压药使用本身高度反映急性病情严重程度；粗死亡率差异不能直接作为因果结论。

## 研究问题
在成人首次 ICU 入室且有脓毒症编码的住院中，ICU 入室后 6 小时内启动升压药与院内死亡风险之间的调整后关联是什么？

## 比较组定义

| 组成部分 | 定义 |
|---|---|
| Time zero | ICU 入室时间 |
| 治疗组 | ICU 入室后 0-6 小时内首次启动升压药 |
| 对照组 | ICU 入室后 0-6 小时内没有启动升压药 |
| 结局 | 院内死亡 |
| 关键说明 | 对照组不是 never treated，而是 6 小时 landmark 前未启动治疗 |

## 数据来源与队列
项目使用需要 credentialed access 的 MIMIC-IV ICU 数据。仓库不包含原始 MIMIC-IV 数据，也不包含患者级别衍生数据。代码可以基于本地 MIMIC-IV 文件夹或 BigQuery 导出的 cohort 运行。队列构建逻辑限制为成人、首次 ICU stay、脓毒症编码、并具有足够基线信息的住院记录。

## 方法
分析采用 target-trial-style 的观察性研究流程：

1. 定义 eligibility、time zero、treatment strategy、comparator strategy 和 outcome。
2. 从 ICU inputevents 中识别首次升压药启动时间。
3. 提取早期 ICU 窗口前可获得的基线特征。
4. 估计早期升压药启动的 propensity score。
5. 使用 SMD 评估 IPTW 前后的协变量平衡。
6. 拟合 propensity-score-weighted logistic regression，并用 AIPW 估计风险差。
7. 将粗死亡率与调整后结果分开解释。

最终版本将早期的 extraction + analysis 脚本和 analysis-only 脚本合并为一个主入口。后续 analysis-only 版本中的关键修复被保留下来：在 categorical dummy encoding 后，显式将设计矩阵转换成 numeric/float 类型，避免 `statsmodels` 因 bool/object 类型列报错。

## 结果

| 指标 | 结果 |
|---|---:|
| 分析队列 | 9,502 ICU stays |
| 早期升压药组 | 3,772 |
| 对照组 | 5,730 |
| 早期升压药比例 | 39.7% |
| 总体院内死亡率 | 28.3% |
| 早期升压药组粗死亡率 | 32.9% |
| 对照组粗死亡率 | 25.2% |
| IPTW 前最大绝对 SMD | 0.288 |
| IPTW 后最大绝对 SMD | 0.046 |
| IPTW 后达标数值协变量 | 13/13 个 |SMD| < 0.1 |
| AIPW 风险差 | 1.17 个百分点 |
| 95% CI | -0.52 至 2.87 个百分点 |

## 解释
早期升压药患者通常病情更重，因此粗死亡率更高不能直接解释为“升压药导致死亡率升高”。调整后风险差较小且统计不确定，说明粗差异很大程度上反映了基线病情严重程度和临床选择。该项目重点展示 RWE 中如何定义比较组、构建 time zero、进行 PS/IPTW 平衡诊断，并区分 crude association 与 adjusted estimand。

## 局限性
- 脓毒症和升压药暴露定义是 portfolio 项目的简化版本。
- 治疗启动由临床病情和医生判断决定，残余混杂很可能存在。
- landmark 和 immortal-time 相关设计选择需要进一步敏感性分析。
- 对照组不是 never-treated group。
- 结果不能作为升压药疗效或安全性的正式临床证据。

## 可复现性
运行方法见 `RUN_ORDER.md` 和 `notebooks/01_mimic_iv_early_vasopressor_pipeline.Rmd`。MIMIC-IV 原始数据因数据使用协议不包含在仓库中；仓库仅保留代码、说明文档和汇总级输出。
