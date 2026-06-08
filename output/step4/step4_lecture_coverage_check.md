# 第四讲覆盖核对与执行说明

## PDF核对说明

`讲义/经济计量方法导论 第四讲.pdf` 可渲染封面，标题为“第四讲 简单截面数据的经济计量方法（三）”。该 PDF 的文本层无法通过本地 `textutil`、`strings` 有效提取；本地二进制页标记显示约 30 页。因此，本文件基于封面主题、`what_we_chat_before/第三步_第四讲.md` 中列出的第四讲知识点，以及第四讲通常围绕“异方差诊断、稳健推断与效率改进”的结构进行核对。

## 知识点覆盖

1. 高斯—马尔可夫假定4：同方差性。已在 `step4_model_decision_notes.txt` 中说明为什么行业截面可能不满足同方差。
2. 异方差导致的问题。已在说明中明确：OLS点估计可保留，但普通标准误、t检验和F检验可能不可靠。
3. 残差图诊断。已输出 residuals vs fitted、absolute residuals vs fitted、residuals vs RateVulnerability、四分位残差方差、scale-location 图。
4. BP检验。已输出 `step4_heteroskedasticity_tests.csv`。
5. White检验。已输出基于拟合值及平方的改进White检验，并补充完整White检验。
6. 异方差稳健标准误。已输出普通SE、HC0、HC1、HC2、HC3 对比表。
7. 稳健联合检验。已输出基于 HC3 的控制变量联合 Wald/F 检验。
8. WLS。已用 `1 / HistVol^2` 构造经济含义明确的权重，并作为稳健性检验。
9. FGLS。已用 `log(u_hat^2)` 辅助方程估计方差函数，并用 `1 / h_hat` 加权。
10. FGLS/WLS 后再诊断。已对 WLS 和 FGLS 输出改进 White 检验。
11. 最终处理原则。主结论以 OLS + HC3 为主，WLS/FGLS 作为稳健性与效率改进，不替代主模型。

## 当前结果判断

- BP 检验 p = 0.3465，改进 White 检验 p = 0.1956，未拒绝同方差原假设。
- 考虑到样本只有 49 个行业，仍采用 HC3 作为保守稳健推断。
- M2 主模型中 `RateVulnerability_z` 在 HC3 下系数为 -1.0979，p = 0.0472，仍在 5% 水平显著。
- WLS 中核心系数仍为负但不显著，说明按历史波动率降权后统计证据变弱。
- FGLS-HC3 中核心系数为 -1.2446，p = 0.0253，方向和显著性支持主结论。
- 因此第四讲结论应写成：异方差检验未显示强烈证据；在更保守的 HC3 稳健标准误下，利率脆弱度的负向效应仍成立；WLS/FGLS 作为补充稳健性而非替代基准模型。

## 修改意见补充记录

1. 已修正完整 White 检验交互项生成逻辑，8 个解释变量下自由度从 45 修正为 44；完整 White 检验 p 值更新为 0.3101。完整 White 检验变量较多，仍不作为主文重点。
2. 已将辅助回归 LM 检验自由度从 `length(coef)-1` 改为 `aux_model$rank - 1`，避免共线或 alias 变量影响自由度。
3. 已新增影响点诊断 `step4_influence_diagnostics.csv`，包含 Cook's distance、hat value 和 `RateVulnerability_z` 的 DFBETAS。当前 Cook's distance 最大的是 Gold。
4. 已新增 FGLS-White：用 `fitted + fitted^2` 估计方差函数，对应改进 White 检验逻辑。HC3 下核心系数为 -1.1307，p = 0.0404。
5. 已新增 groupwise WLS：按 `RateVulnerability_z` 四分位组内残差方差加权，对应残差四分位图。HC3 下核心系数为 -1.0483，p = 0.0311。
6. 各补充加权模型中，除 `1 / HistVol^2` 的简单 WLS 外，核心系数方向均保持为负且在 HC3 下显著。最终口径仍以 OLS-HC3 为主，FGLS-White 和 groupwise WLS 作为补充稳健性。

## Gold 影响点补充处理

1. Gold 的 Cook's distance 最大，为 0.6207，说明它对 M2 回归结果有较强影响。但 Gold 在货币政策事件中具有明确经济含义：黄金行业同时受利率、通胀预期、避险需求和美元走势影响，因此不应简单视为错误观测或直接删除。
2. 已新增 `step4_influence_sensitivity.csv`：
   - 全样本 OLS-HC3：系数 -1.0979，p = 0.0472；
   - 剔除 Gold：系数 -0.7381，p = 0.1127；
   - 剔除 Cook's distance 前三行业 Gold、FabPr、Autos：系数 -0.5635，p = 0.2286。
3. 已新增 `step4_leave_one_out_M2_HC3.csv`：49 次逐一剔除后核心系数 49/49 次为负，30/49 次在 5% 水平显著，47/49 次在 10% 水平显著。
4. 论文口径应写为：影响点诊断显示 Gold 对估计精度有明显影响。剔除高影响行业后，核心系数方向仍为负，但显著性下降，说明负向关系的方向较稳定，而统计显著性部分依赖包含 Gold 等经济上特殊行业的全样本信息。因此 Gold 不应被机械删除，影响点检验作为敏感性分析报告。
