# 第三讲覆盖核对

基于 `what_we_chat_before/第二部_第三讲.md` 与第三讲“简单截面数据的经济计量方法（二）”主题核对，本阶段已经覆盖第三讲核心知识点：

1. 多元线性回归的必要性：从一元模型升级为控制可观测遗漏变量的多元截面模型。
2. 偏效应解释：`RateVulnerability_z` 的系数解释为控制其他行业特征后的边际/偏相关。
3. 高斯—马尔可夫假定扩展：强调 `E(u|x1,...,xk)=0` 与可观测遗漏变量控制。
4. FWL 思想：输出 `step3_fwl_check.csv`，验证残差化后系数与完整模型一致。
5. 多重共线性：输出 `step3_vif.csv`，报告 tolerance 和 VIF。
6. 变量重要性：输出 `step3_standardized_coefficients.csv` 和标准化系数图。
7. t 检验和 F 检验：输出常规表、HC1 稳健表和嵌套模型 F 检验。
8. 虚拟变量与交互项：构造 TechGrowth、Finance、Defensive，并估计行业组斜率差异。
9. 模型设定诊断：用 RESET 型检验与利率脆弱度平方项检验函数形式误设。
10. 过度设定提醒：主模型保留 M0/M1/M2，交互模型作为扩展，不作为主模型。

执行层面说明：第三讲脚本重新构造了“条件利率脆弱度”，即在时间序列第一阶段同时控制 Fama-French 五因子和 `d_DGS2` 后得到的利率敏感度。因此第三讲的 M0 与第一步只用 `d_DGS2` 得到的简单利率脆弱度模型不完全相同。

## 修改意见补充记录

1. `M3_interactions_exploratory` 不作为主结论，因为它不是 M2 的完整扩展，且行业组样本量较小。
2. 已新增 `M3_full_interactions`，在 M2 全部控制变量基础上加入 TechGrowth、Finance、Defensive 及其与 `RateVulnerability_z` 的交互项。
3. 已新增 M2 leave-one-out 稳健性检验。49 次逐一剔除行业后，`RateVulnerability_z` 系数全部为负，48 次在 5% 水平显著，49 次在 10% 水平显著；最弱情形为剔除 Gold。
4. 已新增 `step3_coef_path.csv`，记录 M0、M1、M2 中核心系数随控制变量加入的变化，用于解释遗漏变量和偏效应。
5. 已新增 `step3_robust_wald_interaction.csv`。完整 M3 的 HC1 Wald 检验显示交互项联合显著，但常规 F 检验仅在 10% 附近，因此异质性结果应作为探索性补充，不作为主结论。
