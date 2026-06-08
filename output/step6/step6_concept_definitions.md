# Step 6：第六讲概念定义与本文变量对应关系

## 1. 为什么需要单独定义内生变量和外生变量？

第六讲的联立方程模型不是普通多元回归。它关注的是系统内部变量相互影响、共同决定的问题。因此，在估计之前必须说明哪些变量是系统内生决定的，哪些变量是系统外部给定的，哪些变量虽然出现在右侧但可能与结构误差相关。

## 2. 本文主收益结构方程

`CAR_3day = alpha0 + alpha1 RateVulnerability_z + controls + u1`

- `CAR_3day`：主结构方程的被解释变量，也是探索性 SEM 中的内生变量。
- `RateVulnerability_z`：主收益方程的右侧内生变量，因为它可能与遗漏行业特征相关。
- included exogenous controls：MktBeta_z, SMBBeta_z, HMLBeta_z, RMWBeta_z, CMABeta_z, HistVol_z, PreMomentum_z。
- excluded exogenous instruments：Z_precovid_2Y_z, Z_early_nonoverlap_2Y_z, Z_verylong_2Y_z。

## 3. 探索性反向方程

`RateVulnerability_z = gamma0 + gamma1 CAR_3day + historical_IVs + u2`

这个方程只用于展示第六讲中双向或联立系统的概念，不作为主识别。原因是 `RateVulnerability_z` 是事件前历史窗口估计得到的行业利率脆弱度，而 `CAR_3day` 是 FOMC 事件窗口的收益反应；从时间顺序看，当期 CAR 不能反过来决定历史 RV。

## 4. 结论

本文主结论仍应基于收益结构方程的 OLS-HC3、2SLS 辅助检验和稳健性结果。完整 `CAR_3day <-> RateVulnerability_z` 联立系统只作为第六讲教学性和探索性模块。
