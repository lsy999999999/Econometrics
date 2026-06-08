# Step 5 IV feasibility summary

## Candidate IVs tested

1. `Z_presample_2Y_z`: pre-sample 2-year rate sensitivity estimated before the main event-estimation window.
2. `Z_presample_10Y_z`: pre-sample 10-year rate sensitivity, used together with `Z_presample_2Y_z` for overidentification.
3. `Z_presample_simple_2Y_z`: pre-sample simple 2-year rate sensitivity without factor controls.
4. `Z_precovid_2Y_z`: 2017-2019 factor-adjusted 2Y rate sensitivity, added in Step 5B.
5. `Z_verylong_2Y_z`, `Z_early_nonoverlap_2Y_z`, and `Z_termspread_z`: additional distinct Step 5B IV candidates.

## Duplicate-name correction

A duplicate alias for the 2017-2019 IV was previously generated from the same window and construction. It was not a distinct IV. The code now keeps only `Z_precovid_2Y_z`.

## Relevance

- Original `Z_presample_2Y_z`: first-stage F = 5.85, partial R² = 0.128.
- Original `Z_presample_2Y_z + Z_presample_10Y_z`: joint first-stage F = 3.24, partial R² = 0.143.
- Original `Z_presample_simple_2Y_z`: first-stage F = 4.68, partial R² = 0.105.
- Step 5B `Z_precovid_2Y_z`: first-stage F = 16.42, partial R² = 0.291.
- Step 5B `Z_verylong_2Y_z`: first-stage F = 9.75, partial R² = 0.196.

`Z_precovid_2Y_z` is the strongest current IV by relevance.

## 2SLS results

- `Z_precovid_2Y_z`: beta = -1.812, p = 0.302.
- `Z_verylong_2Y_z`: beta = -1.076, p = 0.507.
- Original one-IV 2Y estimate: beta = -0.420, p = 0.792.
- Original two-IV 2Y+10Y estimate: beta = -0.861, p = 0.621.

The stronger Step 5B IVs produce negative 2SLS coefficients, but they are not statistically significant.

## DWH and overidentification diagnostics

- DWH using `Z_precovid_2Y_z`: p = 0.286, so we do not reject exogeneity of `RateVulnerability_z`.
- Overidentified distinct-IV combinations do not reject overidentifying restrictions: Sargan/Hansen-style p-values range from 0.109 to 0.237.

These diagnostics are supportive but not conclusive because N = 49 and the exclusion restriction is an economic assumption.

## Final verdict

`Z_precovid_2Y_z` is economically motivated and relevant enough for auxiliary IV evidence. It is temporally predetermined, so reverse causality from the event CAR is unlikely. However, it may still capture persistent industry risk characteristics, duration, safe-haven exposure, or macro sensitivity. Therefore IV/2SLS should be treated as exploratory robustness, not the paper's primary identification strategy. The main conclusion should rely on OLS with controls, HC3, proxy checks, measurement-error robustness, and influence diagnostics.
