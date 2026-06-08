# Step 5B additional IV search summary

## Duplicate-name audit

A duplicate alias for the 2017-2019 IV was identified and removed. The code now keeps a single name: `Z_precovid_2Y_z`.

## What was added

Additional candidate IVs were constructed from non-overlapping historical rate-sensitivity windows:

- `Z_recent_nonoverlap_2Y_z`: event -520 to -261 trading days.
- `Z_early_nonoverlap_2Y_z`: event -760 to -521 trading days.
- `Z_verylong_2Y_z`: event -1500 to -521 trading days.
- `Z_precovid_2Y_z`: 2017-01-03 to 2019-12-31.
- `Z_termspread_z`: sensitivity to changes in 10Y-2Y term spread.
- `RV_early_2Y_z` and `RV_late_2Y_z`: split-sample measurement-error IV candidates.

## First-stage ranking

The strongest candidates are:

1. `Z_precovid_2Y_z`: F = 16.42, partial R² = 0.291.
2. `Z_verylong_2Y_z`: F = 9.75, partial R² = 0.196.
3. `Z_early_nonoverlap_2Y_z`: F = 5.85, partial R² = 0.128.

`Z_recent_nonoverlap_2Y_z` is very weak: F = 0.08. The split-sample IVs are also weak: F = 0.28.

## 2SLS results

- `Z_precovid_2Y_z`: beta = -1.812, p = 0.302.
- `Z_verylong_2Y_z`: beta = -1.076, p = 0.507.
- `Z_early_nonoverlap_2Y_z`: beta = -0.420, p = 0.792.

The strongest IVs produce negative 2SLS coefficients, consistent with the OLS direction, but the estimates are not statistically significant.

## Durbin-Wu-Hausman check for strongest IV

Using residual-inclusion DWH with `Z_precovid_2Y_z`, the residual coefficient p-value is 0.286. The test does not reject the null that `RateVulnerability_z` can be treated as exogenous in the M2 equation. This supports keeping OLS-HC3 as the main specification. It does not prove exogeneity, especially with N = 49.

## Overidentification tests

Because `Z_precovid_2Y_z` is a single IV for one endogenous variable, it is exactly identified by itself and cannot be used for an overidentification test alone.

For genuinely distinct IV combinations, the following tests were run:

1. `Z_precovid_2Y_z + Z_verylong_2Y_z`:
   - Sargan p = 0.172;
   - Hansen-style J p = 0.208.
2. `Z_precovid_2Y_z + Z_early_nonoverlap_2Y_z`:
   - Sargan p = 0.159;
   - Hansen-style J p = 0.237.
3. `Z_precovid_2Y_z + Z_termspread_z`:
   - Sargan p = 0.138;
   - Hansen-style J p = 0.109.

These tests do not reject the null that the overidentifying restrictions are valid. However, this should not be interpreted as proof of IV exogeneity because N = 49 is small and exclusion restrictions are primarily economic assumptions.

## Verdict

`Z_precovid_2Y_z` is the best available IV candidate by relevance, with first-stage F above 10. It can be reported as auxiliary IV evidence. However, its exclusion restriction is still contestable because 2017-2019 historical rate sensitivity may capture persistent industry duration, macro-risk, or safe-haven characteristics that also affect 2022 FOMC-window returns. The IV section should state that stronger auxiliary IV relevance was found, but 2SLS remains imprecise and should not replace the main OLS-HC3 identification strategy.
