# Proxy and IV distinction in Step 5

## Proxy results

Additional proxy models were added for `FinanceDummy` and `DefensiveDummy`.

Core coefficient `RateVulnerability_z` after adding each proxy:

- `GoldDummy`: beta = -0.738, p = 0.113.
- `HighGrowthProxy`: beta = -1.000, p = 0.057.
- `HighVolProxy`: beta = -1.089, p = 0.074.
- `FinanceDummy`: beta = -1.119, p = 0.047.
- `DefensiveDummy`: beta = -0.948, p = 0.060.

The coefficient remains negative across all proxy specifications. Significance weakens for some proxies, especially `GoldDummy`, which is consistent with Gold being an economically special industry rather than a bad observation.

## Mechanism interactions

- `RateVulnerability_z × HighGrowthProxy`: beta = -1.026, p = 0.345.
- `RateVulnerability_z × HighVolProxy`: beta = -1.677, p = 0.118.

The high-volatility interaction has the expected negative sign and weak directional evidence, but it is not statistically significant at conventional levels. It should be presented as suggestive mechanism evidence only.

## Why these proxies are not IVs

`GoldDummy`, `HighGrowthProxy`, `HighVolProxy`, `FinanceDummy`, and `DefensiveDummy` are not valid instruments because they plausibly affect event-window CAR directly:

- Gold has safe-haven, inflation-hedge, and dollar-exposure channels.
- Growth industries are directly sensitive to discount-rate changes.
- High-volatility industries may experience direct risk repricing during FOMC events.
- Finance has balance-sheet and net-interest-margin channels.
- Defensive industries have stable-demand and low-cyclicality channels.

Therefore they are proxy/control variables for omitted industry characteristics, not IVs satisfying the exclusion restriction.
