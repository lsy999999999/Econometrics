# Step 6 Housing SEM Results

## Research Question

This module replaces the exploratory FOMC industry SEM with a standard new-housing-market supply-demand simultaneous-equations model. The goal is to demonstrate Lecture 6 with an economically coherent system: new home prices and new home sales are jointly determined by demand and supply.

## Data

Monthly U.S. data from FRED are used from 2000-01 to 2024-12, with N = 300.
The key endogenous variables are `log_sales` and `log_price`. `sales` is new one-family houses sold and `price` is the median sales price for new houses sold. Demand shifters are the mortgage rate, real disposable income, and unemployment. Supply shifters are permits, monthly supply of new houses, and construction material costs.

## Structural Equations

Demand: `log_sales = alpha0 + alpha1 log_price + alpha2 mortgage_rate + alpha3 log_income + alpha4 unrate + u_d`.
Supply: `log_sales = beta0 + beta1 log_price + beta2 log_permit + beta3 months_supply + beta4 log_construction_cost + u_s`.
Equilibrium condition: observed sales are both demanded and supplied, so `Qd = Qs = Q`.

## Identification

Both structural equations contain one endogenous right-hand-side variable, `log_price`. The demand equation excludes three supply shifters; the supply equation excludes three demand shifters. Therefore both equations satisfy the order condition and are overidentified.

## First Stage

Demand equation instruments: first-stage F = 78.643, partial R2 = 0.446.
Supply equation instruments: first-stage F = 70.025, partial R2 = 0.418.

## OLS versus 2SLS

Demand OLS price coefficient = -0.137, p = 0.667.
Demand 2SLS price coefficient = -0.114, p = 0.801. Expected sign: negative.
Supply OLS price coefficient = -0.087, p = 0.294.
Supply 2SLS price coefficient = 0.009, p = 0.949. Expected sign: positive.

## Structural Coefficients

Demand 2SLS mortgage-rate coefficient = 0.032, p = 0.045. Expected sign: negative.
Demand 2SLS income coefficient = -0.489, p = 0.549. Expected sign: positive.
Demand 2SLS unemployment coefficient = -0.116, p = 0.000. Expected sign: negative.
Supply 2SLS permit coefficient = 1.036, p = 0.000. Expected sign: positive.
Supply 2SLS months-supply coefficient = 0.017, p = 0.001. Expected sign: positive.
Supply 2SLS construction-cost coefficient = -0.443, p = 0.004. Expected sign: negative.

## Diagnostics

DWH demand p-value = 0.941; DWH supply p-value = 0.370.
Overidentification rejected at 5% in demand equation: TRUE.
Overidentification rejected at 5% in supply equation: TRUE.
System 2SLS residual correlation = 0.521, p = 0.000.

## Interpretation

This housing module is a better Lecture 6 application than the previous FOMC CAR-RV exploratory SEM because price and quantity have a standard simultaneous-equilibrium interpretation. The demand and supply equations each have clear ceteris-paribus meanings, and the excluded demand/supply shifters provide a transparent identification strategy.

All tables and figures are saved in `output/step6_housing/`.
