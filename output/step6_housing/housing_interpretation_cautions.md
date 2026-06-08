# Housing SEM Interpretation Cautions

1. The static levels model is retained as the Lecture 6 baseline, but monthly housing variables are trending and serially correlated.
2. `months_supply` is not used in the clean IV set for the demand equation because it has a mechanical relationship with sales.
3. Newey-West standard errors address serial correlation in inference; they do not solve simultaneity or invalid instruments.
4. The preferred empirical discussion should compare levels, trend/month fixed effects, year-over-year changes, lagged IVs, and the pre-COVID sample.
5. Overidentification rejections are treated as warnings about exclusion restrictions, not as mechanical failures of the code.
6. The module is strongest as a Lecture 6 supply-demand SEM demonstration; causal claims should rely on the model-selection grid and diagnostic tables.
