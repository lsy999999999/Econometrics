# Step 6 lecture coverage check

This module maps Lecture 6, simultaneous equations models, into the empirical project while keeping the economic identification limits explicit.

Covered topics:
1. Simultaneous-equation motivation: `CAR_3day` and `RateVulnerability_z` are placed in a two-equation exploratory system.
2. Structural equations: return equation and exploratory vulnerability equation are both estimated.
3. Structural coefficients: key coefficients are saved in `step6_model_comparison.csv`.
4. Structural errors: residual correlation is saved in `step6_system_residual_correlation.csv`.
5. Endogenous variables: `CAR_3day` and `RateVulnerability_z` are listed in `step6_sem_variable_roles.csv`.
6. Exogenous and predetermined variables: historical rate-exposure IVs and factor controls are classified in the variable-role table.
7. Predetermined variables: pre-event windows such as `Z_precovid_2Y_z` are treated as predetermined external shifters.
8. OLS inconsistency under simultaneity: OLS is reported only as a baseline comparison.
9. Recursive-system logic: `step6_recursive_check.csv` notes that a valid recursive system requires one-way ordering and uncorrelated errors.
10. ILS: `step6_ils_reduced_form_demo.csv` gives a just-identified reduced-form/Wald-ratio teaching example.
11. 2SLS: manual HC3 2SLS results are saved for both equations.
12. System estimation: `systemfit` 2SLS and optional 3SLS outputs are saved when available.
13. Order condition: `step6_identification_conditions.csv` classifies under/just/over identification by excluded exogenous variables.
14. Rank-condition intuition: `step6_first_stage_relevance.csv` reports first-stage F and partial R-squared as empirical relevance checks.
15. Overidentification: Sargan and Hansen-style J tests are saved in `step6_overid_tests.csv`.
16. 2SLS fit-statistic warning: `step6_2sls_fit_warning.csv` records the Lecture 6 warning that ordinary R-squared/F are not reliable for 2SLS.

Important limitation:
The reverse equation is not used as main evidence because same-event CAR cannot literally determine a historical pre-event RV measure, and the excluded controls for CAR plausibly affect RV directly.
