# BELEX financial time-series study

This project analyzes real BELEX price data in MATLAB. The analysis converts
prices to returns, fits an autoregressive (AR) model for the conditional mean,
and checks the AR residuals for ARCH effects. When the ARCH effect is
significant, a GARCH model describes the conditional variance. The final stage
forecasts log returns and conditional variances, then reconstructs a price
path from the return forecasts.

The pipeline runs seven stages:

1. load and validate BELEX prices, then calculate gross, net, and log returns;
2. compare autocorrelation in prices and log returns;
3. calculate descriptive statistics and Normality tests for log returns;
4. select, estimate, and check an AR model for log returns;
5. test the squared AR residuals for an ARCH effect and read their PACF for
   the ARCH order they support;
6. fit GARCH(1,1) to the AR residuals when the ARCH effect is significant; and
7. forecast log returns (AR), conditional variances (GARCH), and prices.

The two models answer different questions:

```text
AR    -> conditional mean     -> log-return and price forecasts
GARCH -> conditional variance -> variance forecasts and return intervals
```

## Repository structure

```text
data/                  BELEX price data
scripts/               One script per stage
src/
  returns/             Return calculations and price reconstruction
  diagnostics/         ACF and descriptive-statistics functions
  ar/                  AR estimation, order selection, adequacy, forecasts
  garch/               ARCH test, GARCH estimation, adequacy, forecasts
results/
  figures/             Generated plots
  tables/              Generated tables and interpretations
report/analysis.md     Analysis
run_pipeline.m         Project entry point
```

The source data are in `data/belex.mat`. The functions under `src/` are the
supplied time-series functions. `ARestimate.m` and `statsig.m` compute the Normal
critical value with `sqrt(2)*erfcinv(0.05)` instead of `-norminv(0.025)`, the
same equivalent form already used in `basicstat1.m` and `tratio2.m`, so the
Statistics and Machine Learning Toolbox is not needed.

Stage 4 selects the AR order with `PACFAR`, which follows the PACF
criterion: the order is the largest candidate up to `ceil(log(T))` whose
highest-lag coefficient is significant. The order that a
stop-at-the-first-insignificant rule would give is reported alongside it for
comparison, and `ar_adequacy.csv` reports the Ljung-Box degrees of freedom
both as `adequateAR` computes them and with only the significant parameters
counted.

## Prerequisites

- MATLAB (tested with R2025b)
- Optimization Toolbox (`GARCHcoef` uses `fmincon`)

## Run the pipeline

In MATLAB, set the Current Folder to the repository root and enter:

```matlab
run_pipeline
```

`scripts/garch_model.m` runs only when the ARCH test is significant. Without an
ARCH effect, the pipeline still produces the AR return and price forecasts and
reports that conditional-variance modeling was not supported.

## Outputs

| Stage | Tables (`results/tables/`) | Figures (`results/figures/`) |
| --- | --- | --- |
| 1 | | `belex_prices.png`, `belex_log_returns.png` |
| 2 | `acf_significance.csv`, `acf_interpretation.txt` | `price_acf.png`, `log_return_acf.png` |
| 3 | `basic_statistics_rt.csv`, `basic_statistics_rt_interpretation.txt` | |
| 4 | `ar_order_selection.csv`, `ar_coefficients.csv`, `ar_adequacy.csv`, `ar_model_interpretation.txt` | `ar_residuals.png`, `ar_residual_acf.png` |
| 5 | `arch_effect.csv`, `arch_order_pacf.csv`, `arch_effect_interpretation.txt` | `ar_squared_residual_acf.png` |
| 6 | `garch_coefficients.csv`, `garch_summary.csv`, `garch_model_interpretation.txt` | `garch_conditional_variance.png` |
| 7 | `forecasts.csv`, `forecast_interpretation.txt` | `forecast_prices.png`, `forecast_log_returns.png`, `forecast_conditional_variance.png` |
