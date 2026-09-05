# BELEX financial time-series study

This project uses MATLAB and methods from the Time Series course to analyze
real BELEX price data. The full analysis converts prices to returns, fits an
autoregressive (AR) model for the conditional mean, and checks the AR residuals
for ARCH effects. If the data support a GARCH model, the study also estimates
conditional variance. The final stages forecast returns and volatility, then
reconstruct a price path from the return forecasts.

The current pipeline covers three parts of the analysis:

1. loading and validating BELEX prices, then calculating gross, net, and log
   returns;
2. comparing autocorrelation in prices and log returns; and
3. calculating descriptive statistics and Normality tests for log returns.

## Repository structure

```text
data/                  BELEX price data
scripts/               Implemented analysis
src/
  returns/             Return calculations
  diagnostics/         ACF and descriptive-statistics functions
results/
  figures/             Generated plots
  tables/              Generated tables and interpretations
report/analysis.md     Analysis
run_pipeline.m         Project entry point
```

The source data are in `data/belex.mat`.

## Prerequisites

- MATLAB
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

## Run the pipeline

In MATLAB, set the Current Folder to the repository root and enter:

```matlab
run_pipeline
```
