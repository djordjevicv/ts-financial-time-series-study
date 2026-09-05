# Real Financial Time-Series Analysis Using BELEX Data

This MATLAB project applies the Time Series course methods to real BELEX price data.

## Repository structure

```text
data/                  BELEX price data
scripts/               Analysis stages
src/
  returns/             Return calculations
  diagnostics/         ACF and descriptive-statistics functions
results/
  figures/             Generated plots
  tables/              Generated tables
run_pipeline.m         Project entry point
```

## Data

`data/belex.mat` is the source for the analysis.

The files under `src/` are the original course implementations.

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

## Run

Open MATLAB, set the Current Folder to the repository root, and run:

```matlab
run_pipeline
```
