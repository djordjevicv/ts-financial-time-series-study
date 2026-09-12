# Changes on top of `main` (after PR #8)

`main` already contains Stages 4–7 as merged in PR #8. Everything below is the
uncommitted work that follows from checking the implementation against the
course script.

## Headline: no result changed

Every estimate, test statistic and forecast is identical to `main`:

- AR(8), same coefficients
- ARCH test `Q = 831.49` against `15.51`
- GARCH(1,1): `alpha_0 = 4.97e-06`, `alpha_1 = 0.2439`, `beta_1 = 0.7461`
- All 20 forecast rows

`forecasts.csv`, `garch_summary.csv`, `garch_coefficients.csv` and
`ar_coefficients.csv` are byte-identical. Nothing under `src/` changed in this
round. The changes add evidence, a new diagnostic and clearer wording.

## 1. AR order selection is now documented as the PACF criterion

**Files:** `scripts/ar_model.m`, `results/tables/ar_order_selection.csv`,
`README.md`, `report/analysis.md`

`PACFAR` scans every candidate order up to `ceil(log(T)) = 8` and keeps the
largest one with a significant highest-lag coefficient, rather than stopping
at the first insignificant one. That is not an oversight: the criterion is
`PACF(p) != 0` and `PACF(p+s) = 0` for every `s > 0`, and on BELEX the PACF is
significant again at lag 8 after four insignificant orders. A rule that stopped
early would pick AR(3), which does not satisfy the second half of the
criterion.

We tested that alternative before rejecting it. AR(3) fails the adequacy test
(`Q = 16.88` against `12.59`), so it would not have survived the model check
either.

What changed:

- The comment block states the criterion instead of just describing the loop.
- The order a stopping rule would give is computed and reported, so the
  alternative is on the record: "Stopping instead at the first insignificant
  highest-lag coefficient, AR(4) (t = 1.734), would have given AR(3)."
- `ar_order_selection.csv` gains a `Selected` column.

## 2. Adequacy degrees of freedom are reported both ways

**Files:** `scripts/ar_model.m`, `results/tables/ar_adequacy.csv`

The Ljung-Box test uses `chi2` with `m - g` degrees of freedom, where `g` is
the number of parameters statistically different from zero. `adequateAR`
instead counts every nonzero coefficient, which for AR(8) means 9 instead of 4
and leaves just 1 degree of freedom.

`adequateAR` is unchanged — it is the supplied implementation and its decision
is still the result. We only report both counts, so the adequacy of AR(8) does
not appear to hang on a single degree of freedom:

| Count of fitted parameters | g | df | Critical value | Decision |
| --- | ---: | ---: | ---: | --- |
| `adequateAR`: all nonzero coefficients | 9 | 1 | 3.84 | adequate |
| Parameters significantly different from zero | 4 | 6 | 12.59 | adequate |

`Q = 3.295` passes both. `ar_adequacy.csv` gains four columns:
`SignificantParameters`, `DegreesOfFreedomSignificantOnly`,
`CriticalValueSignificantOnly`, `AdequateSignificantOnly`.

## 3. The ARCH order is now estimated, which justifies GARCH(1,1)

**Files:** `scripts/arch_effect.m`, new `results/tables/arch_order_pacf.csv`

The squared residuals follow an AR-like relation, so the PACF of `a.^2` gives
the ARCH order the data support. We previously jumped straight to GARCH(1,1)
as a project choice with no evidence behind it.

`PACFAR(a.^2)` points to **ARCH(8)** — the highest-lag coefficient is
significant at orders 1, 2, 3, 5, 6 and 8. A pure ARCH(8) needs nine
parameters; GARCH(1,1) represents the same persistence with three. The fitted
model does not change; the choice is now argued rather than asserted.

## 4. The forecast intervals say what they actually are

**File:** `scripts/forecasting.m`

The intervals `rForecast +/- 1.96*sqrt(sigma2Forecast)` use the conditional
variance of the return at each step. A true h-step interval accumulates the
forecast error of the intermediate steps (`D(e(1)) = sigma_a^2`,
`D(e(2)) = (1 + phi_1^2) sigma_a^2`, and so on), so ours is exact at `h = 1`
and too narrow afterwards.

The arithmetic is unchanged, because the project notes name this formula. The
generated text now states the limitation. A second sentence was added that
prints only if the AR model fails its adequacy test (it does not at present).

## 5. Wording

**Files:** `scripts/ar_model.m`, `scripts/garch_model.m`, `README.md`,
`report/analysis.md`

References to the course and the script were removed from the repository. The
statements stand on their own: "by this criterion", "which is how g is
defined", "the supplied time-series functions".

## Files touched

| File | Change |
| --- | --- |
| `scripts/ar_model.m` | Criterion comment, stopping-rule comparison, `Selected` column, both degrees-of-freedom counts |
| `scripts/arch_effect.m` | PACF of squared residuals, new table, ARCH-order sentence |
| `scripts/forecasting.m` | Interval caveat, AR-adequacy caveat |
| `scripts/garch_model.m` | Wording only |
| `README.md`, `report/analysis.md` | Stage 5 description, new outputs, both df counts, ARCH-order paragraph, interval caveat, wording |
| `results/tables/arch_order_pacf.csv` | New file — needs `git add` |
| `results/tables/*` (others) | Regenerated; only the new columns differ |
| `results/figures/*.png` | Re-rendered by the latest run; identical plots |

## Verification

`run_pipeline` was run end to end in MATLAB R2025b after every change, with no
errors or warnings. No file under `src/` was modified in this round.
