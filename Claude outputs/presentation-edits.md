# Slide edits — BELEX15 deck

Copy-paste text for each change. Slide numbers refer to the PDF you sent.
Skipped by request: the "by the course criterion" wording on slide 16, and the
figure styling.

---

## Slide 2 — Project idea (fix: the price forecast does not come from GARCH)

**Replace the last bullet:**

> Finally, a forecast was produced from the GARCH model and mapped back to the
> price level.

**With:**

> Finally, the AR model produced log-return forecasts, which were mapped back
> to the price level, while GARCH produced the conditional-variance forecasts.

Optionally add a one-line reminder under the bullets, since this is the
distinction the whole project rests on:

> AR gives the conditional mean, GARCH gives the conditional variance. Only the
> return forecasts are turned back into prices.

---

## Slide 4 — Data (fix: no stationarity test was run)

**Replace the second bullet:**

> The series is not stationary, so it cannot be modeled directly at the price
> level.

**With:**

> The ACF of the price level decays almost not at all (0.9996 at lag 1 to
> 0.9823 at lag 20), the signature of a non-stationary level, so the series is
> not modeled directly at the price level.

If you know them, add the sample period to the first bullet: "2,400 daily
closing prices of the BELEX15 index, [start] to [end]." It explains the 2008
spike without a word from you.

---

## Slide 5 — From prices to returns (note the index convention)

The formula `Gt = Pt+1 / Pt` matches the code but not the usual
`G_t = P_t / P_{t-1}`. Add a small caption under the formulas:

> Indexing follows `returns.m`: element k of the return vectors corresponds to
> the move from price k to price k+1, so 2,400 prices give 2,399 returns.

---

## Slide 11 — Adequacy of AR(8) (add: the test is not as weak as "1 d.f." looks)

**Add to the "The model is adequate" box, after the first line:**

> The degrees of freedom depend on how the fitted parameters are counted.
> `adequateAR` counts all 9 nonzero coefficients, leaving 1 d.f. (critical
> value 3.84). Counting only the 4 that are statistically significant leaves
> 6 d.f. (critical value 12.59). Q = 3.295 passes either way.

---

## Slide 12 — ARCH effect test (fix the overstatement, add the ARCH order)

**Replace:**

> There is no more usable structure in the mean of the residuals...

**With:**

> No significant autocorrelation is left in the first 10 residual lags...

**Add a fourth stat box, or a line under the box:**

> **ARCH(8)** — the order the PACF of the squared residuals points to

**And the supporting sentence:**

> The squared residuals follow an AR-like relation, so the PACF of `a²` gives
> the ARCH order the data support: it is significant at orders 1, 2, 3, 5, 6
> and 8, so the order is 8 (coefficient 0.0450, t = 2.19). A pure ARCH(8)
> would need nine parameters; GARCH(1,1) represents the same persistence with
> three, which is why it is fitted next.

This is the answer to "why GARCH(1,1)?", which the deck currently assumes.

---

## Slide 17 — Forecast table (label what the intervals are)

**Replace the third bullet:**

> The intervals use a Gaussian approximation: Normality was rejected and GARCH
> is marginally inadequate, so the intervals are approximate.

**With:**

> The intervals use a Gaussian approximation and the conditional variance of
> the return at each step, so they are one-step intervals evaluated at every
> horizon rather than full h-step forecast intervals. Normality was rejected
> and GARCH is marginally inadequate, so they are indicative.

---

## Slides 10 and 14 — merge the two code slides

Two full code dumps are a lot of airtime. Merge into one "Implementation"
slide and use the freed slot for the ARCH-order evidence above.

**Title:** Implementation in MATLAB

**Left panel — `ARLScoef.m` (AR, closed form):**

```matlab
function res=ARLScoef(r,p)
M=tmatrix(r,p);          % lags + column of ones
y=r(p+1:end,1);
x=linsolve(M'*M,M'*y);   % normal equations
res=x;
```

**Right panel — `GARCHcoef.m` (GARCH, constrained optimization):**

```matlab
function res=GARCHcoef(a,m,s)
x0=[zeros(m,1);0.1;zeros(s,1)];
A=[ones(1,m) 0 ones(1,s)];  b=1-eps;   % alpha+beta < 1
lb=[zeros(m,1);eps;zeros(s,1)];        % non-negative
[x,f]=fmincon(@(x)fmlGARCH(a,x,m,s),x0,A,b,[],[],lb,[]);
res=x;
```

**Two bullets:**

> - The AR coefficients have a closed form, so they come straight from the
>   normal equations, with no numerical optimization.
> - GARCH has no closed form: `fmincon` minimizes the negative log-likelihood
>   from `fmlGARCH.m` under the stationarity constraint and non-negativity.

---

## Backup slide — questions worth having an answer to

Not for the main flow; keep it after "Thank you".

- **Why keep lags 4 to 7 when they are insignificant?** The PACF criterion
  tests only the highest lag. A subset model (lags 1, 2, 3, 8 only) is the
  natural alternative and was not fitted.
- **What would you do about GARCH failing the adequacy test?** The failure is
  in the level of the standardized residuals, not their squares, so the mean
  model is the place to look — lags 13 to 17 are outside AR(8)'s reach.
  t-distributed innovations or EGARCH would address the tails and asymmetry
  instead.
- **Is the 13 to 17 lag pattern seasonal?** It is roughly a three-week cycle.
  The project did not test for seasonality.
- **How was stationarity of the returns established?** It was not tested
  formally; the evidence is the ACF comparison in Stage 2.
