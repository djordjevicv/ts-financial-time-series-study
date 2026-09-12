# BELEX price and return analysis

## Transforming prices into returns

The data begin as a vector of BELEX prices, `P`. Each price is a level, so it
depends partly on the previous price. For time-series modeling, `returns.m`
converts consecutive prices into three types of return:

\[
G_t = \frac{P_{t+1}}{P_t}, \qquad
R_t = G_t - 1, \qquad
r_t = \log(G_t).
\]

`Gt` contains gross returns, `Rt` contains net returns, and `rt` contains
continuously compounded log returns. Returns answer a different question from
prices: instead of showing the level of the index, they show its change between
two observations. The later AR model will use `rt`. The price vector `P` is
kept for comparison and for reconstructing forecast prices.

## Autocorrelation of prices and returns

The sample ACF was calculated for the first 20 lags of both series. Every price
lag is statistically significant. The autocorrelation falls only slightly,
from about `0.999580` at lag 1 to `0.982316` at lag 20, which shows how strongly
the price level depends on its recent history.

Seven of the 20 log-return lags are significant:

```text
[1 2 13 14 15 16 17]
```

Changing from prices to log returns reduces the number of significant lags
from 20 to 7. Much of the persistence in the price series disappears, although
the significant return lags show that some serial dependence remains.

The ACF is useful evidence when choosing an AR order, but it does not choose
the order by itself. In particular, seven significant lags do not imply
`p = 7`. The AR stage must select and check the order separately, using `rt` as
the model input.

## Descriptive statistics of log returns

The Stage 3 results for `rt` are:

| Statistic | Value |
| --- | ---: |
| Mean | -0.000130650 |
| Variance | 0.000190575 |
| Standard deviation | 0.0138049 |
| Skewness | 0.134847 |
| Kurtosis | 16.6398 |
| Excess kurtosis | 13.6398 |
| Jarque-Bera statistic | 18604.0 |
| Jarque-Bera critical value | 5.99 |

The mean is close to zero when compared with the standard deviation. Positive
skewness points to a slight right asymmetry, and the skewness test rejects
symmetry at the 5% level. Kurtosis is `16.6398`, far above the Normal value of
3. The log returns are therefore leptokurtic, with much heavier tails than a
Normal distribution.

The Jarque-Bera statistic of `18604.0` is well above the critical value of
`5.99`, so the Normality test is rejected at the 5% level. The transformation
removes most of the persistence seen in prices, but the log returns still have
detectable autocorrelation and a clearly non-Normal distribution.

## AR model for the conditional mean

The AR model is fitted to the log returns `rt`. `PACFAR` fits AR(1), AR(2),
... by least squares and reads the PACF as the highest-lag coefficient of each
fit. With 2399 returns it tests up to AR(8), because `ceil(log(2399)) = 8`.
The criterion for the order is that the PACF is nonzero at the order and
zero at every larger lag, so the order is the largest candidate with a
significant highest-lag coefficient:

| Candidate order | Highest-lag coefficient | t-ratio | Significant |
| ---: | ---: | ---: | :---: |
| 1 | 0.322197 | 16.66 | yes |
| 2 | 0.048658 | 2.38 | yes |
| 3 | -0.050432 | -2.47 | yes |
| 4 | 0.035472 | 1.73 | no |
| 5 | -0.015019 | -0.73 | no |
| 6 | 0.016963 | 0.83 | no |
| 7 | 0.023915 | 1.17 | no |
| 8 | 0.056899 | 2.78 | **yes, the largest** |

The selected model is therefore **AR(8)**. The PACF is significant again at
lag 8 after four insignificant orders, so a rule that stopped at the first
insignificant order would have given AR(3). That AR(3) does not satisfy the
criterion, because its PACF is not zero at every larger lag. The selection
table is saved in `ar_order_selection.csv` with both orders marked.

The least-squares estimates of AR(8) are:

| Term | Estimate | Standard error | t-ratio | Significant |
| --- | ---: | ---: | ---: | :---: |
| `phi_1` | 0.309607 | 0.020455 | 15.14 | yes |
| `phi_2` | 0.060085 | 0.021421 | 2.80 | yes |
| `phi_3` | -0.059407 | 0.021459 | -2.77 | yes |
| `phi_4` | 0.038296 | 0.021491 | 1.78 | no |
| `phi_5` | -0.018381 | 0.021492 | -0.86 | no |
| `phi_6` | 0.006055 | 0.021462 | 0.28 | no |
| `phi_7` | 0.006375 | 0.021442 | 0.30 | no |
| `phi_8` | 0.056899 | 0.020482 | 2.78 | yes |
| `phi_0` | -0.0000894 | 0.0002667 | -0.34 | no |

The first lag dominates: about 31% of the latest log return carries over to the
next conditional mean. Lags 2, 3, and 8 are also significant, while lags 4 to
7 are not. The PACF criterion tests only the highest lag, so the full AR(8) is
kept and the insignificant coefficients are reported rather than removed. The
constant is not significant, which matches the Stage 3 finding that the mean
return is close to zero.

`adequateAR` tests the residuals `a` with a Ljung-Box statistic over 10 lags.
The statistic is `Q = 3.295`, below the critical value, so **AR(8) is
adequate** for the conditional mean. The degrees of freedom depend on how the
fitted parameters are counted:

| Count of fitted parameters | g | Degrees of freedom | Critical value | Decision |
| --- | ---: | ---: | ---: | --- |
| `adequateAR`: all nonzero coefficients | 9 | 1 | 3.84 | adequate |
| Parameters significantly different from zero | 4 | 6 | 12.59 | adequate |

Both counts give the same decision, and counting only the significant
parameters leaves a wide margin. The test covers only the first 10 residual
lags. Individually, one of the first 20 residual ACF lags (lag 14) is
significant, and the log-return ACF lags 13 to 17 from Stage 2 lie beyond the
largest order that `PACFAR` considers. The AR(8) model captures the short-lag
dependence but not that longer-lag pattern.

## ARCH effect in the AR residuals

`ARCHeffect` applies the Ljung-Box test to `a.^2` with 8 lags. The statistic
is `Q = 831.49`, far above the critical value `15.51`, so the ARCH effect is
statistically significant. Twelve of the first 20 squared-residual ACF lags are
individually significant (lags 1 to 9, 14, 17, and 20). Large residuals tend to
follow large residuals, which is the volatility clustering visible in the
log-return plot. Volatility clustering is also consistent with the heavy tails
found in Stage 3. Because the effect is significant, the pipeline models the
conditional variance.

The squared residuals satisfy an AR-like relation, `a2(t) ~ alpha_0 +
sum(alpha_i a2(t-i))`, so the PACF of `a.^2` indicates the ARCH order the data
support. Applied to the squared AR(8) residuals it points to **ARCH(8)**: the
highest-lag coefficient is significant at orders 1, 2, 3, 5, 6 and 8, the
largest being 8 (`0.0450`, t = 2.19). A pure ARCH(8) would need nine
parameters. GARCH(1,1) represents the same persistence with three, which is
why it is the model fitted in Stage 6. The candidate orders are in
`arch_order_pacf.csv`.

## GARCH(1,1) conditional variance

GARCH(1,1) was fitted to the AR(8) residuals with `GARCHcoef`:

\[
\sigma_t^2 = \alpha_0 + \alpha_1 a_{t-1}^2 + \beta_1 \sigma_{t-1}^2 .
\]

| Parameter | Estimate |
| --- | ---: |
| `alpha_0` | 4.97277e-06 |
| `alpha_1` | 0.243878 |
| `beta_1` | 0.746063 |

`GARCHcoef` starts `fmincon` at `alpha_0 = 0.1`, several orders of magnitude
above the residual variance, so the script also refits the model on the
percent residuals `100*a` and maps `alpha_0` back. The GARCH likelihood is
unchanged by this rescaling apart from a constant. The direct fit reached a
negative log-likelihood of `-9827.1863`, and the rescaled fit reached
`-9827.1871`. The difference is within tolerance, so the direct `GARCHcoef`
estimate is kept. The two fits agree to about three significant digits.

The persistence is `alpha_1 + beta_1 = 0.9899`. A volatility shock therefore
decays slowly, with a half-life of about 69 observations. The implied
unconditional variance is `alpha_0 / (1 - alpha_1 - beta_1) = 4.94e-4`, about
2.9 times the sample variance of the residuals (`1.69e-4`). When persistence
is this close to 1, the ratio is very sensitive to small changes in the
estimates, and the long-run level reflects the high-volatility episodes in the
sample. Treat it as an imprecise figure.

`adequateGARCH` applies the Ljung-Box test to the standardized residuals
`a./sqrt(s2)` with 8 lags. The statistic is `Q = 15.582`, just above the
critical value `15.51`, so GARCH(1,1) is **not adequate** by this
criterion. The margin is small (about 0.5%). As a supplementary check, the
script applies `ARCHeffect` to the standardized residuals. It gives
`Q = 11.99`, below `15.51`, compared with `831.49` before the GARCH fit. The
GARCH model removes the ARCH effect; the marginal failure comes from
autocorrelation in the level of the standardized residuals, not from remaining
volatility clustering. This is consistent with the longer-lag dependence that
the AR(8) mean does not capture. The variance forecasts below are still
reported, with this caveat.

## Forecasts

The pipeline forecasts `k = 20` observations from the last observed price,
`729.84`. The AR model supplies conditional-mean log returns, and GARCH
supplies conditional variances. The price path is reconstructed with `task4`
from the AR log-return forecasts only.

| h | AR log return | GARCH variance | GARCH std. dev. | Gaussian 95% return interval | Price |
| ---: | ---: | ---: | ---: | :---: | ---: |
| 1 | -0.002865 | 1.365e-4 | 0.01168 | [-0.02576, 0.02003] | 727.75 |
| 2 | 0.001304 | 1.401e-4 | 0.01184 | [-0.02189, 0.02450] | 728.70 |
| 5 | 0.000968 | 1.507e-4 | 0.01227 | [-0.02309, 0.02503] | 730.49 |
| 10 | -0.000164 | 1.676e-4 | 0.01295 | [-0.02554, 0.02521] | 729.41 |
| 20 | -0.000143 | 1.990e-4 | 0.01411 | [-0.02779, 0.02751] | 728.46 |

The first return forecast is negative (`-0.00286`) because the last two
observed returns were negative and `phi_1` carries them forward. The forecasts
then oscillate for a few steps and settle near the AR long-run mean,
`phi_0 / (1 - sum(phi_i)) = -0.000149`. The reconstructed price moves from
`727.75` at `h = 1` to a high of `731.06` at `h = 6`, then ends at `728.46` at
`h = 20`. The cumulative change is -0.19%. The path compounds conditional-mean
log returns, so it is a central path, not a price interval. The mean forecasts
are small next to the forecast standard deviation: about a quarter of it at
`h = 1` and roughly a tenth or less afterward. The AR model offers little
directional information beyond the first few steps.

The GARCH conditional variance at `h = 1` is `1.365e-4`, a standard deviation
of 1.17% per observation. This is slightly below the full-sample residual
standard deviation of 1.30%, because the last residual was small. The forecast
variance then rises steadily toward the unconditional level, reaching
`1.990e-4` at `h = 20`. Given the 69-observation half-life, it is still far
from that level. The return intervals use
`rForecast +/- 1.96*sqrt(sigma2Forecast)` and are a **Gaussian conditional
approximation**. Each one uses only the conditional variance of the return at
that step, so it is a one-step interval evaluated at every horizon: it does
not accumulate the forecast error of the steps in between, and for `h > 1` it
is narrower than a full h-step forecast interval. Stage 3 also rejected
Normality, and the GARCH model narrowly failed its adequacy test, so these
bounds are indicative rather than exact.

## Summary

- Log returns remove most of the price persistence. Short-lag dependence
  remains and is modeled by AR(8), the order the PACF criterion selects. It
  passes the adequacy test for the conditional mean under either way of
  counting the fitted parameters.
- The squared AR residuals show a strong ARCH effect (`Q = 831.49` against
  `15.51`).
- GARCH(1,1) captures this effect: no significant ARCH effect remains in the
  standardized residuals. Volatility is highly persistent (`0.9899`). The model
  narrowly fails the adequacy test on the standardized residuals
  (`Q = 15.582` against `15.51`).
- AR forecasts give conditional-mean log returns and a nearly flat
  reconstructed price path. GARCH forecasts give conditional variances that
  rise slowly toward the long-run level.
