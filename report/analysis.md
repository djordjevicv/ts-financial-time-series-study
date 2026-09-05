# BELEX price and return analysis

## Transforming prices into returns

The data begin as a vector of BELEX prices, `P`. Each price is a level, so it
depends partly on the previous price. For time-series modeling, the course
function converts consecutive prices into three types of return:

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
