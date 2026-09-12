% Forecast log returns (AR), conditional variances (GARCH), and prices.

if ~exist('P', 'var') || ~exist('rt', 'var') || ~exist('xAR', 'var') || ...
        ~exist('p', 'var') || ~exist('a', 'var') || ~exist('hasARCH', 'var')
    error('BELEX:MissingModel', ...
        'Run scripts/ar_model.m and scripts/arch_effect.m before this script.');
end

if hasARCH && ~exist('xGARCH', 'var')
    error('BELEX:MissingGARCH', ...
        'The ARCH effect is significant; run scripts/garch_model.m first.');
end

% Resolve output folders from this file's location, independent of MATLAB's
% current working directory.
scriptPath = mfilename('fullpath');
repositoryRoot = fileparts(fileparts(scriptPath));
figuresDirectory = fullfile(repositoryRoot, 'results', 'figures');
tablesDirectory = fullfile(repositoryRoot, 'results', 'tables');
if exist(figuresDirectory, 'dir') ~= 7
    mkdir(figuresDirectory);
end
if exist(tablesDirectory, 'dir') ~= 7
    mkdir(tablesDirectory);
end

% Forecast horizon, in observations.
k = 20;

% AR -> conditional mean of the log return.
rForecast = kpredictionsAR(rt, xAR, k);
rForecast = rForecast(:);

% Prices are reconstructed from the predicted log returns only:
% P(T+h) = P(T+h-1) * exp(rForecast(h)).
priceForecast = task4(P(end), rForecast);
priceForecast = priceForecast(:);

if numel(rForecast) ~= k || any(~isfinite(rForecast)) || ...
        numel(priceForecast) ~= k + 1 || any(~isfinite(priceForecast))
    error('BELEX:InvalidForecast', ...
        'The AR and price forecasts must contain finite values for %d steps.', k);
end

% GARCH -> conditional variance of the log return. It is never used to
% build the price path.
if hasARCH
    sigma2Forecast = kpredictionsGARCH(a, xGARCH, mG, sG, k);
    sigma2Forecast = sigma2Forecast(:);
    if numel(sigma2Forecast) ~= k || any(~isfinite(sigma2Forecast)) || ...
            any(sigma2Forecast <= 0)
        error('BELEX:InvalidVarianceForecast', ...
            'The GARCH variance forecast must be finite and positive.');
    end
    % Gaussian conditional approximation of a 95% log-return interval.
    returnLower95 = rForecast - 1.96 * sqrt(sigma2Forecast);
    returnUpper95 = rForecast + 1.96 * sqrt(sigma2Forecast);
else
    sigma2Forecast = [];
    returnLower95 = [];
    returnUpper95 = [];
    % Remove GARCH outputs left by an earlier run so results/ matches this run.
    staleGARCHOutputs = { ...
        fullfile(figuresDirectory, 'garch_conditional_variance.png'), ...
        fullfile(figuresDirectory, 'forecast_conditional_variance.png'), ...
        fullfile(tablesDirectory, 'garch_coefficients.csv'), ...
        fullfile(tablesDirectory, 'garch_summary.csv'), ...
        fullfile(tablesDirectory, 'garch_model_interpretation.txt')};
    for staleIndex = 1:numel(staleGARCHOutputs)
        if exist(staleGARCHOutputs{staleIndex}, 'file') == 2
            delete(staleGARCHOutputs{staleIndex});
        end
    end
    clear staleGARCHOutputs staleIndex;
end

nPrices = numel(P);
Horizon = (1:k)';
Observation = nPrices + Horizon;
LogReturnForecast = rForecast;
PriceForecast = priceForecast(2:end);
if hasARCH
    ConditionalVarianceForecast = sigma2Forecast;
    ConditionalStdForecast = sqrt(sigma2Forecast);
    GaussianLower95 = returnLower95;
    GaussianUpper95 = returnUpper95;
    forecastTable = table(Horizon, Observation, LogReturnForecast, ...
        ConditionalVarianceForecast, ConditionalStdForecast, ...
        GaussianLower95, GaussianUpper95, PriceForecast);
else
    forecastTable = table(Horizon, Observation, LogReturnForecast, ...
        PriceForecast);
end
writetable(forecastTable, fullfile(tablesDirectory, 'forecasts.csv'));

% Long-run level implied by the AR model: phi_0 / (1 - sum(phi_i)).
arLagSum = sum(xAR(1:p));
if abs(1 - arLagSum) > eps
    arLongRunMean = xAR(p + 1) / (1 - arLagSum);
else
    arLongRunMean = NaN;
end
cumulativeLogReturn = sum(rForecast);

% Figures: recent history followed by the k-step forecasts. The colors stay
% visible with both the light and the dark MATLAB figure themes.
observedColor = [0 0.447 0.741];
referenceColor = [0.929 0.694 0.125];
priceHistoryLength = min(250, nPrices);
priceHistory = (nPrices - priceHistoryLength + 1):nPrices;
priceFigure = figure('Visible', 'off');
plot(priceHistory, P(priceHistory), 'Color', observedColor);
hold on;
plot(nPrices:(nPrices + k), priceForecast, 'r.-');
hold off;
title(sprintf('BELEX Price Path Reconstructed from AR(%d) Log-Return Forecasts', p));
xlabel('Observation');
ylabel('Price');
legend('Observed price', 'Reconstructed forecast', 'Location', 'best');
grid on;
saveas(priceFigure, fullfile(figuresDirectory, 'forecast_prices.png'));
close(priceFigure);

returnHistoryLength = min(100, numel(rt));
returnHistory = (numel(rt) - returnHistoryLength + 1):numel(rt);
returnFigure = figure('Visible', 'off');
plot(returnHistory + 1, rt(returnHistory), 'Color', [0.6 0.6 0.6]);
hold on;
plot(Observation, rForecast, 'r.-');
if hasARCH
    plot(Observation, returnLower95, '--', 'Color', referenceColor);
    plot(Observation, returnUpper95, '--', 'Color', referenceColor);
    legend('Observed log return', 'AR conditional mean', ...
        'Gaussian approx. 95% (GARCH variance)', 'Location', 'best');
else
    legend('Observed log return', 'AR conditional mean', 'Location', 'best');
end
hold off;
title(sprintf('AR(%d) Log-Return Forecasts', p));
xlabel('Observation');
ylabel('Log return');
grid on;
saveas(returnFigure, fullfile(figuresDirectory, 'forecast_log_returns.png'));
close(returnFigure);

if hasARCH
    varianceHistoryLength = min(100, numel(s2));
    varianceHistory = (numel(s2) - varianceHistoryLength + 1):numel(s2);
    varianceFigure = figure('Visible', 'off');
    plot(varianceHistory + p + 1, s2(varianceHistory), 'Color', observedColor);
    hold on;
    plot(Observation, sigma2Forecast, 'r.-');
    if isfinite(garchUnconditionalVariance)
        plot([varianceHistory(1) + p + 1, Observation(end)], ...
            garchUnconditionalVariance * [1 1], '--', 'Color', referenceColor);
        legend('Fitted conditional variance', 'GARCH forecast', ...
            'Unconditional variance', 'Location', 'best');
    else
        legend('Fitted conditional variance', 'GARCH forecast', ...
            'Location', 'best');
    end
    hold off;
    title(sprintf('GARCH(%d,%d) Conditional-Variance Forecasts', mG, sG));
    xlabel('Observation');
    ylabel('Conditional variance');
    grid on;
    saveas(varianceFigure, ...
        fullfile(figuresDirectory, 'forecast_conditional_variance.png'));
    close(varianceFigure);
end

% Describe the forecasts.
meanText = sprintf([ ...
    'AR(%d) conditional-mean log-return forecasts: %.6g at h = 1 and ' ...
    '%.6g at h = %d.'], p, rForecast(1), rForecast(k), k);
if isfinite(arLongRunMean)
    meanText = sprintf([ ...
        '%s The AR recursion moves toward the long-run mean ' ...
        'phi_0 / (1 - sum(phi_i)) = %.6g.'], meanText, arLongRunMean);
end
if exist('arAdequate', 'var') && arAdequate ~= 1
    meanText = sprintf([ ...
        '%s AR(%d) failed its adequacy test, so these conditional-mean ' ...
        'forecasts should be read with caution.'], meanText, p);
end

priceText = sprintf([ ...
    'Price path reconstructed with task4 from the last observed price ' ...
    '%.6g: %.6g at h = 1 and %.6g at h = %d (cumulative forecast log return ' ...
    '%.6g, a %.4g%% change). It compounds the conditional-mean log returns ' ...
    'and is a central path, not a price interval.'], ...
    P(end), priceForecast(2), priceForecast(end), k, cumulativeLogReturn, ...
    100 * (priceForecast(end) / P(end) - 1));

if hasARCH
    varianceText = sprintf([ ...
        'GARCH(%d,%d) conditional-variance forecasts: %.6g at h = 1 ' ...
        '(sigma %.6g) and %.6g at h = %d (sigma %.6g).'], mG, sG, ...
        sigma2Forecast(1), sqrt(sigma2Forecast(1)), sigma2Forecast(k), k, ...
        sqrt(sigma2Forecast(k)));
    if isfinite(garchUnconditionalVariance)
        varianceText = sprintf([ ...
            '%s The forecasts move toward the unconditional variance %.6g.'], ...
            varianceText, garchUnconditionalVariance);
    end
    intervalText = sprintf([ ...
        'Gaussian conditional approximation of the 95%% log-return interval, ' ...
        'rForecast +/- 1.96 sqrt(sigma2Forecast): [%.6g, %.6g] at h = 1 and ' ...
        '[%.6g, %.6g] at h = %d.'], returnLower95(1), returnUpper95(1), ...
        returnLower95(k), returnUpper95(k), k);
    intervalText = sprintf([ ...
        '%s Each interval uses only the conditional variance of the return ' ...
        'at that step, so it is a one-step interval evaluated at every ' ...
        'horizon. It does not accumulate the forecast error of the steps in ' ...
        'between, and is therefore narrower than a full h-step forecast ' ...
        'interval for h > 1.'], intervalText);
    if exist('jarqueBeraRejectedRt', 'var') && jarqueBeraRejectedRt
        intervalText = sprintf([ ...
            '%s The Jarque-Bera test rejected Normality of the log returns, ' ...
            'so these intervals are an approximation.'], intervalText);
    end
    if exist('garchAdequate', 'var') && ~garchAdequate
        intervalText = sprintf([ ...
            '%s The GARCH model failed the adequacy test, so the variance ' ...
            'forecasts should be read with caution.'], intervalText);
    end
else
    varianceText = [ ...
        'The ARCH test did not support conditional-variance modeling, so no ' ...
        'GARCH variance forecasts or return intervals were produced.'];
    intervalText = '';
end

reportText = sprintf('%s\n%s\n%s\n', meanText, priceText, varianceText);
if ~isempty(intervalText)
    reportText = sprintf('%s%s\n', reportText, intervalText);
end
fprintf('%s', reportText);

reportFile = fullfile(tablesDirectory, 'forecast_interpretation.txt');
reportId = fopen(reportFile, 'w');
if reportId == -1
    error('BELEX:ReportWriteFailed', ...
        'Could not open the forecast interpretation file for writing: %s', ...
        reportFile);
end
fprintf(reportId, '%s', reportText);
fclose(reportId);

% Keep k, rForecast, sigma2Forecast, priceForecast, the Gaussian return
% bounds, and forecastTable.
clear ConditionalStdForecast ConditionalVarianceForecast GaussianLower95 ...
    GaussianUpper95 Horizon LogReturnForecast Observation PriceForecast ...
    arLagSum figuresDirectory intervalText meanText nPrices observedColor ...
    priceFigure referenceColor ...
    priceHistory priceHistoryLength priceText reportFile reportId ...
    reportText repositoryRoot returnFigure returnHistory ...
    returnHistoryLength scriptPath tablesDirectory varianceFigure ...
    varianceHistory varianceHistoryLength varianceText;
