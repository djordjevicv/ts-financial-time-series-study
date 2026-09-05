% Compare serial dependence in BELEX prices and log returns.

if ~exist('P', 'var') || ~exist('rt', 'var')
    error('BELEX:MissingSeries', ...
        'Run scripts/load_transform.m before this script.');
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

% Use the same maximum lag for a direct price-versus-return comparison.
m = 20;

% acfgraf computes the price ACF. tratio2 redraws it and circles lags that
% are statistically significant at its 5% level.
priceFigure = figure('Visible', 'off');
acfPrices = acfgraf(P, m);
testPrices = tratio2(P, m);
title('BELEX Price ACF');
xlabel('Lag');
ylabel('Sample ACF');
xlim([1 m]);
grid on;
saveas(priceFigure, fullfile(figuresDirectory, 'price_acf.png'));
close(priceFigure);

% Repeat the identical diagnostic procedure for log returns.
returnFigure = figure('Visible', 'off');
acfReturns = acfgraf(rt, m);
testReturns = tratio2(rt, m);
title('BELEX Log-Return ACF');
xlabel('Lag');
ylabel('Sample ACF');
xlim([1 m]);
grid on;
saveas(returnFigure, fullfile(figuresDirectory, 'log_return_acf.png'));
close(returnFigure);

% tratio2 rows are: lag, ACF, standard error, t-ratio, and significance.
% Combine both series into one table so every reported lag is reproducible.
Series = [repmat({'Price'}, m, 1); repmat({'Log return'}, m, 1)];
Lag = [testPrices(1, :)'; testReturns(1, :)'];
ACF = [acfPrices(:); acfReturns(:)];
StandardError = [testPrices(3, :)'; testReturns(3, :)'];
TRatio = [testPrices(4, :)'; testReturns(4, :)'];
Significant = logical([testPrices(5, :)'; testReturns(5, :)']);
acfSignificance = table(Series, Lag, ACF, StandardError, TRatio, Significant);
writetable(acfSignificance, ...
    fullfile(tablesDirectory, 'acf_significance.csv'));

% Row 5 equals one for significant lags and zero otherwise.
significantPriceLags = find(testPrices(5, :) ~= 0);
significantReturnLags = find(testReturns(5, :) ~= 0);
nSignificantPrice = numel(significantPriceLags);
nSignificantReturn = numel(significantReturnLags);

% Describe the observed change instead of assuming that returns must have
% weaker serial dependence.
if nSignificantReturn < nSignificantPrice
    transformationSummary = sprintf([ ...
        'Converting prices to log returns reduces the significant-lag count ' ...
        'by %d, from %d to %d.'], ...
        nSignificantPrice - nSignificantReturn, ...
        nSignificantPrice, nSignificantReturn);
elseif nSignificantReturn == nSignificantPrice
    transformationSummary = sprintf([ ...
        'Converting prices to log returns leaves the significant-lag count ' ...
        'unchanged at %d.'], nSignificantReturn);
else
    transformationSummary = sprintf([ ...
        'Converting prices to log returns increases the significant-lag count ' ...
        'by %d, from %d to %d.'], ...
        nSignificantReturn - nSignificantPrice, ...
        nSignificantPrice, nSignificantReturn);
end

% Significant return ACF lags are initial AR candidates. PACFAR selects the
% final AR order in the conditional-mean modeling stage.
if nSignificantReturn == 0
    arOrderGuidance = [ ...
        'No return lags from 1 to 20 are significant. Use AR(0) as a benchmark ' ...
        'and confirm the final order with PACFAR.'];
else
    arOrderGuidance = sprintf([ ...
        'For AR-order selection, examine significant return lags %s and ' ...
        'confirm the final order with PACFAR.'], ...
        mat2str(significantReturnLags));
end

% Print the conclusions and save the same text with the tabular results.
reportText = sprintf([ ...
    'Price: %d of %d lags are significant; positions %s.\n' ...
    'Log return: %d of %d lags are significant; positions %s.\n' ...
    '%s Thus, %d of %d tested return lags (%.1f%%) retain significant ' ...
    'serial dependence.\n%s\n'], ...
    nSignificantPrice, m, mat2str(significantPriceLags), ...
    nSignificantReturn, m, mat2str(significantReturnLags), ...
    transformationSummary, nSignificantReturn, m, ...
    100 * nSignificantReturn / m, arOrderGuidance);

fprintf('%s', reportText);
reportFile = fullfile(tablesDirectory, 'acf_interpretation.txt');
reportId = fopen(reportFile, 'w');
if reportId == -1
    error('BELEX:ReportWriteFailed', ...
        'Could not open the ACF interpretation file for writing: %s', reportFile);
end
fprintf(reportId, '%s', reportText);
fclose(reportId);

% Keep the diagnostic outputs in the workspace for later stages, while
% removing variables used only to construct files and messages.
clear ACF Lag Series Significant StandardError TRatio arOrderGuidance ...
    figuresDirectory priceFigure reportFile reportId reportText ...
    repositoryRoot returnFigure scriptPath tablesDirectory ...
    transformationSummary;
