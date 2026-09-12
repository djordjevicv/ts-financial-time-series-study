% Test the AR residuals for ARCH effects.

if ~exist('a', 'var') || ~exist('p', 'var')
    error('BELEX:MissingResiduals', ...
        'Run scripts/ar_model.m before this script.');
end

if ~isnumeric(a) || ~isreal(a) || size(a, 2) ~= 1 || any(~isfinite(a))
    error('BELEX:InvalidResiduals', ...
        'AR residuals a must be a finite, real numeric column vector.');
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

% ARCHeffect applies the Ljung-Box test to a.^2 with ceil(log(T)) lags and
% returns [lagCount; alpha; Q; criticalValue], not a decision.
archResult = ARCHeffect(a);
if numel(archResult) ~= 4 || any(~isfinite(archResult))
    error('BELEX:InvalidARCHResult', ...
        'ARCHeffect must return four finite values.');
end
hasARCH = archResult(3) > archResult(4);

archLagCount = archResult(1);
archAlpha = archResult(2);
archQ = archResult(3);
archCriticalValue = archResult(4);
archTest = table(p, archLagCount, archAlpha, archQ, archCriticalValue, ...
    hasARCH, 'VariableNames', {'AROrder', 'LagCount', 'Alpha', 'Q', ...
    'CriticalValue', 'HasARCH'});
writetable(archTest, fullfile(tablesDirectory, 'arch_effect.csv'));

% The squared residuals satisfy a2(t) ~ s2(t) = alpha_0 + sum(alpha_i
% a2(t-i)), an AR-like model, so the PACF of a.^2 indicates the ARCH order
% the data support. This is evidence about the order, not a fitted model.
squaredResiduals = a.^2;
archPacfOrder = PACFAR(squaredResiduals);
maxArchCandidate = ceil(log(numel(squaredResiduals)));
ArchCandidateOrder = (1:maxArchCandidate)';
HighestLagCoefficient = zeros(maxArchCandidate, 1);
HighestLagStandardError = zeros(maxArchCandidate, 1);
HighestLagTRatio = zeros(maxArchCandidate, 1);
HighestLagSignificant = false(maxArchCandidate, 1);
for candidateOrder = 1:maxArchCandidate
    candidateEstimate = ARestimate(squaredResiduals, candidateOrder);
    HighestLagCoefficient(candidateOrder) = candidateEstimate(1, 1);
    HighestLagStandardError(candidateOrder) = candidateEstimate(2, 1);
    HighestLagTRatio(candidateOrder) = ...
        candidateEstimate(1, 1) / candidateEstimate(2, 1);
    HighestLagSignificant(candidateOrder) = candidateEstimate(3, 1) == 1;
end
Selected = ArchCandidateOrder == archPacfOrder;
archOrderPacf = table(ArchCandidateOrder, HighestLagCoefficient, ...
    HighestLagStandardError, HighestLagTRatio, HighestLagSignificant, ...
    Selected, 'VariableNames', {'CandidateOrder', ...
    'HighestLagCoefficient', 'HighestLagStandardError', 'HighestLagTRatio', ...
    'HighestLagSignificant', 'Selected'});
writetable(archOrderPacf, ...
    fullfile(tablesDirectory, 'arch_order_pacf.csv'));

% Plot the ACF of the squared residuals. Significant lags here are the
% autocorrelation that the ARCH test summarizes.
squaredAcfLags = 20;
squaredAcfFigure = figure('Visible', 'off');
acfgraf(a.^2, squaredAcfLags);
squaredResidualTRatio = tratio2(a.^2, squaredAcfLags);
title(sprintf('ACF of Squared AR(%d) Residuals', p));
xlabel('Lag');
ylabel('Sample ACF');
xlim([1 squaredAcfLags]);
grid on;
saveas(squaredAcfFigure, ...
    fullfile(figuresDirectory, 'ar_squared_residual_acf.png'));
close(squaredAcfFigure);
significantSquaredLags = find(squaredResidualTRatio(5, :) ~= 0);

if hasARCH
    archText = sprintf([ ...
        'ARCH test on the squared AR(%d) residuals: Q = %.6g with %d lags, ' ...
        'above the %.6g critical value at alpha = %.2g. The squared residuals ' ...
        'are autocorrelated, so the ARCH effect is statistically significant ' ...
        'and the conditional variance is modeled with GARCH.'], ...
        p, archQ, archLagCount, archCriticalValue, archAlpha);
else
    archText = sprintf([ ...
        'ARCH test on the squared AR(%d) residuals: Q = %.6g with %d lags, ' ...
        'not above the %.6g critical value at alpha = %.2g. The ARCH effect ' ...
        'is not statistically significant, so a GARCH model is not ' ...
        'supported and the pipeline skips conditional-variance modeling.'], ...
        p, archQ, archLagCount, archCriticalValue, archAlpha);
end

if isempty(significantSquaredLags)
    squaredAcfText = sprintf( ...
        'No squared-residual ACF lag from 1 to %d is individually significant.', ...
        squaredAcfLags);
else
    squaredAcfText = sprintf( ...
        '%d of %d squared-residual ACF lags are individually significant: %s.', ...
        numel(significantSquaredLags), squaredAcfLags, ...
        mat2str(significantSquaredLags));
end

if archPacfOrder == 0
    archOrderText = sprintf([ ...
        'The PACF of the squared residuals is not significant at any order ' ...
        'up to %d, so it points to no ARCH order.'], maxArchCandidate);
elseif archPacfOrder >= 3
    archOrderText = sprintf([ ...
        'The PACF of the squared residuals points to ARCH(%d) (orders up to ' ...
        '%d tested). A pure ARCH model of that order needs %d parameters, ' ...
        'so GARCH(1,1) is used instead: it represents the same persistence ' ...
        'with three.'], archPacfOrder, maxArchCandidate, archPacfOrder + 1);
else
    archOrderText = sprintf([ ...
        'The PACF of the squared residuals points to ARCH(%d) (orders up to ' ...
        '%d tested). GARCH(1,1) is fitted as the project choice; a pure ' ...
        'ARCH(%d) is a compact alternative worth comparing.'], ...
        archPacfOrder, maxArchCandidate, archPacfOrder);
end

reportText = sprintf('%s\n%s\n%s\n', archText, squaredAcfText, ...
    archOrderText);
fprintf('%s', reportText);

reportFile = fullfile(tablesDirectory, 'arch_effect_interpretation.txt');
reportId = fopen(reportFile, 'w');
if reportId == -1
    error('BELEX:ReportWriteFailed', ...
        'Could not open the ARCH interpretation file for writing: %s', ...
        reportFile);
end
fprintf(reportId, '%s', reportText);
fclose(reportId);

% Keep archResult, hasARCH, and the result table for later stages.
clear ArchCandidateOrder HighestLagCoefficient HighestLagSignificant ...
    HighestLagStandardError HighestLagTRatio Selected archAlpha ...
    archCriticalValue archLagCount archOrderText archQ archText ...
    candidateEstimate candidateOrder figuresDirectory reportFile reportId ...
    reportText repositoryRoot scriptPath squaredAcfFigure squaredAcfLags ...
    squaredAcfText squaredResidualTRatio squaredResiduals ...
    tablesDirectory;
