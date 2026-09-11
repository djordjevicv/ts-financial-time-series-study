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

reportText = sprintf('%s\n%s\n', archText, squaredAcfText);
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
clear archAlpha archCriticalValue archLagCount archQ archText ...
    figuresDirectory reportFile reportId reportText repositoryRoot ...
    scriptPath squaredAcfFigure squaredAcfLags squaredAcfText ...
    squaredResidualTRatio tablesDirectory;
