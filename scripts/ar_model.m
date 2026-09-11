% Fit the AR conditional-mean model to BELEX log returns.

if ~exist('rt', 'var')
    error('BELEX:MissingSeries', ...
        'Run scripts/load_transform.m before this script.');
end

if ~isnumeric(rt) || ~isreal(rt) || size(rt, 2) ~= 1 || any(~isfinite(rt))
    error('BELEX:InvalidSeries', ...
        'Log returns rt must be a finite, real numeric column vector.');
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

nReturns = numel(rt);

% PACFAR fits AR(1), ..., AR(ceil(log(T))) by least squares and keeps the
% largest order whose highest-lag coefficient is significant at the 5% level.
% Repeat the same candidate tests with ARestimate so the choice is auditable.
maxCandidateOrder = ceil(log(nReturns));
CandidateOrder = (1:maxCandidateOrder)';
HighestLagCoefficient = zeros(maxCandidateOrder, 1);
HighestLagStandardError = zeros(maxCandidateOrder, 1);
HighestLagTRatio = zeros(maxCandidateOrder, 1);
HighestLagSignificant = false(maxCandidateOrder, 1);
for candidateOrder = 1:maxCandidateOrder
    candidateEstimate = ARestimate(rt, candidateOrder);
    HighestLagCoefficient(candidateOrder) = candidateEstimate(1, 1);
    HighestLagStandardError(candidateOrder) = candidateEstimate(2, 1);
    HighestLagTRatio(candidateOrder) = ...
        candidateEstimate(1, 1) / candidateEstimate(2, 1);
    HighestLagSignificant(candidateOrder) = candidateEstimate(3, 1) == 1;
end
arOrderSelection = table(CandidateOrder, HighestLagCoefficient, ...
    HighestLagStandardError, HighestLagTRatio, HighestLagSignificant);
writetable(arOrderSelection, ...
    fullfile(tablesDirectory, 'ar_order_selection.csv'));

p = PACFAR(rt);

significantCandidateOrders = find(HighestLagSignificant)';
if isempty(significantCandidateOrders)
    expectedOrder = 0;
else
    expectedOrder = max(significantCandidateOrders);
end
if p ~= expectedOrder
    error('BELEX:OrderSelectionMismatch', ...
        'PACFAR selected AR(%d), but the candidate table implies AR(%d).', ...
        p, expectedOrder);
end

% Estimate the selected model. The supplied ordering is
% xAR = [phi_p; ...; phi_1; phi_0].
xAR = ARLScoef(rt, p);
arInfo = ARestimate(rt, p);
[arAdequate, a] = adequateAR(rt, xAR);
a = a(:);

if numel(xAR) ~= p + 1 || size(arInfo, 2) ~= p + 1
    error('BELEX:InvalidARCoefficients', ...
        'AR(%d) must have exactly %d coefficients.', p, p + 1);
end
if numel(a) ~= nReturns - p || any(~isfinite(a))
    error('BELEX:InvalidARResiduals', ...
        'AR(%d) must produce %d finite residuals.', p, nReturns - p);
end
if ~isscalar(arAdequate) || ~any(arAdequate == [0 1])
    error('BELEX:InvalidAdequacyResult', ...
        'adequateAR must return a scalar 0 or 1.');
end

% arInfo rows are: estimate, standard error, significance indicator, and
% the estimate with insignificant coefficients set to zero.
Lag = [(p:-1:1)'; 0];
Term = cell(p + 1, 1);
for coefficientIndex = 1:(p + 1)
    Term{coefficientIndex} = sprintf('phi_%d', Lag(coefficientIndex));
end
Estimate = xAR(:);
StandardError = arInfo(2, :)';
TRatio = Estimate ./ StandardError;
Significant = arInfo(3, :)' == 1;
RestrictedEstimate = arInfo(4, :)';
arCoefficients = table(Term, Lag, Estimate, StandardError, TRatio, ...
    Significant, RestrictedEstimate);
writetable(arCoefficients, fullfile(tablesDirectory, 'ar_coefficients.csv'));

% adequateAR returns only its decision. Recompute the Ljung-Box statistic that
% it passes to ljungbox3 (10 lags; degrees of freedom reduced by the number
% of nonzero coefficients) so the report can show Q next to the decision.
arLjungBoxLags = 10;
arFittedParameters = sum(abs(sign(xAR)));
arLjungBoxDf = arLjungBoxLags - arFittedParameters;
chiSquare95 = [3.84 5.99 7.81 9.48 11.07 12.59 14.07 15.51 16.92 18.31];
arLjungBoxCriticalValue = chiSquare95(arLjungBoxDf);
nResiduals = numel(a);
ljungBoxSum = 0;
for residualLag = 1:arLjungBoxLags
    ljungBoxSum = ljungBoxSum + ...
        acf(a, residualLag)^2 / (nResiduals - residualLag);
end
arLjungBoxQ = nResiduals * (nResiduals + 2) * ljungBoxSum;
if (arLjungBoxQ <= arLjungBoxCriticalValue) ~= (arAdequate == 1)
    error('BELEX:AdequacyMismatch', ...
        'The recomputed Ljung-Box decision disagrees with adequateAR.');
end

arAdequacy = table(p, arLjungBoxLags, arFittedParameters, arLjungBoxDf, ...
    arLjungBoxQ, arLjungBoxCriticalValue, arAdequate == 1, ...
    'VariableNames', {'AROrder', 'LjungBoxLags', 'FittedParameters', ...
    'DegreesOfFreedom', 'Q', 'CriticalValue', 'Adequate'});
writetable(arAdequacy, fullfile(tablesDirectory, 'ar_adequacy.csv'));

% The residual a(t) belongs to the price observation p + 1 + t, because rt(1)
% is the return from observation 1 to observation 2.
residualObservations = (p + 2):(nReturns + 1);
residualFigure = figure('Visible', 'off');
plot(residualObservations, a);
title(sprintf('AR(%d) Residuals of BELEX Log Returns', p));
xlabel('Observation');
ylabel('Residual');
grid on;
saveas(residualFigure, fullfile(figuresDirectory, 'ar_residuals.png'));
close(residualFigure);

residualAcfLags = 20;
residualAcfFigure = figure('Visible', 'off');
acfgraf(a, residualAcfLags);
residualTRatio = tratio2(a, residualAcfLags);
title(sprintf('ACF of AR(%d) Residuals', p));
xlabel('Lag');
ylabel('Sample ACF');
xlim([1 residualAcfLags]);
grid on;
saveas(residualAcfFigure, fullfile(figuresDirectory, 'ar_residual_acf.png'));
close(residualAcfFigure);
significantResidualLags = find(residualTRatio(5, :) ~= 0);

% Describe the order selection.
if isempty(significantCandidateOrders)
    orderText = sprintf([ ...
        'PACFAR tested AR(1) to AR(%d). No highest-lag coefficient was ' ...
        'significant, so the selected model is AR(0), a constant mean.'], ...
        maxCandidateOrder);
else
    orderText = sprintf([ ...
        'PACFAR tested AR(1) to AR(%d). The highest-lag coefficient was ' ...
        'significant for orders %s, so the largest of them, AR(%d), was ' ...
        'selected.'], maxCandidateOrder, ...
        mat2str(significantCandidateOrders), p);
end

if exist('significantReturnLags', 'var') && ...
        any(significantReturnLags > maxCandidateOrder)
    orderText = sprintf([ ...
        '%s Significant log-return ACF lags %s from the price-versus-return ' ...
        'stage lie beyond the largest candidate order and are not modeled ' ...
        'by this AR specification.'], orderText, ...
        mat2str(significantReturnLags(significantReturnLags > maxCandidateOrder)));
end

% Describe the coefficients, listing significant lags separately from the
% constant.
significantLagTerms = Lag(Significant & Lag > 0)';
insignificantLagTerms = Lag(~Significant & Lag > 0)';
coefficientLines = '';
for coefficientIndex = 1:(p + 1)
    if Significant(coefficientIndex)
        significanceLabel = 'significant';
    else
        significanceLabel = 'not significant';
    end
    coefficientLines = sprintf('%s  %s = %.6g (SE %.6g, t = %.4g, %s)\n', ...
        coefficientLines, Term{coefficientIndex}, Estimate(coefficientIndex), ...
        StandardError(coefficientIndex), TRatio(coefficientIndex), ...
        significanceLabel);
end
if p == 0
    coefficientSummary = 'The model contains only the constant phi_0.';
elseif isempty(insignificantLagTerms)
    coefficientSummary = sprintf( ...
        'All %d lag coefficients are significant at the 5%% level.', p);
else
    coefficientSummary = sprintf([ ...
        'Lag coefficients %s are significant at the 5%% level; lags %s are ' ...
        'not. PACFAR tests only the highest lag, so the full AR(%d) is kept ' ...
        'and the insignificant coefficients are reported rather than removed.'], ...
        mat2str(significantLagTerms), mat2str(insignificantLagTerms), p);
end

if arAdequate == 1
    adequacyText = sprintf([ ...
        'Ljung-Box test on the residuals: Q = %.6g with %d lags and ' ...
        '%d degree(s) of freedom, not above the critical value %.6g. The ' ...
        'residuals show no significant remaining autocorrelation, so AR(%d) ' ...
        'is adequate for the conditional mean.'], arLjungBoxQ, arLjungBoxLags, ...
        arLjungBoxDf, arLjungBoxCriticalValue, p);
else
    adequacyText = sprintf([ ...
        'Ljung-Box test on the residuals: Q = %.6g with %d lags and ' ...
        '%d degree(s) of freedom, above the critical value %.6g. Significant ' ...
        'autocorrelation remains, so AR(%d) is not adequate for the ' ...
        'conditional mean.'], arLjungBoxQ, arLjungBoxLags, arLjungBoxDf, ...
        arLjungBoxCriticalValue, p);
end

if isempty(significantResidualLags)
    residualAcfText = sprintf( ...
        'No residual ACF lag from 1 to %d is individually significant.', ...
        residualAcfLags);
else
    residualAcfText = sprintf( ...
        'Individually significant residual ACF lags (1 to %d): %s.', ...
        residualAcfLags, mat2str(significantResidualLags));
end

reportText = sprintf('%s\nSelected AR(%d) coefficients:\n%s%s\n%s\n%s\n', ...
    orderText, p, coefficientLines, coefficientSummary, adequacyText, ...
    residualAcfText);
fprintf('%s', reportText);

reportFile = fullfile(tablesDirectory, 'ar_model_interpretation.txt');
reportId = fopen(reportFile, 'w');
if reportId == -1
    error('BELEX:ReportWriteFailed', ...
        'Could not open the AR interpretation file for writing: %s', ...
        reportFile);
end
fprintf(reportId, '%s', reportText);
fclose(reportId);

% Keep the documented Stage 4 outputs (p, xAR, arInfo, arAdequate, a) and
% the result tables; remove variables used only to build files and messages.
clear CandidateOrder Estimate HighestLagCoefficient HighestLagSignificant ...
    HighestLagStandardError HighestLagTRatio Lag RestrictedEstimate ...
    Significant StandardError TRatio Term adequacyText candidateEstimate ...
    candidateOrder chiSquare95 coefficientIndex coefficientLines ...
    coefficientSummary expectedOrder figuresDirectory ...
    insignificantLagTerms ljungBoxSum nResiduals nReturns orderText ...
    reportFile reportId reportText repositoryRoot residualAcfFigure ...
    residualAcfLags residualAcfText residualFigure residualLag ...
    residualObservations residualTRatio scriptPath significanceLabel ...
    significantLagTerms tablesDirectory;
