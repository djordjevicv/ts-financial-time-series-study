% Model the conditional variance of the AR residuals with GARCH.

if ~exist('a', 'var') || ~exist('hasARCH', 'var')
    error('BELEX:MissingARCHTest', ...
        'Run scripts/ar_model.m and scripts/arch_effect.m before this script.');
end

if ~hasARCH
    error('BELEX:NoARCHEffect', ...
        ['The ARCH test on the AR residuals is not significant. GARCH is ' ...
        'fitted only when hasARCH is true.']);
end

if exist('fmincon', 'file') == 0
    error('BELEX:MissingOptimizationToolbox', ...
        'GARCHcoef requires fmincon from the Optimization Toolbox.');
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

% Project choice: start with GARCH(1,1). The supplied parameter ordering is
% xGARCH = [alpha_m; ...; alpha_1; alpha_0; beta_s; ...; beta_1].
mG = 1;
sG = 1;
nGARCHParameters = mG + 1 + sG;

xGARCHDirect = GARCHcoef(a, mG, sG);

% Convergence check. GARCHcoef starts fmincon at alpha_0 = 0.1, which is
% several orders of magnitude above the variance of daily log-return
% residuals, so the optimization on a is badly scaled and can stop early.
% The GARCH likelihood is scale-equivariant: replacing a by c*a multiplies
% alpha_0 by c^2, leaves every alpha_i and beta_j unchanged, and shifts
% fmlGARCH by a constant. Refit on percent residuals (c = 100), map alpha_0
% back, and keep the estimate with the lower negative log-likelihood on a.
garchScale = 100;
xGARCHScaled = GARCHcoef(garchScale * a, mG, sG);
xGARCHRescaled = xGARCHScaled;
xGARCHRescaled(mG + 1) = xGARCHScaled(mG + 1) / garchScale^2;

negLogLikelihoodDirect = fmlGARCH(a, xGARCHDirect, mG, sG);
negLogLikelihoodRescaled = fmlGARCH(a, xGARCHRescaled, mG, sG);

% The rescaled estimate must also satisfy GARCHcoef's constraints on a.
garchConstraintTolerance = 1e-8;
garchLowerBounds = [zeros(mG, 1); 1e-6; zeros(sG, 1)];
rescaledFeasible = all(xGARCHRescaled >= ...
    garchLowerBounds - garchConstraintTolerance) && ...
    sum(xGARCHRescaled([1:mG, mG + 2:end])) <= ...
    1 - 1e-6 + garchConstraintTolerance;

likelihoodTolerance = 1e-6 * max(1, abs(negLogLikelihoodDirect));
if rescaledFeasible && isfinite(negLogLikelihoodRescaled) && ...
        (~isfinite(negLogLikelihoodDirect) || ...
        negLogLikelihoodRescaled < negLogLikelihoodDirect - likelihoodTolerance)
    xGARCH = xGARCHRescaled;
    garchEstimateSource = 'Percent-scaled refit';
else
    xGARCH = xGARCHDirect;
    garchEstimateSource = 'Direct GARCHcoef';
end
xGARCH = xGARCH(:);

s2 = s2sequenceGARCH(a, xGARCH, mG, sG);
s2 = s2(:);
if numel(s2) ~= numel(a) || any(~isfinite(s2)) || any(s2 <= 0)
    error('BELEX:InvalidConditionalVariance', ...
        'The fitted conditional variance must be finite and positive.');
end

% adequateGARCH applies ljungbox2 to a./sqrt(s2) and returns
% [lagCount; alpha; Q; criticalValue], not a decision.
garchResult = adequateGARCH(a, xGARCH, mG, sG);
garchAdequate = garchResult(3) <= garchResult(4);

% Supplementary check with the supplied ARCH test: are the squared
% standardized residuals still autocorrelated after the GARCH fit?
standardizedResiduals = a ./ sqrt(s2);
standardizedArchResult = ARCHeffect(standardizedResiduals);
remainingARCH = standardizedArchResult(3) > standardizedArchResult(4);

% Parameter table in the supplied ordering.
Term = cell(nGARCHParameters, 1);
for alphaIndex = 1:mG
    Term{alphaIndex} = sprintf('alpha_%d', mG - alphaIndex + 1);
end
Term{mG + 1} = 'alpha_0';
for betaIndex = 1:sG
    Term{mG + 1 + betaIndex} = sprintf('beta_%d', sG - betaIndex + 1);
end
Estimate = xGARCH;
DirectEstimate = xGARCHDirect(:);
RescaledEstimate = xGARCHRescaled(:);
garchCoefficients = table(Term, Estimate, DirectEstimate, RescaledEstimate);
writetable(garchCoefficients, ...
    fullfile(tablesDirectory, 'garch_coefficients.csv'));

% Persistence sum(alpha_i) + sum(beta_j) governs how slowly volatility shocks
% decay; below one it implies a finite unconditional variance.
garchPersistence = sum(xGARCH([1:mG, mG + 2:end]));
if garchPersistence < 1
    garchUnconditionalVariance = xGARCH(mG + 1) / (1 - garchPersistence);
else
    garchUnconditionalVariance = Inf;
end
if garchPersistence > 0 && garchPersistence < 1
    garchHalfLife = log(0.5) / log(garchPersistence);
else
    garchHalfLife = Inf;
end
residualSampleVariance = var(a);

garchSummary = table(mG, sG, {garchEstimateSource}, ...
    negLogLikelihoodDirect, negLogLikelihoodRescaled, garchPersistence, ...
    garchUnconditionalVariance, residualSampleVariance, garchHalfLife, ...
    garchResult(1), garchResult(2), garchResult(3), garchResult(4), ...
    garchAdequate, standardizedArchResult(3), standardizedArchResult(4), ...
    remainingARCH, ...
    'VariableNames', {'m', 's', 'EstimateSource', ...
    'NegLogLikelihoodDirect', 'NegLogLikelihoodRescaled', 'Persistence', ...
    'UnconditionalVariance', 'ResidualSampleVariance', 'HalfLife', ...
    'AdequacyLagCount', 'AdequacyAlpha', 'AdequacyQ', ...
    'AdequacyCriticalValue', 'Adequate', 'StandardizedSquaredQ', ...
    'StandardizedSquaredCriticalValue', 'RemainingARCH'});
writetable(garchSummary, fullfile(tablesDirectory, 'garch_summary.csv'));

% Plot the residuals with the Gaussian 95% band implied by the fitted
% conditional standard deviation, and the conditional variance itself.
residualObservations = (p + 2):(p + 1 + numel(a));
volatilityFigure = figure('Visible', 'off');
subplot(2, 1, 1);
plot(residualObservations, a, 'Color', [0.6 0.6 0.6]);
hold on;
plot(residualObservations, 1.96 * sqrt(s2), 'r');
plot(residualObservations, -1.96 * sqrt(s2), 'r');
hold off;
title(sprintf('AR(%d) Residuals and GARCH(%d,%d) +/-1.96 sigma', p, mG, sG));
xlabel('Observation');
ylabel('Residual');
legend('Residual a', '+/-1.96 sqrt(s2)', 'Location', 'best');
grid on;
subplot(2, 1, 2);
plot(residualObservations, s2);
title(sprintf('GARCH(%d,%d) Conditional Variance', mG, sG));
xlabel('Observation');
ylabel('Conditional variance');
grid on;
saveas(volatilityFigure, ...
    fullfile(figuresDirectory, 'garch_conditional_variance.png'));
close(volatilityFigure);

% Describe the observed results.
parameterLines = '';
for parameterIndex = 1:nGARCHParameters
    parameterLines = sprintf('%s  %s = %.6g\n', parameterLines, ...
        Term{parameterIndex}, Estimate(parameterIndex));
end

if strcmp(garchEstimateSource, 'Direct GARCHcoef')
    estimationText = sprintf([ ...
        'GARCHcoef on a reached negative log-likelihood %.8g; the ' ...
        'percent-scaled refit reached %.8g on the same scale. The direct ' ...
        'estimate is kept because the refit did not improve the likelihood.'], ...
        negLogLikelihoodDirect, negLogLikelihoodRescaled);
else
    estimationText = sprintf([ ...
        'GARCHcoef on a stopped at negative log-likelihood %.8g, while the ' ...
        'percent-scaled refit reached %.8g on the same scale. The refit is ' ...
        'kept because it has the higher likelihood; the direct estimate is ' ...
        'recorded in garch_coefficients.csv.'], ...
        negLogLikelihoodDirect, negLogLikelihoodRescaled);
end

if garchPersistence < 1
    persistenceText = sprintf([ ...
        'Persistence sum(alpha_i) + sum(beta_j) = %.6g, so volatility ' ...
        'shocks decay with a half-life of about %.3g observations. The ' ...
        'implied unconditional variance is %.6g, compared with the sample ' ...
        'variance %.6g of the AR residuals.'], garchPersistence, ...
        garchHalfLife, garchUnconditionalVariance, residualSampleVariance);
else
    persistenceText = sprintf([ ...
        'Persistence sum(alpha_i) + sum(beta_j) = %.6g, so the model does ' ...
        'not imply a finite unconditional variance.'], garchPersistence);
end

if garchAdequate
    adequacyText = sprintf([ ...
        'adequateGARCH: Ljung-Box Q = %.6g on the standardized residuals ' ...
        'a./sqrt(s2) with %d lags, not above the %.6g critical value. The ' ...
        'GARCH(%d,%d) model is adequate by the course criterion.'], ...
        garchResult(3), garchResult(1), garchResult(4), mG, sG);
else
    adequacyText = sprintf([ ...
        'adequateGARCH: Ljung-Box Q = %.6g on the standardized residuals ' ...
        'a./sqrt(s2) with %d lags, above the %.6g critical value. The ' ...
        'GARCH(%d,%d) model is not adequate by the course criterion; ' ...
        'its variance forecasts should be read with caution.'], ...
        garchResult(3), garchResult(1), garchResult(4), mG, sG);
end

if remainingARCH
    remainingText = sprintf([ ...
        'Supplementary ARCHeffect on the standardized residuals: Q = %.6g, ' ...
        'above %.6g, so some ARCH effect remains after the fit.'], ...
        standardizedArchResult(3), standardizedArchResult(4));
else
    remainingText = sprintf([ ...
        'Supplementary ARCHeffect on the standardized residuals: Q = %.6g, ' ...
        'not above %.6g (compared with Q = %.6g before the GARCH fit), so ' ...
        'no significant ARCH effect remains.'], ...
        standardizedArchResult(3), standardizedArchResult(4), archResult(3));
end

reportText = sprintf('GARCH(%d,%d) estimates (%s):\n%s%s\n%s\n%s\n%s\n', ...
    mG, sG, garchEstimateSource, parameterLines, estimationText, ...
    persistenceText, adequacyText, remainingText);
fprintf('%s', reportText);

reportFile = fullfile(tablesDirectory, 'garch_model_interpretation.txt');
reportId = fopen(reportFile, 'w');
if reportId == -1
    error('BELEX:ReportWriteFailed', ...
        'Could not open the GARCH interpretation file for writing: %s', ...
        reportFile);
end
fprintf(reportId, '%s', reportText);
fclose(reportId);

% Keep mG, sG, xGARCH, s2, garchResult, garchAdequate, and the result tables.
clear DirectEstimate Estimate RescaledEstimate Term adequacyText ...
    alphaIndex betaIndex estimationText figuresDirectory ...
    garchConstraintTolerance garchLowerBounds garchScale ...
    likelihoodTolerance parameterIndex parameterLines persistenceText ...
    remainingText reportFile reportId reportText repositoryRoot ...
    rescaledFeasible residualObservations scriptPath tablesDirectory ...
    volatilityFigure xGARCHScaled;
