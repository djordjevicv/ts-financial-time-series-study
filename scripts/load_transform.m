% Load and transform BELEX data.

scriptPath = mfilename('fullpath');
repositoryRoot = fileparts(fileparts(scriptPath));
dataFile = fullfile(repositoryRoot, 'data', 'belex.mat');

if exist(dataFile, 'file') ~= 2
    error('BELEX:MissingFile', ...
        'BELEX data file was not found: %s', dataFile);
end

try
    loadedData = load(dataFile);
catch loadError
    error('BELEX:LoadFailed', ...
        'Could not load BELEX data file "%s": %s', ...
        dataFile, loadError.message);
end

if ~isfield(loadedData, 'belex')
    error('BELEX:MissingVariable', ...
        'BELEX data file "%s" does not contain the required variable "belex".', ...
        dataFile);
end

belex = loadedData.belex;

if ~isnumeric(belex) || ~isreal(belex) || ~isvector(belex)
    error('BELEX:InvalidSeries', ...
        'Variable "belex" must be a real numeric vector.');
end

if numel(belex) < 2
    error('BELEX:TooFewObservations', ...
        'Variable "belex" must contain at least two price observations.');
end

if any(~isfinite(belex(:)))
    error('BELEX:NonfinitePrice', ...
        'Variable "belex" contains NaN or Inf values; every price must be finite.');
end

if any(belex(:) <= 0)
    error('BELEX:NonpositivePrice', ...
        'Variable "belex" contains zero or negative values; every price must be greater than zero.');
end

P = belex(:);
fprintf('Loaded %d BELEX price observations.\n', numel(P));

returnsFigure = figure('Visible', 'off');
[Rt, Gt, rt] = returns(P);
close(returnsFigure);

expectedReturnCount = numel(P) - 1;
if any([numel(Rt), numel(Gt), numel(rt)] ~= expectedReturnCount)
    error('BELEX:InvalidReturnLength', ...
        'Each return series must contain exactly %d observations.', ...
        expectedReturnCount);
end

if any(~isfinite(Rt)) || any(~isfinite(Gt)) || any(~isfinite(rt))
    error('BELEX:NonfiniteReturn', ...
        'Calculated return series contain NaN or Inf values.');
end

figuresDirectory = fullfile(repositoryRoot, 'results', 'figures');
if exist(figuresDirectory, 'dir') ~= 7
    mkdir(figuresDirectory);
end

priceFigure = figure;
plot(1:numel(P), P);
title('BELEX Prices');
xlabel('Observation');
ylabel('Price');
grid on;
saveas(priceFigure, fullfile(figuresDirectory, 'belex_prices.png'));

logReturnFigure = figure;
plot(2:numel(P), rt);
title('BELEX Log Returns');
xlabel('Observation');
ylabel('Log return');
grid on;
saveas(logReturnFigure, ...
    fullfile(figuresDirectory, 'belex_log_returns.png'));

clear belex dataFile expectedReturnCount figuresDirectory loadedData ...
    logReturnFigure priceFigure repositoryRoot returnsFigure scriptPath;
