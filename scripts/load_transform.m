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

clear belex loadedData dataFile repositoryRoot scriptPath;
