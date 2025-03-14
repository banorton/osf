function [vals, freq, pos] = wvd(signal, varargin)
    p = inputParser;
    addParameter(p, 'sampleRate', 1000, @isnumeric);
    addParameter(p, 'numFrequencyPoints', 1000, @isinteger);

    parse(p, varargin{:});
    sampleRate = p.Results.sampleRate;
    numFrequencyPoints = p.Results.numFrequencyPoints;

    [vals, freq, pos] = wvd(signal, sampleRate, 'smoothedPseudo', MinThreshold=0, NumFrequencyPoints=numFrequencyPoints);
end
