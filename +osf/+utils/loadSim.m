function sim = loadSim(filename)

    if nargin < 1 || isempty(filename)
        error('loadSim:NoFilename', 'Filename must be specified.');
    end
    if ~endsWith(filename, '.mat', 'IgnoreCase', true)
        filename = [filename '.mat'];
    end
    if exist(filename, 'file') ~= 2
        error('loadSim:FileNotFound', 'File not found: %s', filename);
    end
    S = load(filename, 'sim');
    if ~isfield(S, 'sim')
        error('loadSim:MissingSimVar', 'No variable ''sim'' found in file: %s', filename);
    end
    sim = S.sim;
    if ~isa(sim, 'osf.Sim')
        error('loadSim:InvalidType', 'Loaded object is not an instance of osf.Sim');
    end

end
