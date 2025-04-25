function saveSim(sim, filename)

    if ~isa(sim, 'osf.Sim')
        error('saveSim:InvalidType', 'sim must be an instance of osf.Sim');
    end

    if nargin < 2 || isempty(filename)
        ts = datestr(now, 'yymmdd_HHMMSS');
        filename = [ts '_sim'];
    else
        filename = [filename '.mat'];
    end

    save(filename, 'sim');

end
