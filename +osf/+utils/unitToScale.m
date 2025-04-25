function scale = unitToScale(unitStr)
    % Convert a length unit (e.g. 'um','µm','micron','millimeter', etc.) to its numeric scale factor.

    u = lower(strtrim(string(unitStr)));
    switch u
    case {'nm','nanometer','nanometre'}
        scale = 1e-9;
    case {'um','µm','μm','micron','micrometer','micrometre'}
        scale = 1e-6;
    case {'mm','millimeter','millimetre'}
        scale = 1e-3;
    case {'cm','centimeter','centimetre'}
        scale = 1e-2;
    case {'m','meter','metre'}
        scale = 1;
    case {'km','kilometer','kilometre'}
        scale = 1e3;
    otherwise
        error('unitToScale:UnknownUnit', 'Unknown unit: %s', unitStr);
    end
end
