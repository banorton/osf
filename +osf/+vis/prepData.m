function data = prepData(obj, plotType, varargin)
    p = inputParser;
    addRequired(p, 'obj');
    addRequired(p, 'plotType');
    addParameter(p, 'unwrap', true, @islogical);
    addParameter(p, 'wdfLimit', 1, @(x) isnumeric(x) && (x <= 1 && x > 0));
    addParameter(p, 'overlayDetector', false, @(x) isa(x, 'osf.Detector'));
    addParameter(p, 'title', 'default', @ischar);
    addParameter(p, 'zoom', 1, @isnumeric);
    parse(p, obj, plotType, varargin{:});

    plotType = p.Results.plotType;

    if strcmp(plotType, 'default')
        plotType = sprintf('%s.default', class(obj));
    end

    data.meta.objectType = class(obj);
    data.meta.plotType = plotType;

    if isa(obj, 'osf.Field')
        [rows, cols] = size(obj.amplitude);
    elseif isnumeric(obj) && isreal(obj)
        [rows, cols] = size(obj);
    end

    switch plotType
    case {'osf.Field.default', 'amplitude', 'a'}
        title = p.Results.title;
        data.xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols) * 1e3;
        data.yAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows) * 1e3;
        data.imageData = obj.amplitude;
        data.cmap = 'gray';
        if strcmp(title, 'default')
            data.title = title;
        else
            data.title = 'Field Amplitude';
        end
        data.xlabel = 'x (mm)';
        data.ylabel = 'y (mm)';
        data.colorbarLabel = 'Amplitude';

    case {'phase', 'p'}
        title = p.Results.title;
        unwrapFlag = p.Results.unwrap;
        data.xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols) * 1e3;
        data.yAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows) * 1e3;
        if unwrapFlag
            imageData = osf.utils.phaseUnwrap(obj.phase);
        else
            imageData = obj.phase;
        end
        data.imageData = imageData;
        data.cmap = obj.cmap;
        if strcmp(title, 'default')
            data.title = title;
        else
            data.title = 'Field Phase';
        end
        data.xlabel = 'x (mm)';
        data.ylabel = 'y (mm)';
        data.colorbarLabel = 'Phase';

    case {'cross', 'a.cross', 'amplitude.cross'}
        title = p.Results.title;
        xpos = round(rows/2);
        ypos = round(cols/2);
        [data.xAxis, data.xCross] = obj.getAmplitudeCross('x', xpos);
        [data.yAxis, data.yCross] = obj.getAmplitudeCross('y', ypos);
        data.xAxis = data.xAxis * 1e3;
        data.yAxis = data.yAxis * 1e3;
        data.xlabel = 'Position (mm)';
        data.ylabel = 'Amplitude';
        if strcmp(title, 'default')
            data.title = sprintf('Amplitude Cross Section (x = %.2fmm, y = %.2fmm)', ...
            (floor(rows/2) - xpos) * obj.resolution * 1e3, (floor(cols/2) - ypos) * obj.resolution * 1e3);
        else
            data.title = title;
        end
        data.zoom = p.Results.zoom;

    case {'phase.cross', 'p.cross'}
        xpos = round(rows/2);
        ypos = round(cols/2);
        [data.xAxis, data.xCross] = obj.getPhaseCross('x', xpos);
        [data.yAxis, data.yCross] = obj.getPhaseCross('y', ypos);
        data.xAxis = data.xAxis * 1e3;
        data.yAxis = data.yAxis * 1e3;
        data.xCross = unwrap(data.xCross);
        data.yCross = unwrap(data.yCross);
        data.xlabel = 'Position (mm)';
        data.ylabel = 'Phase';
        data.title = sprintf('Phase Cross Section (x = %.2fmm, y = %.2fmm)', ...
        (floor(rows/2) - xpos) * obj.resolution * 1e3, (floor(cols/2) - ypos) * obj.resolution * 1e3);
        data.zoom = p.Results.zoom;

    case {'wdf', 'amplitude.wdf', 'a.wdf', 'wdf.x', 'amplitude.wdf.x', 'a.wdf.x'}
        title = p.Results.title;
        wdfLimit = p.Results.wdfLimit;
        if isa(p.Results.overlayDetector, 'osf.Detector')
            d = p.Results.overlayDetector;

            width = d.resolution(1) * d.pixelPitch;
            x0 = -width/2;
            x1 = width/2;
            data.xBox = [x0, x1, x1, x0];

            [fmin, fmax] = d.bandwidth();
            data.yBox = [fmin.x, fmin.x, fmax.x, fmax.x];
        end
        [data.xAxis, cross] = obj.getAmplitudeCross();
        [imageData, yAxis, ~] = osf.utils.wvd(cross, sampleRate=1/obj.resolution);
        [imgRows, ~] = size(imageData);
        data.imageData = imageData(1:round(imgRows*wdfLimit),:);
        data.yAxis = yAxis(1:round(imgRows*wdfLimit));
        data.colormap = 'default';
        data.xlabel = 'x (m)';
        data.ylabel = 'Local Spatial Frequency (m^{-1})';
        data.colorbarLabel = 'Wigner-Ville Intensity';
        if strcmp(title, 'default')
            data.title = '';
        end

    case {'phase.wdf', 'p.wdf'}
        wdfLimit = p.Results.wdfLimit;
        [data.xAxis, cross] = obj.getPhaseCross();
        [imageData, yAxis, ~] = osf.utils.wvd(cross, sampleRate=obj.fieldLength/obj.resolution);
        [imgRows, ~] = size(imageData);
        data.imageData = imageData(1:round(imgRows*wdfLimit),:);
        data.yAxis = yAxis(1:round(imgRows*wdfLimit));
        data.colormap = 'default';
        data.xlabel = 'x (m)';
        data.ylabel = 'Local Spatial Frequency (m^{-1})';
        data.colorbarLabel = 'Wigner-Ville Intensity';

    case {'osf.Sim.default'}
        data = prepSim(obj, data);

    case {'double.default', 'single.default', 'uint8.default', 'uint16.default'}
        title = p.Results.title;
        data.imageData = obj;
        data.xAxis = 1:size(obj, 2);
        data.yAxis = 1:size(obj, 1);
        data.cmap = 'gray';
        data.xlabel = '';
        data.ylabel = '';
        data.colorbarLabel = 'Value';
        if strcmp(title, 'default')
            data.title = '';
        end

    otherwise
        error('Unhandled plotType: %s', plotType);
    end

    function data = prepSim(obj, data)
        % Prepares simulation data for plotting and ray tracing.
        %
        % Required
        %   data : Structure to populate with simulation information.
        %
        % This function collects information about elements (positions, names,
        % types, etc.) and computes global settings. If there is more than one
        % element, it performs the paraxial ray tracing; otherwise, it skips it.

        if isempty(obj.elements)
            data.isempty = true;
            return;
        else
            data.isempty = false;
        end

        elementHeight = 0.04;
        labelOffset = 0.02;
        currentX = 0;
        componentCenters = [];

        for i = 1:length(obj.elements)
            currentX = currentX + obj.distances(i);
            element = obj.elements{i};
            data.elements(i).position = currentX;
            data.elements(i).name = element.name;
            data.elements(i).type = lower(element.elementType);
            componentCenters = [componentCenters, currentX];
        end

        data.distances = obj.distances;
        data.axisLimits = [min(0, -0.5*currentX), currentX + 0.5*currentX, -0.1, 0.1];
        data.elementHeight = elementHeight;
        data.labelOffset = labelOffset;
        data.componentCenters = componentCenters;
        data.baselineY = -0.03;
        data.tickHeight = 0.002;

        % Only perform paraxial ray tracing if there is more than one element.
        if isscalar(obj.elements)
            data.paraxial = [];
            data.axisLimits = [-.02, .02, -0.1, 0.1];
        else
            parax = osf.parax.ParaxialSystem(obj);
            parax = osf.parax.ParaxialSystem.solveSystem(parax);
            data.paraxial.dist = cumsum(parax.distances);
            data.paraxial.marginalRay = parax.marginalRay.heights * (0.8 * elementHeight / max(2 * abs(parax.marginalRay.heights)));
            data.paraxial.chiefRay = parax.chiefRay.heights * (0.8 * elementHeight / max(2 * abs(parax.chiefRay.heights)));
        end
    end
end
