function data = prepData(obj, plotType, varargin)
    p = inputParser;
    addRequired(p, 'obj');
    addRequired(p, 'plotType');
    addParameter(p, 'unwrap', false, @islogical);
    addParameter(p, 'wdfLimit', 1, @(x) isnumeric(x) && (x <= 1 && x > 0));
    addParameter(p, 'overlayDetector', false, @(x) isa(x, 'osf.Detector'));
    addParameter(p, 'title', '', @ischar);
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
        if ~isempty(title)
            data.title = title;
        else
            data.title = 'Field Intensity';
        end
        data.xlabel = 'x (mm)';
        data.ylabel = 'y (mm)';
        data.colorbarLabel = 'Intensity';

    case {'phase', 'p', 'phase.unwrap', 'p.unwrap'}
        title = p.Results.title;
        unwrapFlag = p.Results.unwrap;
        data.xAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, cols) * 1e3;
        data.yAxis = linspace(-obj.fieldLength/2, obj.fieldLength/2, rows) * 1e3;
        if unwrapFlag || ismember(plotType, {'phase.unwrap', 'p.unwrap'})
            imageData = osf.utils.phaseUnwrap(obj.phase);
        else
            imageData = obj.phase;
        end
        data.imageData = imageData;
        data.cmap = obj.cmap;
        if ~isempty(title)
            data.title = title;
        else
            data.title = 'Field Phase';
        end
        data.xlabel = 'x (mm)';
        data.ylabel = 'y (mm)';
        data.colorbarLabel = 'Phase';

    case {'cross', 'a.cross', 'amplitude.cross'}
        xpos = round(rows/2);
        ypos = round(cols/2);
        [data.xAxis, data.xCross] = obj.getAmplitudeCross('x', xpos);
        [data.yAxis, data.yCross] = obj.getAmplitudeCross('y', ypos);
        data.xAxis = data.xAxis * 1e3;
        data.yAxis = data.yAxis * 1e3;
        data.xlabel = 'Position (mm)';
        data.ylabel = 'Amplitude';
        data.title = sprintf('Amplitude Cross Section (x = %.2fmm, y = %.2fmm)', ...
        xpos * obj.resolution * 1e3, ypos * obj.resolution * 1e3);

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
        xpos * obj.resolution * 1e3, ypos * obj.resolution * 1e3);

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
        data.title = title;

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
        data.title = title;
        data.xlabel = 'Column';
        data.ylabel = 'Row';
        data.colorbarLabel = 'Value';

    otherwise
        error('Unhandled plotType: %s', plotType);
    end

    function data = prepSim(obj, data)
        if isempty(obj.elements)
            data.isempty = true;
            return
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
        data.axisLimits = [min(0, -.5*currentX), currentX + .5*currentX, -0.1, 0.1];
        data.elementHeight = elementHeight;
        data.labelOffset = labelOffset;
        data.componentCenters = componentCenters;
        data.baselineY = -0.03;
        data.tickHeight = 0.002;

        parax = osf.parax.ParaxialSystem(obj);
        parax = osf.parax.ParaxialSystem.solveSystem(parax);
        data.paraxial.dist = cumsum(parax.distances);
        data.paraxial.marginalRay = parax.marginalRay.heights * ...
        (0.8 * elementHeight / max(2 * abs(parax.marginalRay.heights)));
        data.paraxial.chiefRay = parax.chiefRay.heights * ...
        (0.8 * elementHeight / max(2 * abs(parax.chiefRay.heights)));
    end
end
