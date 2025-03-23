classdef Detector < osf.Element
    properties
        name
        elementType
        apertureType
        apertureParams
        dim

        pixelPitch
        chipSize
        resolution
        bitDepth = 8
    end

    methods
        function obj = Detector(resolution, pixelPitch, chipSize, varargin)
            p = inputParser;
            addParameter(p, 'name', '', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'bitDepth', 8, @isnumeric);
            addRequired(p, 'pixelPitch', @isnumeric);
            addRequired(p, 'chipSize', @isnumeric);
            addRequired(p, 'resolution', @isnumeric);
            parse(p, pixelPitch, chipSize, resolution, varargin{:});

            obj.elementType = 'detector';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            obj.pixelPitch = p.Results.pixelPitch;
            obj.chipSize = p.Results.chipSize;
            obj.resolution = p.Results.resolution;
            obj.bitDepth = p.Results.bitDepth;
        end

        function field = apply(obj, field)
            intensity = field.amplitude .^ 2;

            % Field physical size
            [fieldRows, fieldCols] = size(intensity);
            fieldSizeX = fieldCols * field.resolution;
            fieldSizeY = fieldRows * field.resolution;

            % Detector resolution
            detCols = obj.resolution(1);
            detRows = obj.resolution(2);
            detSizeX = detCols * obj.pixelPitch;
            detSizeY = detRows * obj.pixelPitch;

            % Coordinates
            xField = linspace(-fieldSizeX/2, fieldSizeX/2, fieldCols);
            yField = linspace(-fieldSizeY/2, fieldSizeY/2, fieldRows);
            xDet = linspace(-detSizeX/2, detSizeX/2, detCols);
            yDet = linspace(-detSizeY/2, detSizeY/2, detRows);

            % Interpolate
            [XField, YField] = meshgrid(xField, yField);
            [XDet, YDet] = meshgrid(xDet, yDet);
            projected = interp2(XField, YField, intensity, XDet, YDet, 'linear', 0);

            % Quantize
            levels = 2^obj.bitDepth;
            normalized = projected / max(projected(:));
            quantized = round(normalized * (levels - 1)) / (levels - 1);

            osf.show(quantized, title="Detector");
        end

        function field = phaseFunction(obj, varargin)
        end
    end
end
