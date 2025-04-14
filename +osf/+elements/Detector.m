classdef Detector < osf.elements.Element
    properties
        name
        elementType
        apertureType
        apertureParams
        dim
        id

        pixelPitch
        resolution
        bitDepth = 8
        show
    end

    methods
        function obj = Detector(resolution, pixelPitch, varargin)
            p = inputParser;
            addRequired(p, 'pixelPitch', @isnumeric);
            addRequired(p, 'resolution', @isnumeric);
            addParameter(p, 'name', 'Detector', @ischar);
            addParameter(p, 'dim', 2, @(x) isnumeric(x) && ismember(x, [1, 2]));
            addParameter(p, 'bitDepth', 8, @isnumeric);
            addParameter(p, 'show', true, @islogical);
            addParameter(p, 'id', 0, @isnumeric);
            parse(p, pixelPitch, resolution, varargin{:});

            obj.elementType = 'detector';
            obj.apertureType = 'none';
            obj.apertureParams = struct();
            obj.name = p.Results.name;
            obj.dim = p.Results.dim;
            obj.pixelPitch = p.Results.pixelPitch;
            obj.resolution = p.Results.resolution;
            obj.bitDepth = p.Results.bitDepth;
            obj.show = p.Results.show;
            obj.id = p.Results.id;
        end

        function field = apply(obj, field)
            % Field physical size
            [fieldRows, fieldCols] = size(field.intensity);
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
            projected = interp2(XField, YField, field.intensity, XDet, YDet, 'linear', 0);

            % Quantize
            levels = 2^obj.bitDepth;
            normalized = projected / max(projected(:));
            quantized = round(normalized * (levels - 1)) / (levels - 1);

            if obj.show
                osf.show(quantized, title="Detector");
            end
        end

        function [fmin, fmax] = bandwidth(obj)
            pitch = obj.pixelPitch;
            res = obj.resolution;

            fmax.x = 1 / (2 * pitch);
            fmax.y = 1 / (2 * pitch);

            fmin.x = 1 / (res(1) * pitch);
            fmin.y = 1 / (res(2) * pitch);
        end

        function field = phaseFunction(obj, varargin)
        end
    end
end
