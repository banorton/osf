classdef (Abstract) Element < handle
    % Abstract class for optical elements

    properties (Abstract)
        name           % Name of the optical element
        phaseFunction  % Function handle for phase modification (if applicable)
        elementType    % Type of optical element (e.g., 'plane', 'lens', 'custom')
        apertureType   % Type of aperture ('none', 'circ', 'rect', etc.)
        apertureParams % Parameters for defining the aperture (radius, length, width, etc.)
        dim            % Dimensionality of the element (1 or 2)
    end

    methods (Abstract)
        apply(field);
    end

    methods

        function distanceMatrix = getDistanceMatrix(obj, field)
            % Get field dimensions and resolution
            [Ny, Nx] = size(field.amplitude);
            dx = field.resolution;

            % Create meshgrid for calculating distances from center
            if obj.dim == 2
                % Create a grid of physical coordinates centered at zero
                [X, Y] = meshgrid(linspace(-Nx/2, Nx/2, Nx) * dx, ...
                                  linspace(-Ny/2, Ny/2, Ny) * dx);
                % Calculate the distance from the center for each point
                distanceMatrix = struct('X', X, 'Y', Y);

            elseif obj.dim == 1
                % Create a vector of physical coordinates centered at zero
                X = linspace(-Nx/2, Nx/2, Nx) * dx;
                distanceMatrix = struct('X', X);

            else
                error('Dimensionality must be either 1 or 2.');
            end
        end

        function field = applyAperture(obj, field)
            % Apply the aperture of the optical element to the input field
            fieldLength = field.fieldLength;
            if obj.dim == 2
                fieldSize = size(field.amplitude);

                if strcmp(obj.apertureType, 'circ')
                    radius = obj.apertureParams.radius;
                    [X, Y] = meshgrid(linspace(-fieldLength/2, fieldLength/2, fieldSize(2)), ...
                                      linspace(-fieldLength/2, fieldLength/2, fieldSize(1)));
                    apertureMask = sqrt(X.^2 + Y.^2) <= radius;
                elseif strcmp(obj.apertureType, 'rect')
                    w = obj.apertureParams.width;
                    l = obj.apertureParams.length;
                    [X, Y] = meshgrid(linspace(-fieldLength/2, fieldLength/2, fieldSize(2)), ...
                                      linspace(-fieldLength/2, fieldLength/2, fieldSize(1)));
                    apertureMask = (abs(X) <= w/2) & (abs(Y) <= l/2);
                else
                    error('Unsupported aperture type for 2D case.');
                end

            elseif obj.dim == 1
                fieldSize = length(field.amplitude);

                if strcmp(obj.apertureType, 'circ')
                    radius = obj.apertureParams.radius;
                    x = linspace(-fieldLength/2, fieldLength/2, fieldSize);
                    apertureMask = abs(x) <= radius;
                elseif strcmp(obj.apertureType, 'rect')
                    l = obj.apertureParams.length;
                    x = linspace(-fieldLength/2, fieldLength/2, fieldSize);
                    apertureMask = abs(x) <= l/2;
                else
                    error('Unsupported aperture type for 1D case.');
                end

            else
                error('Dimensionality must be either 1 or 2.');
            end

            field.amplitude = field.amplitude .* apertureMask;
            field.phase = field.phase .* apertureMask;
        end

        function obj = addCircAperture(obj, radius)
            % Add a circular aperture to the element
            obj.apertureType = 'circ';
            obj.apertureParams.radius = radius;
        end

        function obj = addRectAperture(obj, width, length)
            % Add a rectangular aperture to the element
            obj.apertureType = 'rect';
            obj.apertureParams.width = width;
            obj.apertureParams.length = length;
        end

        function wrappedPhi = wrap(~, phi)
            % Wrap phase to the range [-pi, pi]
            wrappedPhi = atan2(sin(phi), cos(phi));
        end

    end
end
