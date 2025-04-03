function mask = genMask(arraySize, type, regionSize, position)
    if nargin < 4
        position = [];
    end

    % Determine dimensionality
    if arraySize(1) == 1 || arraySize(2) == 1
        dim = 1;
        arrayLength = max(arraySize);
    elseif ismatrix(arraySize)
        dim = 2;
    else
        error('Unsupported arraySize format.');
    end

    if isempty(position)
        position = zeros(1, dim);
    end

    if dim == 1
        x = 1:arrayLength;
        cx = round(arrayLength / 2 + position(1));

        switch type
            case 'rectangle'
                hw = round(regionSize(1) / 2);
                mask = false(1, arrayLength);
                mask(max(1, cx - hw):min(arrayLength, cx + hw)) = true;
            otherwise
                error('Unsupported shape type for 1D: %s', type);
        end

    elseif dim == 2
        [h, w] = deal(arraySize(1), arraySize(2));
        [X, Y] = meshgrid(1:w, 1:h);
        cx = round(w / 2 + position(1));
        cy = round(h / 2 + position(2));

        switch type
            case 'rectangle'
                hw = round(regionSize(1) / 2);
                hh = round(regionSize(2) / 2);
                mask = false(h, w);
                mask(max(1, cy - hh):min(h, cy + hh), max(1, cx - hw):min(w, cx + hw)) = true;

            case 'circle'
                radius = regionSize(1);
                mask = (X - cx).^2 + (Y - cy).^2 <= radius^2;

            case 'annulus'
                r1 = regionSize(1); % inner
                r2 = regionSize(2); % outer
                d2 = (X - cx).^2 + (Y - cy).^2;
                mask = (d2 <= r2^2) & (d2 > r1^2);

            otherwise
                error('Unsupported shape type for 2D: %s', type);
        end
    end
end
