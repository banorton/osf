function X_norm = normalize(X)
    % NORMALIZE - Scales an array or vector so that its values are between 0 and 1
    %
    % Usage:
    %   X_norm = normalize(X)
    %
    % Input:
    %   X - Input vector or matrix (numeric)
    %
    % Output:
    %   X_norm - Normalized output, with min(X_norm) = 0 and max(X_norm) = 1
    %
    % Description:
    %   - The function rescales the values of X so that the minimum value becomes 0
    %     and the maximum value becomes 1.
    %   - If all values in X are the same, it returns an array of zeros.

    % Ensure X is numeric
    if ~isnumeric(X)
        error('Input must be a numeric array or vector.');
    end

    % Find minimum and maximum values
    X_min = min(X(:));
    X_max = max(X(:));

    % Prevent division by zero when all values are the same
    if X_max == X_min
        X_norm = zeros(size(X));
    else
        X_norm = (X - X_min) / (X_max - X_min);
    end
end
