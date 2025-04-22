function out = minMaxNorm(img, range)
    if nargin < 2
        range = [0 1];
    end

    img = double(img);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    out = range(1) + img * (range(2) - range(1));
end
