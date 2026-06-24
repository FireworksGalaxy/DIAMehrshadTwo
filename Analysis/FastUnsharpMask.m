function out = FastUnsharpMask(A, radius, amount)
% FastUnsharpMask - Unsharp masking equivalent to imsharpen, but computing the
% Gaussian blur in the frequency domain (FFT). This is dramatically faster than
% imsharpen for large radii, where imsharpen convolves a huge kernel spatially.
%
%   out = FastUnsharpMask(A, radius, amount)
%
% Mirrors imsharpen(A,'Radius',radius,'Amount',amount):
%   out = A + amount * (A - gaussianBlur(A, sigma = radius))
% computed in double precision with replicate padding and a saturating cast
% back to the input class (matching imsharpen's behaviour with Threshold = 0).
%
% NOTE: This reproduces imsharpen to within a small border-truncation
% difference (imsharpen's internal Gaussian uses a different, undocumented
% truncation), so the detected object count can differ by ~1% versus imsharpen.
%
% Works for grayscale and multi-channel (RGB) images.

    classA = class(A);
    Ad = double(A);

    blurred = zeros(size(Ad));
    for c = 1:size(Ad, 3)
        blurred(:, :, c) = imgaussfilt(Ad(:, :, c), radius, ...
            'FilterDomain', 'frequency', 'Padding', 'replicate');
    end

    sharpened = Ad + amount * (Ad - blurred);
    out = cast(sharpened, classA);   % saturating cast, like imsharpen
end
