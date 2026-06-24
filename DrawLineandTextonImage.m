function DrawLineandTextonImage(BW,I2,config,B,L)
% DrawLineandTextonImage - Diagnostic overlay: draw each detected object's
% boundary on the image and label it with its index at the centroid.
%
% Usage:
%   DrawLineandTextonImage(BW, I2, config, B, L)
%
% Input:
%   BW     - binary mask (used only if B/L are not supplied).
%   I2     - grayscale image shown alongside the mask.
%   config - configuration struct (drawing is skipped if show_figures is false).
%   B, L   - (Optional) precomputed boundaries/label image to avoid recomputing.
%
% This function is for visualization only; it has no return value.

% Reuse precomputed boundaries/labels when supplied; otherwise compute them.
if nargin < 5 || isempty(B) || isempty(L)
    [B, L] = bwboundaries(BW,'noholes');
end

if config.show_figures
    figure('Name','AfterCounting');
    imshowpair(I2,BW,'montage');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'red', 'LineWidth', 4)
    end
    title('Overlayed images')
    hold off

    % Label each object with its index at the centroid
    s = regionprops(L, 'Centroid');
    hold on
    for k = 1:numel(s)
        c = s(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
    hold off
end
end
