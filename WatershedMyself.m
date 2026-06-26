function watershedLabels=WatershedMyself(new_image,config)
% WatershedMyself - Split touching/overlapped bubbles using a watershed transform.
%
% Usage:
%   watershedLabels = WatershedMyself(new_image, config)
%
% Input:
%   new_image - RGB image containing only the overlapped clusters (background
%               flattened by MaskOverlappedBubblesBeforeWatershed).
%   config    - configuration struct from BubbleDetectionConfig().
%
% Output:
%   watershedLabels - label image where 0 marks watershed ridges (separation
%                     lines) and each positive value is one segmented region.
%
% Method: sharpen -> adaptive binarize -> distance transform -> extended-minima
% markers -> watershed. The extended-minima depth controls how aggressively
% clusters are split (smaller = more splits).
sharped=FastUnsharpMask(new_image,config.sharpen_radius_watershed,config.sharpen_amount_watershed); % Fast (FFT) unsharp mask; radius/amount from config
img=rgb2gray(sharped);
bw = imbinarize(img, 'adaptive', 'foreground', 'dark', 'Sensitivity', config.binarize_sensitivity_watershed);

% Complement the binary image (if needed)
bw_complement = imcomplement(bw);

% figure('Name','Changing the colors');
% imshowpair(bw,bw_complement,'montage');

% Perform distance transform
distanceTransform = bwdist(~bw_complement,'euclidean');

% figure('Name','Distance Transform important');
% imshowpair(bw,distanceTransform,'montage');

% Apply watershed transformation.
% Optionally adapt the extended-minima depth to the distance-transform scale
% so markers track cluster size instead of using a single fixed value.
if isfield(config,'watershed_adaptive_markers') && config.watershed_adaptive_markers
    hmin = config.extended_minima_fraction * max(distanceTransform(:));
else
    hmin = config.extended_minima_value;
end
mask = imextendedmin(distanceTransform, hmin); % Marker sensitivity
modifiedDistance = imimposemin(-distanceTransform, mask);
watershedLabels = watershed(modifiedDistance);

% Remove the background label
segmentedBubbles = bw;
segmentedBubbles(watershedLabels == 0) = 0;

% Visualize results
if config.show_figures
    figure;
    subplot(1, 2, 1); imshow(bw); title('Binary Image');
    subplot(1, 2, 2); imshow(label2rgb(watershedLabels, 'jet', 'k', 'shuffle')); title('Watershed Segmentation');
end


end

% =========================================================================
% HOW THE WATERSHED SPLIT WORKS - a short story with one example
% =========================================================================
%
% Imagine two bubbles that touch and overlap, like this binary blob:
%
%        ____        ____
%       /    \      /    \
%      | bub  |    | bub  |
%      |  A   |----|  B   |       (the "----" is the neck where they touch)
%       \____/      \____/
%
% We want one label for A and one label for B, with a dividing line at the
% neck. The four key lines above do exactly that. Here is the journey.
%
% -------------------------------------------------------------------------
% STEP 1:  distanceTransform = bwdist(~bw_complement,'euclidean');
% -------------------------------------------------------------------------
% For every pixel, measure the distance to the nearest edge of the shape.
% Pixels deep inside a bubble are FAR from any edge (large value); pixels
% near the rim are CLOSE to an edge (small value). A profile drawn through
% both bubble centres looks like two hills with a dip (the neck) between:
%
%   value
%    3 |   *           *          <- the two bubble "centres" (deepest inside)
%    2 |  * *    v    * *         <- v = neck dip between the bubbles
%    1 | *   *  * *  *   *
%    0 |*_____*___*_____*____
%        bubble A  neck  bubble B
%
% Note: for a perfectly round bubble the top is a single point. For an
% elongated/deformed bubble the top is a whole RIDGE (the medial axis /
% skeleton), not a single centroid point.
%
% -------------------------------------------------------------------------
% STEP 2:  mask = imextendedmin(distanceTransform, hmin);
% -------------------------------------------------------------------------
% This finds the marker "seeds" - one blob of TRUE pixels per bubble centre.
% It returns a binary image (the LOCATIONS), not a number; counting the
% separate TRUE blobs tells you HOW MANY seeds were found.
%
%   Ideally: 1 TRUE blob per bubble  ->  here that means 2 seeds (A and B).
%
% hmin is the depth filter (how deep a basin must be to count as a seed):
%   - hmin too small -> extra seeds inside one bubble  -> OVER-split.
%   - hmin too large -> A and B merge into one seed     -> UNDER-split.
% So tuning hmin really means "get exactly one marker per bubble".
% (hmin = extended_minima_value, or a fraction of the max distance if
%  watershed_adaptive_markers is enabled.)
%
% -------------------------------------------------------------------------
% STEP 3:  modifiedDistance = imimposemin(-distanceTransform, mask);
% -------------------------------------------------------------------------
% Watershed always floods UPWARD starting from the lowest points (minima).
% Two things happen here:
%   (a) The negative sign (-distanceTransform) flips the terrain so the
%       bubble centres (were the highest hills) become the deepest VALLEYS,
%       and the rims/neck become high RIDGES.
%   (b) imimposemin does NOT search for the valleys - we hand it the seed
%       locations via `mask`. It forces the terrain to have valleys ONLY at
%       those seed spots and digitally fills in every other dip, so random
%       noise cannot create false extra valleys.
%
%   after this step: exactly 2 valleys (A and B), one neck ridge, nothing else.
%
% -------------------------------------------------------------------------
% STEP 4:  watershedLabels = watershed(modifiedDistance);
% -------------------------------------------------------------------------
% Picture rain filling the two valleys at the same time. Each valley grows
% its own pool. Where the two rising pools meet - at the neck ridge - a dam
% (the watershed line) is built. That dam is the split between the bubbles.
%
%        A A A | B B B          | = watershed ridge line (label 0)
%        A A A | B B B          A = region label 1
%        A A A | B B B          B = region label 2
%
% Output `watershedLabels`:
%     0  -> the dividing line(s) between regions,
%     1  -> all pixels of bubble A,
%     2  -> all pixels of bubble B,
%   (the surrounding background also gets its own label, which the rest of
%    the pipeline removes as an over-large region later).
%
% -------------------------------------------------------------------------
% ONE-LINE SUMMARY
% -------------------------------------------------------------------------
%   bwdist        -> find how deep each pixel sits inside the shape,
%   imextendedmin -> mark WHERE the bubble centres are (the seeds),
%   imimposemin   -> make ONLY those seeds the valleys to flood from,
%   watershed     -> flood from the seeds and draw a dam at the neck.
% =========================================================================
