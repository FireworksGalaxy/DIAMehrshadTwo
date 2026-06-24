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
