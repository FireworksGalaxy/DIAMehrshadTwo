function ShowLabels(LL,figTitle,config)
% ShowLabels - Display labeled image with region numbers at centroids
% Usage: ShowLabels(LL, figTitle) or ShowLabels(LL, figTitle, config)

% Use default config if not provided (for backwards compatibility)
if nargin < 3
    config = BubbleDetectionConfig();
end

% Only display if show_figures is enabled
if ~config.show_figures
    return;
end

img_test=LL;
% Display the labeled image
figure('Name',figTitle);
imshow(img_test, []);
colormap jet; % Optional: You can apply a colormap to make regions easier to distinguish

% Hold the current figure to overlay text
hold on;

% Get properties of the labeled regions (centroids, etc.)
stats_label = regionprops(img_test, 'Centroid');

% Loop through each region and add the label number at its centroid
for k = 1:length(stats_label)
    % Get the centroid of the region
    centroid = stats_label(k).Centroid;
    
    % Use the label number as text
    labelNumber = k;  % You can change this if the labels aren't 1, 2, 3, ...
    
    % Place the label number at the centroid
    text(centroid(1), centroid(2), num2str(labelNumber), ...
        'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Release the hold
hold off;
end
