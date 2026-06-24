function FinalExtration=FindRealOverlappedObjectsAmongExtracted(Extracted_object,DeformedObject,areaThreshold,config)
% FindRealOverlappedObjectsAmongExtracted - Keep only the deformed candidates
% that are large enough to actually be overlapped bubbles.
%
% Usage:
%   FinalExtration = FindRealOverlappedObjectsAmongExtracted(...
%                        Extracted_object, DeformedObject, areaThreshold, config)
%
% A deformed object is treated as a genuine overlap only if its area exceeds
% areaThreshold (config.area_threshold_overlapped). Smaller deformed objects
% are assumed to be single (merely non-round) bubbles and are left out.
%
% Output:
%   FinalExtration - binary mask of the confirmed overlapped objects.

L=Extracted_object;
stats = regionprops(L, 'Area', 'Eccentricity', 'Solidity'); % Get properties

% Identify objects with multiple bubbles (area-based heuristic)
areaValues = [stats.Area];
s = regionprops(logical(L), 'Centroid');
selectedLabels = find(areaValues > areaThreshold);

% Create a binary mask with only selected objects
selectedObjects = ismember(L, selectedLabels);

% Display results
if config.show_figures
    figure('Name', 'The First Image after the image selection');
    subplot(1,2,1); imshow(L); title('Original Image');
    subplot(1,2,2); imshow(selectedObjects); title('Extracted Multi-Bubble Objects');
    hold on;
    for k = 1:size(DeformedObject,1)
        boundary = DeformedObject{k};
        plot(boundary(:,2), boundary(:,1), 'blue', 'LineWidth', 4)
        if k <= numel(s)
            c = s(k).Centroid;
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle','Color','red');
        end
    end
    hold off;
end

FinalExtration = selectedObjects;
end
