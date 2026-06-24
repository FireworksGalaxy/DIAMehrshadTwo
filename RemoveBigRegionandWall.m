function LL=RemoveBigRegionandWall(LL,area_tolerance,BBLength)
% Remove large regions (likely walls / merged blobs) whose area exceeds
% area_tolerance. Label re-ordering is done by the caller (ReorderingLabels).
stats = regionprops(LL, 'Area');

% Consider only the first BBLength labels (matches original scope)
nLabels = min(BBLength, numel(stats));
areas = [stats(1:nLabels).Area];

% Labels to remove (no index gaps, safe when none qualify)
labelsToRemove = find(areas > area_tolerance);

for lbl = labelsToRemove
    LL(LL == lbl) = 0; % Set all pixels of this label to background (0)
end
end
