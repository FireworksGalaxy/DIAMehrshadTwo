function [FinalSegmentedObjects, separationInfo] = SeparateOverlappedObjectsByStrategy(BW, I_Formask, RealOverlap, L, overlapPoints, config)
% SeparateOverlappedObjectsByStrategy - Separate overlapped bubble clusters.
%
% Two-bubble overlaps are separated with the normal watershed workflow.
% Multi-overlap clusters can also use watershed, or be skipped while tuning
% the two-bubble/red-region behavior first.

FinalSegmentedObjects = zeros(size(BW));
separationInfo = struct('Label', {}, 'DefectCount', {}, 'Method', {}, ...
    'CurveCoefficients', {}, 'CurveCenter', {}, 'CurveScale', {}, 'UsedFallback', {});

if ~any(RealOverlap(:))
    return;
end

overlapLabels = unique(L(RealOverlap));
overlapLabels(overlapLabels == 0) = [];
if isempty(overlapLabels)
    return;
end

switch lower(config.separation_defect_count_source)
    case 'curvature'
        defectLabels = [];
        if ~isempty(overlapPoints)
            defectLabels = [overlapPoints.Label];
        end

    case 'convexity'
        convexityDefects = DetectConvexityDefects(BW, ...
            config.convexity_defect_depth_threshold, false);
        defectLabels = [];
        if ~isempty(convexityDefects)
            defectLabels = [convexityDefects.Label];
        end

    otherwise
        error('SeparateOverlappedObjectsByStrategy:InvalidDefectSource', ...
            'separation_defect_count_source must be ''curvature'' or ''convexity''.');
end

watershedMask = false(size(BW));

for labelIndex = 1:numel(overlapLabels)
    labelNumber = overlapLabels(labelIndex);
    defectCount = nnz(defectLabels == labelNumber);

    infoIndex = numel(separationInfo) + 1;
    separationInfo(infoIndex).Label = labelNumber;
    separationInfo(infoIndex).DefectCount = defectCount;
    separationInfo(infoIndex).CurveCoefficients = [];
    separationInfo(infoIndex).CurveCenter = [];
    separationInfo(infoIndex).CurveScale = [];
    separationInfo(infoIndex).UsedFallback = false;

    if defectCount == config.two_bubble_defect_count
        watershedMask = watershedMask | (L == labelNumber);
        separationInfo(infoIndex).Method = 'watershed_two_defects';
    elseif defectCount >= config.multiple_bubble_min_defects
        switch lower(config.multiple_overlap_separation_method)
            case 'watershed'
                watershedMask = watershedMask | (L == labelNumber);
                separationInfo(infoIndex).Method = 'watershed_multiple_defects';
            case 'skip'
                separationInfo(infoIndex).Method = 'skipped_multiple_defects';
            otherwise
                error('SeparateOverlappedObjectsByStrategy:InvalidMultipleMethod', ...
                    'multiple_overlap_separation_method must be ''watershed'' or ''skip''.');
        end
    else
        % Overlapped clusters with an ambiguous defect count (e.g. 0 or 1
        % detected neck points) are still real overlapped objects. Send them
        % through watershed so they are not lost from the final result.
        watershedMask = watershedMask | (L == labelNumber);
        separationInfo(infoIndex).Method = 'watershed_unclassified';
    end
end

if any(watershedMask(:))
    newImage = MaskOverlappedBubblesBeforeWatershed(I_Formask, watershedMask, config);
    save('newImage.mat', 'newImage'); % save the masked cluster image as a matrix
    watershedLabels = WatershedMyself(newImage, config); % This function return only the overlapped objects.
    watershedObjects = CombineWatershedSegments(watershedLabels, config);
    FinalSegmentedObjects = appendLabels(FinalSegmentedObjects, watershedObjects);
end

FinalSegmentedObjects = ReorderingLabels(FinalSegmentedObjects);

end

function combinedLabels = appendLabels(combinedLabels, labelsToAdd)
if isempty(labelsToAdd) || ~any(labelsToAdd(:))
    return;
end

labelsToAdd = ReorderingLabels(labelsToAdd);
offset = max(combinedLabels(:));
addMask = labelsToAdd > 0;
combinedLabels(addMask) = labelsToAdd(addMask) + offset;
end
