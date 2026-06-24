function [FinalSegmentedObjects, separationInfo] = SeparateOverlappedObjectsByStrategy(BW, I_Formask, RealOverlap, L, overlapPoints, config)
% SeparateOverlappedObjectsByStrategy - Route overlapped clusters to split methods.
%
% Two-defect clusters are split by watershed. Multiple-defect clusters are
% split by config.multiple_overlap_separation_method: 'watershed' or the
% experimental implicit cubic curve fit.

FinalSegmentedObjects = zeros(size(BW));
separationInfo = struct('Label', {}, 'DefectCount', {}, 'Method', {}, ...
    'CurveCoefficients', {}, 'CurveCenter', {}, 'CurveScale', {}, 'UsedFallback', {});
strategyTimer = tic;

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
        defectTimer = tic;
        convexityDefects = DetectConvexityDefects(BW, ...
            config.convexity_defect_depth_threshold, false);
        if config.verbose
            fprintf('Separate strategy: convexity defect counting %.2f s.\n', toc(defectTimer));
        end
        defectLabels = [];
        if ~isempty(convexityDefects)
            defectLabels = [convexityDefects.Label];
        end

    otherwise
        error('SeparateOverlappedObjectsByStrategy:InvalidDefectSource', ...
            'separation_defect_count_source must be ''curvature'' or ''convexity''.');
end

watershedMask = false(size(BW));
curveFitLabels = [];

for labelIndex = 1:numel(overlapLabels)
    labelNumber = overlapLabels(labelIndex);
    defectCount = nnz(defectLabels == labelNumber);
    clusterMask = L == labelNumber;

    infoIndex = numel(separationInfo) + 1;
    separationInfo(infoIndex).Label = labelNumber;
    separationInfo(infoIndex).DefectCount = defectCount;
    separationInfo(infoIndex).CurveCoefficients = [];
    separationInfo(infoIndex).CurveCenter = [];
    separationInfo(infoIndex).CurveScale = [];
    separationInfo(infoIndex).UsedFallback = false;

    if defectCount == config.two_bubble_defect_count
        watershedMask = watershedMask | clusterMask;
        separationInfo(infoIndex).Method = 'watershed_two_defects';
    elseif defectCount >= config.multiple_bubble_min_defects
        switch lower(config.multiple_overlap_separation_method)
            case 'watershed'
                watershedMask = watershedMask | clusterMask;
                separationInfo(infoIndex).Method = 'watershed_multiple_defects';
            case 'curve_fit'
                curveFitLabels(end + 1) = labelNumber;
                separationInfo(infoIndex).Method = 'curve_fit_multiple_defects';
            otherwise
                error('SeparateOverlappedObjectsByStrategy:InvalidMethod', ...
                    'multiple_overlap_separation_method must be ''watershed'' or ''curve_fit''.');
        end
    else
        watershedMask = watershedMask | clusterMask;
        separationInfo(infoIndex).Method = 'watershed_unclassified';
    end
end

if ~isempty(curveFitLabels)
    curveTimer = tic;
    [curveLabels, curveInfo, failedCurveLabels] = SeparateMultipleOverlapsByImplicitCubic(BW, L, curveFitLabels, config);
    if config.verbose
        fprintf('Separate strategy: curve-fit separation %.2f s.\n', toc(curveTimer));
    end
    FinalSegmentedObjects = appendLabels(FinalSegmentedObjects, curveLabels);

    for curveIndex = 1:numel(curveInfo)
        infoIndex = find([separationInfo.Label] == curveInfo(curveIndex).Label, 1);
        if ~isempty(infoIndex)
            separationInfo(infoIndex).CurveCoefficients = curveInfo(curveIndex).Coefficients;
            separationInfo(infoIndex).CurveCenter = curveInfo(curveIndex).Center;
            separationInfo(infoIndex).CurveScale = curveInfo(curveIndex).Scale;
        end
    end

    if config.curve_fit_fallback_to_watershed && ~isempty(failedCurveLabels)
        watershedMask = watershedMask | ismember(L, failedCurveLabels);
        for failedLabel = reshape(failedCurveLabels, 1, [])
            infoIndex = find([separationInfo.Label] == failedLabel, 1);
            if ~isempty(infoIndex)
                separationInfo(infoIndex).Method = 'curve_fit_failed_watershed_fallback';
                separationInfo(infoIndex).UsedFallback = true;
            end
        end
    end
end

if any(watershedMask(:))
    maskTimer = tic;
    newImage = MaskOverlappedBubblesBeforeWatershed(I_Formask, watershedMask, config);
    if config.verbose
        fprintf('Separate strategy: watershed masking %.2f s.\n', toc(maskTimer));
    end

    watershedTimer = tic;
    watershedLabels = WatershedMyself(newImage, config);
    if config.verbose
        fprintf('Separate strategy: watershed transform %.2f s.\n', toc(watershedTimer));
    end

    combineTimer = tic;
    watershedObjects = CombineWatershedSegments(watershedLabels, config);
    if config.verbose
        fprintf('Separate strategy: watershed combine/merge %.2f s.\n', toc(combineTimer));
    end

    appendTimer = tic;
    FinalSegmentedObjects = appendLabels(FinalSegmentedObjects, watershedObjects);
    if config.verbose
        fprintf('Separate strategy: append/reorder watershed objects %.2f s.\n', toc(appendTimer));
    end
end

FinalSegmentedObjects = ReorderingLabels(FinalSegmentedObjects);
if config.verbose
    fprintf('Separate strategy: total %.2f s.\n', toc(strategyTimer));
end

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

function [curveLabels, curveInfo, failedLabels] = SeparateMultipleOverlapsByImplicitCubic(BW, L, labelsToFit, config)
curveLabels = zeros(size(BW));
curveInfo = struct('Label', {}, 'Coefficients', {}, 'Center', {}, 'Scale', {});
failedLabels = [];

for labelNumber = reshape(labelsToFit, 1, [])
    clusterMask = L == labelNumber;
    boundaries = bwboundaries(clusterMask, 'noholes');
    if isempty(boundaries)
        failedLabels(end + 1) = labelNumber;
        continue;
    end

    boundary = boundaries{1};
    if size(boundary, 1) > config.curve_fit_max_boundary_points
        sampleIndex = round(linspace(1, size(boundary, 1), config.curve_fit_max_boundary_points));
        boundary = boundary(sampleIndex, :);
    end

    [rows, cols] = find(clusterMask);
    centerX = mean(cols);
    centerY = mean(rows);
    scaleValue = max([max(cols) - min(cols), max(rows) - min(rows), 1]);

    x = (boundary(:, 2) - centerX) ./ scaleValue;
    y = (boundary(:, 1) - centerY) ./ scaleValue;
    designMatrix = [x.^3, x .* y, y.^3, y, x];
    coefficients = designMatrix \ (-ones(size(x)));

    bbox = regionprops(clusterMask, 'BoundingBox');
    box = ceil(bbox(1).BoundingBox);
    colStart = max(1, box(1));
    rowStart = max(1, box(2));
    colEnd = min(size(BW, 2), colStart + box(3) - 1);
    rowEnd = min(size(BW, 1), rowStart + box(4) - 1);

    [colGrid, rowGrid] = meshgrid(colStart:colEnd, rowStart:rowEnd);
    xGrid = (colGrid - centerX) ./ scaleValue;
    yGrid = (rowGrid - centerY) ./ scaleValue;
    implicitValue = coefficients(1) .* xGrid.^3 + coefficients(2) .* xGrid .* yGrid + ...
        coefficients(3) .* yGrid.^3 + coefficients(4) .* yGrid + coefficients(5) .* xGrid + 1;

    localCluster = clusterMask(rowStart:rowEnd, colStart:colEnd);
    curveCut = abs(implicitValue) <= config.curve_fit_zero_tolerance;
    splitMask = localCluster;
    splitMask(curveCut & localCluster) = false;
    splitMask = bwareaopen(splitMask, config.curve_fit_min_component_area);
    localLabels = bwlabel(splitMask);

    if max(localLabels(:)) < 2
        failedLabels(end + 1) = labelNumber;
    else
        fullLabels = zeros(size(BW));
        fullLabels(rowStart:rowEnd, colStart:colEnd) = localLabels;
        curveLabels = appendLabels(curveLabels, fullLabels);
    end

    infoIndex = numel(curveInfo) + 1;
    curveInfo(infoIndex).Label = labelNumber;
    curveInfo(infoIndex).Coefficients = coefficients;
    curveInfo(infoIndex).Center = [centerX, centerY];
    curveInfo(infoIndex).Scale = scaleValue;
end
end
