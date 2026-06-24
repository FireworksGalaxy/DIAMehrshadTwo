function [overlapPoints, overlayImage, overlapMask] = DetectCurvatureOverlapPoints(inputImage, config, backgroundImage, Bpre, Lpre)
% DetectCurvatureOverlapPoints - Detect overlap necks from sharp concave turns.
%
% Usage:
%   [points, overlay, mask] = DetectCurvatureOverlapPoints(BW, config, I_Formask, B, L)
%
% Output:
%   overlapPoints - Struct array with Label, Coordinates [x y], Depth,
%                   Angle and AngleDegrees.
%   overlayImage  - RGB image with detected overlap points marked in red.
%   overlapMask   - Binary mask containing objects that have at least one point.

if nargin < 2 || isempty(config)
    config = BubbleDetectionConfig();
end

if nargin < 3
    backgroundImage = [];
end

if nargin >= 5 && ~isempty(Bpre) && ~isempty(Lpre)
    boundaries = Bpre;
    labelImage = Lpre;
elseif islogical(inputImage) || max(inputImage(:)) <= 1
    [boundaries, labelImage] = bwboundaries(logical(inputImage), 'noholes');
else
    labelImage = inputImage;
    boundaries = bwboundaries(labelImage > 0, 'noholes');
end

windowSize = max(2, round(config.curvature_window_size));
angleThresholdDegrees = config.curvature_angle_threshold_degrees;
minDepth = config.curvature_min_depth_pixels;
minSeparation = config.curvature_min_point_separation;

overlapPoints = struct('Label', {}, 'Coordinates', {}, 'Depth', {}, ...
    'Angle', {}, 'AngleDegrees', {}, 'BoundaryIndex', {});
overlapLabels = [];
centroids = regionprops(labelImage, 'Centroid');

for labelNumber = 1:numel(boundaries)
    boundary = boundaries{labelNumber};
    if size(boundary, 1) > 1 && isequal(boundary(1, :), boundary(end, :))
        boundary(end, :) = [];
    end

    numberOfPoints = size(boundary, 1);
    if numberOfPoints < (2 * windowSize + 3)
        continue;
    end

    xBoundary = boundary(:, 2);
    yBoundary = boundary(:, 1);

    if labelNumber > numel(centroids) || isempty(centroids(labelNumber).Centroid)
        continue;
    end
    centroid = centroids(labelNumber).Centroid;

    previousIndex = mod((1:numberOfPoints)' - windowSize - 1, numberOfPoints) + 1;
    nextIndex = mod((1:numberOfPoints)' + windowSize - 1, numberOfPoints) + 1;

    previousVector = [xBoundary(previousIndex) - xBoundary, yBoundary(previousIndex) - yBoundary];
    nextVector = [xBoundary(nextIndex) - xBoundary, yBoundary(nextIndex) - yBoundary];
    previousLength = hypot(previousVector(:, 1), previousVector(:, 2));
    nextLength = hypot(nextVector(:, 1), nextVector(:, 2));
    validAngle = previousLength > 0 & nextLength > 0;

    dotValue = previousVector(:, 1) .* nextVector(:, 1) + previousVector(:, 2) .* nextVector(:, 2);
    cosAngle = zeros(numberOfPoints, 1);
    cosAngle(validAngle) = dotValue(validAngle) ./ (previousLength(validAngle) .* nextLength(validAngle));
    cosAngle = max(-1, min(1, cosAngle));
    includedAngle = acosd(cosAngle);

    chordLength = hypot(xBoundary(nextIndex) - xBoundary(previousIndex), ...
        yBoundary(nextIndex) - yBoundary(previousIndex));
    validChord = chordLength > 0;
    depth = zeros(numberOfPoints, 1);
    depth(validChord) = abs((xBoundary(nextIndex(validChord)) - xBoundary(previousIndex(validChord))) .* ...
        (yBoundary(previousIndex(validChord)) - yBoundary(validChord)) - ...
        (xBoundary(previousIndex(validChord)) - xBoundary(validChord)) .* ...
        (yBoundary(nextIndex(validChord)) - yBoundary(previousIndex(validChord)))) ./ chordLength(validChord);

    chordVectorX = xBoundary(nextIndex) - xBoundary(previousIndex);
    chordVectorY = yBoundary(nextIndex) - yBoundary(previousIndex);
    pointSide = chordVectorX .* (yBoundary - yBoundary(previousIndex)) - ...
        chordVectorY .* (xBoundary - xBoundary(previousIndex));
    centroidSide = chordVectorX .* (centroid(2) - yBoundary(previousIndex)) - ...
        chordVectorY .* (centroid(1) - xBoundary(previousIndex));
    isInwardV = sign(pointSide) == sign(centroidSide) & pointSide ~= 0 & centroidSide ~= 0;

    candidateIndex = find(validAngle & isInwardV & includedAngle <= angleThresholdDegrees & depth >= minDepth);
    if isempty(candidateIndex)
        continue;
    end

    [~, sortOrder] = sort(depth(candidateIndex), 'descend');
    candidateIndex = candidateIndex(sortOrder);
    keptIndex = [];
    for candidate = reshape(candidateIndex, 1, [])
        if isempty(keptIndex)
            keptIndex = candidate;
        else
            boundaryDistance = abs(candidate - keptIndex);
            circularDistance = min(boundaryDistance, numberOfPoints - boundaryDistance);
            if all(circularDistance >= minSeparation)
                keptIndex(end + 1) = candidate;
            end
        end
    end

    overlapLabels(end + 1) = labelNumber;
    for pointIndex = reshape(keptIndex, 1, [])
        pointCount = numel(overlapPoints) + 1;
        overlapPoints(pointCount).Label = labelNumber;
        overlapPoints(pointCount).Coordinates = [xBoundary(pointIndex), yBoundary(pointIndex)];
        overlapPoints(pointCount).Depth = depth(pointIndex);
        overlapPoints(pointCount).Angle = deg2rad(includedAngle(pointIndex));
        overlapPoints(pointCount).AngleDegrees = includedAngle(pointIndex);
        overlapPoints(pointCount).BoundaryIndex = pointIndex;
    end
end

if ~isempty(overlapPoints)
    [~, order] = sort([overlapPoints.Depth], 'descend');
    overlapPoints = overlapPoints(order);
end

if nargout > 2
    overlapMask = ismember(labelImage, unique(overlapLabels));
else
    overlapMask = [];
end

if nargout > 1
    overlayImage = createPointOverlay(labelImage, backgroundImage, overlapPoints, config.curvature_overlay_point_radius);
else
    overlayImage = [];
end

end

function overlayImage = createPointOverlay(labelImage, backgroundImage, overlapPoints, radius)
if isempty(backgroundImage)
    baseImage = mat2gray(labelImage > 0);
else
    baseImage = backgroundImage;
end

if ndims(baseImage) == 2
    overlayImage = repmat(mat2gray(baseImage), 1, 1, 3);
else
    overlayImage = im2double(baseImage);
end

if isempty(overlapPoints)
    return;
end

[imageHeight, imageWidth, ~] = size(overlayImage);
[colOffset, rowOffset] = meshgrid(-radius:radius, -radius:radius);
insideDisk = rowOffset.^2 + colOffset.^2 <= radius^2;
rowOffset = rowOffset(insideDisk);
colOffset = colOffset(insideDisk);

allRows = [];
allCols = [];
for pointIndex = 1:numel(overlapPoints)
    point = overlapPoints(pointIndex).Coordinates;
    allRows = [allRows; point(2) + rowOffset(:)];
    allCols = [allCols; point(1) + colOffset(:)];
end

allRows = round(allRows);
allCols = round(allCols);
valid = allRows >= 1 & allRows <= imageHeight & allCols >= 1 & allCols <= imageWidth;
pixelIndices = sub2ind([imageHeight, imageWidth], allRows(valid), allCols(valid));

redChannel = overlayImage(:, :, 1);
greenChannel = overlayImage(:, :, 2);
blueChannel = overlayImage(:, :, 3);
redChannel(pixelIndices) = 1;
greenChannel(pixelIndices) = 0;
blueChannel(pixelIndices) = 0;
overlayImage(:, :, 1) = redChannel;
overlayImage(:, :, 2) = greenChannel;
overlayImage(:, :, 3) = blueChannel;
end
