function [defects, overlayImage] = DetectConvexityDefects(inputImage, depthThreshold, showFigure, backgroundImage, Bpre, Lpre)
% DetectConvexityDefects - Find inward boundary defects from the convex hull.
%
% Usage:
%   [defects, overlay] = DetectConvexityDefects(BW, 15, true, I_Formask)
%   [defects, overlay] = DetectConvexityDefects(L, config.convexity_defect_depth_threshold)
%   [defects, overlay] = DetectConvexityDefects(BW, 15, true, I_Formask, B, L)
%
% Output:
%   defects      - Struct array with Label, Coordinates [x y], Depth,
%                  BoundaryIndex, HullStart, HullEnd and DefectBoundary.
%   overlayImage - RGB image showing boundary, hull, defect arcs and points.

if nargin < 2 || isempty(depthThreshold)
    depthThreshold = 15;
end

if nargin < 3 || isempty(showFigure)
    showFigure = false;
end

if nargin < 4
    backgroundImage = [];
end

if nargin >= 6 && ~isempty(Bpre) && ~isempty(Lpre)
    boundaries = Bpre;
    labelImage = Lpre;
elseif islogical(inputImage) || max(inputImage(:)) <= 1
    [boundaries, labelImage] = bwboundaries(logical(inputImage), 'noholes');
else
    labelImage = inputImage;
    [boundaries, ~] = bwboundaries(labelImage > 0, 'noholes');
end

defects = struct('Label', {}, 'Coordinates', {}, 'Depth', {}, ...
    'BoundaryIndex', {}, 'HullStart', {}, 'HullEnd', {}, 'DefectBoundary', {});
needsOverlay = showFigure || nargout > 1;
if needsOverlay
    overlayImage = createBaseOverlay(labelImage, backgroundImage);
else
    overlayImage = [];
end

for labelNumber = 1:numel(boundaries)
    boundary = boundaries{labelNumber};
    if size(boundary, 1) > 1 && isequal(boundary(1, :), boundary(end, :))
        boundary(end, :) = [];
    end

    if size(boundary, 1) < 4
        continue;
    end

    xBoundary = boundary(:, 2);
    yBoundary = boundary(:, 1);

    try
        hullOrder = convhull(xBoundary, yBoundary);
    catch
        continue;
    end

    hullIndices = sort(hullOrder(1:end-1));
    hullIndices = hullIndices(:);
    if numel(hullIndices) < 3
        continue;
    end

    if needsOverlay
        hullRows = yBoundary([hullIndices; hullIndices(1)]);
        hullCols = xBoundary([hullIndices; hullIndices(1)]);
        overlayImage = drawPolyline(overlayImage, boundary(:, 1), boundary(:, 2), [0 1 0]);
        overlayImage = drawPolyline(overlayImage, hullRows, hullCols, [0 0.35 1]);
    end

    for edgeIndex = 1:numel(hullIndices)
        startIndex = hullIndices(edgeIndex);
        if edgeIndex < numel(hullIndices)
            endIndex = hullIndices(edgeIndex + 1);
            arcIndices = startIndex:endIndex;
        else
            endIndex = hullIndices(1);
            arcIndices = [startIndex:numel(xBoundary), 1:endIndex];
        end

        if numel(arcIndices) < 3
            continue;
        end

        xStart = xBoundary(startIndex);
        yStart = yBoundary(startIndex);
        xEnd = xBoundary(endIndex);
        yEnd = yBoundary(endIndex);
        chordLength = hypot(xEnd - xStart, yEnd - yStart);
        if chordLength == 0
            continue;
        end

        xArc = xBoundary(arcIndices);
        yArc = yBoundary(arcIndices);
        distances = abs((xEnd - xStart) .* (yStart - yArc) - ...
            (xStart - xArc) .* (yEnd - yStart)) ./ chordLength;

        [maxDepth, maxIndex] = max(distances);
        if maxDepth < depthThreshold
            continue;
        end

        boundaryIndex = arcIndices(maxIndex);
        coordinates = [xBoundary(boundaryIndex), yBoundary(boundaryIndex)];
        defectBoundary = [xArc, yArc];

        defectCount = numel(defects) + 1;
        defects(defectCount).Label = labelNumber;
        defects(defectCount).Coordinates = coordinates;
        defects(defectCount).Depth = maxDepth;
        defects(defectCount).BoundaryIndex = boundaryIndex;
        defects(defectCount).HullStart = [xStart, yStart];
        defects(defectCount).HullEnd = [xEnd, yEnd];
        defects(defectCount).DefectBoundary = defectBoundary;

        if needsOverlay
            overlayImage = drawPolyline(overlayImage, yArc, xArc, [1 0.9 0]);
            overlayImage = drawDisk(overlayImage, coordinates(2), coordinates(1), 5, [1 0 0]);
        end
    end
end

if ~isempty(defects)
    [~, order] = sort([defects.Depth], 'descend');
    defects = defects(order);
end

if showFigure && needsOverlay
    figure('Name', 'Convexity Defect Overlap Candidates');
    imshow(overlayImage);
    hold on;
    for defectIndex = 1:numel(defects)
        point = defects(defectIndex).Coordinates;
        text(point(1) + 6, point(2), sprintf('%.1f', defects(defectIndex).Depth), ...
            'Color', 'white', 'FontWeight', 'bold');
    end
    hold off;
end

end

function overlayImage = createBaseOverlay(labelImage, backgroundImage)
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
end

function imageOut = drawPolyline(imageIn, rows, cols, color)
imageOut = imageIn;
if numel(rows) < 2
    return;
end

segmentLengths = max(abs(diff(rows)), abs(diff(cols))) + 1;
allRows = zeros(1, sum(segmentLengths));
allCols = zeros(1, sum(segmentLengths));
insertIndex = 1;
for pointIndex = 1:(numel(rows) - 1)
    segmentLength = segmentLengths(pointIndex);
    rowValues = round(linspace(rows(pointIndex), rows(pointIndex + 1), ...
        segmentLength));
    colValues = round(linspace(cols(pointIndex), cols(pointIndex + 1), segmentLength));
    endIndex = insertIndex + segmentLength - 1;
    allRows(insertIndex:endIndex) = rowValues;
    allCols(insertIndex:endIndex) = colValues;
    insertIndex = endIndex + 1;
end
imageOut = setPixels(imageOut, allRows, allCols, color);
end

function imageOut = drawDisk(imageIn, centerRow, centerCol, radius, color)
imageOut = imageIn;
[colGrid, rowGrid] = meshgrid((centerCol-radius):(centerCol+radius), ...
    (centerRow-radius):(centerRow+radius));
insideDisk = (rowGrid - centerRow).^2 + (colGrid - centerCol).^2 <= radius^2;
imageOut = setPixels(imageOut, rowGrid(insideDisk), colGrid(insideDisk), color);
end

function imageOut = setPixels(imageIn, rows, cols, color)
imageOut = imageIn;
[imageHeight, imageWidth, ~] = size(imageOut);
validPixels = rows >= 1 & rows <= imageHeight & cols >= 1 & cols <= imageWidth;
rows = rows(validPixels);
cols = cols(validPixels);

if isempty(rows)
    return;
end

pixelIndices = sub2ind([imageHeight, imageWidth], rows(:), cols(:));
for channel = 1:3
    channelImage = imageOut(:, :, channel);
    channelImage(pixelIndices) = color(channel);
    imageOut(:, :, channel) = channelImage;
end
end
