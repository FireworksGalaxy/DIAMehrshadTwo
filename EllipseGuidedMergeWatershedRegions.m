function LL = EllipseGuidedMergeWatershedRegions(LL, config)
% EllipseGuidedMergeWatershedRegions - Merge neighboring watershed regions
% when their combined boundary is better explained by one ellipse.
%
% Usage:
%   LL = EllipseGuidedMergeWatershedRegions(LL, config)
%
% The merge score uses boundary points. Each candidate pair is evaluated by
% fitting ellipses to region A, region B, and the merged A+B region. The error
% is a scale-free radial RMSE: points on the fitted ellipse have residual 0.

LL = ReorderingLabels(LL);
showFigures = isfield(config, 'show_ellipse_merge_figures') && config.show_ellipse_merge_figures;
fitFigureCount = 0;
mergeCount = 0;
verbose = isfield(config, 'verbose') && config.verbose;

if verbose
    fprintf('Ellipse merge: starting with %d watershed regions.\n', max(LL(:)));
end

if showFigures
    figure('Name', 'Ellipse merge - original watershed result');
    imshow(label2rgb(LL, 'jet', 'k', 'shuffle'));
end

for iteration = 1:config.merge_max_iterations
    candidatePairs = findNeighboringRegionPairs(LL, config);
    if isempty(candidatePairs)
        if verbose
            fprintf('Ellipse merge iteration %d: no neighboring candidate pairs.\n', iteration);
        end
        break;
    end

    if showFigures && iteration == 1
        showCandidatePairs(LL, candidatePairs, config);
    end

    labelStats = regionprops(LL, 'BoundingBox');
    fitCache = buildFitCache(LL, labelStats, candidatePairs, config);
    bestDecision = struct('Accepted', false, 'Improvement', -inf);
    bestRejectedDecision = struct('Accepted', false, 'Improvement', -inf, ...
        'LabelA', NaN, 'LabelB', NaN, 'SeparateError', inf, 'MergedError', inf);
    for pairIndex = 1:size(candidatePairs, 1)
        decision = evaluateMergeCandidate(LL, candidatePairs(pairIndex, 1), candidatePairs(pairIndex, 2), config, fitCache, labelStats);
        if decision.Accepted && decision.Improvement > bestDecision.Improvement
            bestDecision = decision;
        elseif ~decision.Accepted && decision.Improvement > bestRejectedDecision.Improvement
            bestRejectedDecision = decision;
        end
    end

    if ~bestDecision.Accepted
        if verbose
            fprintf(['Ellipse merge iteration %d: %d candidates, no accepted merge. ', ...
                'Best rejected pair %g/%g: separate %.4f, merged %.4f, improvement %.4f.\n'], ...
                iteration, size(candidatePairs, 1), bestRejectedDecision.LabelA, bestRejectedDecision.LabelB, ...
                bestRejectedDecision.SeparateError, bestRejectedDecision.MergedError, bestRejectedDecision.Improvement);
        end
        break;
    end

    mergeCount = mergeCount + 1;
    if verbose
        fprintf(['Ellipse merge iteration %d: %d candidates, accepted pair %d/%d ', ...
            '(separate %.4f, merged %.4f, improvement %.4f).\n'], ...
            iteration, size(candidatePairs, 1), bestDecision.LabelA, bestDecision.LabelB, ...
            bestDecision.SeparateError, bestDecision.MergedError, bestDecision.Improvement);
    end

    if showFigures && fitFigureCount < config.ellipse_merge_max_fit_figures
        fitFigureCount = fitFigureCount + 1;
        showEllipseMergeFit(LL, bestDecision, iteration);
    end

    LL = MergeLabels(LL, bestDecision.LabelA, bestDecision.LabelB);
    LL = ReorderingLabels(LL);
end

if showFigures
    figure('Name', 'Ellipse merge - final merged result');
    imshow(label2rgb(LL, 'jet', 'k', 'shuffle'));
end

if verbose
    fprintf('Ellipse merge: accepted %d merge(s), final region count %d.\n', mergeCount, max(LL(:)));
end

end

function candidatePairs = findNeighboringRegionPairs(LL, config)
candidatePairs = zeros(0, 2);
pairScores = zeros(0, 1);

radius = 2;
if isfield(config, 'ellipse_merge_neighbor_dilation_radius')
    radius = config.ellipse_merge_neighbor_dilation_radius;
end
neighborhood = true(2 * radius + 1);

ridgePixels = find(LL == 0 & imdilate(LL > 0, neighborhood));
[ridgeRows, ridgeCols] = ind2sub(size(LL), ridgePixels);

for pixelIndex = 1:numel(ridgePixels)
    rowStart = max(1, ridgeRows(pixelIndex) - radius);
    rowEnd = min(size(LL, 1), ridgeRows(pixelIndex) + radius);
    colStart = max(1, ridgeCols(pixelIndex) - radius);
    colEnd = min(size(LL, 2), ridgeCols(pixelIndex) + radius);

    neighborLabels = unique(LL(rowStart:rowEnd, colStart:colEnd));
    neighborLabels(neighborLabels == 0) = [];
    if numel(neighborLabels) < 2
        continue;
    end

    [labelA, labelB] = find(triu(true(numel(neighborLabels)), 1));
    newPairs = [neighborLabels(labelA), neighborLabels(labelB)];
    newScores = ones(size(newPairs, 1), 1);
    candidatePairs = [candidatePairs; newPairs]; %#ok<AGROW>
    pairScores = [pairScores; newScores]; %#ok<AGROW>

    if isfield(config, 'ellipse_merge_max_candidate_pairs') && ...
            size(unique(candidatePairs, 'rows'), 1) >= 3 * config.ellipse_merge_max_candidate_pairs
        break;
    end
end

if ~isempty(candidatePairs)
    [candidatePairs, ~, groupIndex] = unique(candidatePairs, 'rows');
    pairScores = accumarray(groupIndex, pairScores);
    [~, scoreOrder] = sort(pairScores, 'descend');
    candidatePairs = candidatePairs(scoreOrder, :);
    if isfield(config, 'ellipse_merge_max_candidate_pairs')
        candidatePairs = candidatePairs(1:min(end, config.ellipse_merge_max_candidate_pairs), :);
    end
end
end

function fitCache = buildFitCache(LL, labelStats, candidatePairs, config)
maxLabel = max(LL(:));
fitCache = cell(maxLabel, 1);
candidateLabels = unique(candidatePairs(:));
for labelIndex = 1:numel(candidateLabels)
    labelNumber = candidateLabels(labelIndex);
    if labelNumber <= numel(labelStats) && labelStats(labelNumber).BoundingBox(3) > 0
        fitCache{labelNumber} = fitEllipseToLabelCrop(LL, labelNumber, labelStats(labelNumber).BoundingBox, config);
    end
end
end

function decision = evaluateMergeCandidate(LL, labelA, labelB, config, fitCache, labelStats)
fitA = getCachedFit(fitCache, labelA);
fitB = getCachedFit(fitCache, labelB);

fitMerged = fitEllipseToMergedPair(LL, labelA, labelB, labelStats, config);

separateError = inf;
if fitA.IsValid && fitB.IsValid
    separateError = (fitA.Error * fitA.NumPoints + fitB.Error * fitB.NumPoints) / ...
        (fitA.NumPoints + fitB.NumPoints);
end

function fit = getCachedFit(fitCache, labelNumber)
if labelNumber <= numel(fitCache) && ~isempty(fitCache{labelNumber})
    fit = fitCache{labelNumber};
else
    fit = emptyEllipseFit();
end
end

improvement = separateError - fitMerged.Error;
accepted = fitMerged.IsValid && isfinite(separateError) && ...
    fitMerged.Error < separateError * config.ellipse_merge_error_ratio && ...
    improvement >= config.ellipse_merge_min_improvement && ...
    fitMerged.Error <= config.ellipse_merge_max_merged_error;

decision = struct();
decision.LabelA = labelA;
decision.LabelB = labelB;
decision.FitA = fitA;
decision.FitB = fitB;
decision.FitMerged = fitMerged;
decision.SeparateError = separateError;
decision.MergedError = fitMerged.Error;
decision.Improvement = improvement;
decision.Accepted = accepted;
end

function fit = fitEllipseToLabelCrop(LL, labelNumber, boundingBox, config)
[rowStart, rowEnd, colStart, colEnd] = boundingBoxToLimits(boundingBox, size(LL), 1);
localMask = LL(rowStart:rowEnd, colStart:colEnd) == labelNumber;
fit = fitEllipseToMask(localMask, config, colStart - 1, rowStart - 1);
end

function fit = fitEllipseToMergedPair(LL, labelA, labelB, labelStats, config)
if labelA > numel(labelStats) || labelB > numel(labelStats)
    fit = emptyEllipseFit();
    return;
end

boxA = labelStats(labelA).BoundingBox;
boxB = labelStats(labelB).BoundingBox;
colMin = min(boxA(1), boxB(1));
rowMin = min(boxA(2), boxB(2));
colMax = max(boxA(1) + boxA(3), boxB(1) + boxB(3));
rowMax = max(boxA(2) + boxA(4), boxB(2) + boxB(4));
mergedBox = [colMin, rowMin, colMax - colMin, rowMax - rowMin];

[rowStart, rowEnd, colStart, colEnd] = boundingBoxToLimits(mergedBox, size(LL), 2);
localLabels = LL(rowStart:rowEnd, colStart:colEnd);
maskA = localLabels == labelA;
maskB = localLabels == labelB;
bridgeMask = imdilate(maskA, true(3)) & imdilate(maskB, true(3));
mergedMask = maskA | maskB | bridgeMask;
fit = fitEllipseToMask(mergedMask, config, colStart - 1, rowStart - 1);
end

function [rowStart, rowEnd, colStart, colEnd] = boundingBoxToLimits(boundingBox, imageSize, padding)
rowStart = max(1, floor(boundingBox(2)) - padding);
colStart = max(1, floor(boundingBox(1)) - padding);
rowEnd = min(imageSize(1), ceil(boundingBox(2) + boundingBox(4)) + padding);
colEnd = min(imageSize(2), ceil(boundingBox(1) + boundingBox(3)) + padding);
end

function fit = fitEllipseToMask(regionMask, config, xOffset, yOffset)
fit = emptyEllipseFit();

boundaries = bwboundaries(regionMask, 'noholes');
if isempty(boundaries)
    return;
end

points = [];
for boundaryIndex = 1:numel(boundaries)
    boundary = boundaries{boundaryIndex};
    points = [points; boundary(:, 2) + xOffset, boundary(:, 1) + yOffset]; %#ok<AGROW>
end

if size(points, 1) < 6
    return;
end

if isfield(config, 'ellipse_merge_max_boundary_points') && ...
        size(points, 1) > config.ellipse_merge_max_boundary_points
    sampleIndex = round(linspace(1, size(points, 1), config.ellipse_merge_max_boundary_points));
    points = points(sampleIndex, :);
end

center = mean(points, 1);
centered = points - center;
covarianceMatrix = (centered' * centered) / max(size(centered, 1) - 1, 1);
[vectors, values] = eig(covarianceMatrix);
[~, order] = sort(diag(values), 'descend');
vectors = vectors(:, order);

rotated = centered * vectors;
axisA = sqrt(2) * sqrt(mean(rotated(:, 1).^2));
axisB = sqrt(2) * sqrt(mean(rotated(:, 2).^2));

if axisA <= eps || axisB <= eps || ~isfinite(axisA) || ~isfinite(axisB)
    return;
end

radialDistance = sqrt((rotated(:, 1) ./ axisA).^2 + (rotated(:, 2) ./ axisB).^2);
residuals = radialDistance - 1;

fit.IsValid = true;
fit.Center = center;
fit.Axes = [axisA, axisB];
fit.Rotation = atan2(vectors(2, 1), vectors(1, 1));
fit.Error = sqrt(mean(residuals.^2));
fit.NumPoints = size(points, 1);
fit.Points = points;
end

function fit = emptyEllipseFit()
fit = struct('IsValid', false, 'Center', [NaN, NaN], 'Axes', [NaN, NaN], ...
    'Rotation', NaN, 'Error', inf, 'NumPoints', 0, 'Points', zeros(0, 2));
end

function showCandidatePairs(LL, candidatePairs, config)
figure('Name', 'Ellipse merge - candidate neighboring pairs');
imshow(label2rgb(LL, 'jet', 'k', 'shuffle'));
hold on;

stats = regionprops(LL, 'Centroid');
maxPairs = min(size(candidatePairs, 1), config.ellipse_merge_max_visual_pairs);
for pairIndex = 1:maxPairs
    labelA = candidatePairs(pairIndex, 1);
    labelB = candidatePairs(pairIndex, 2);
    if labelA > numel(stats) || labelB > numel(stats)
        continue;
    end
    centroidA = stats(labelA).Centroid;
    centroidB = stats(labelB).Centroid;
    plot([centroidA(1), centroidB(1)], [centroidA(2), centroidB(2)], 'w-', 'LineWidth', 1);
end
hold off;
end

function showEllipseMergeFit(LL, decision, iteration)
localMask = (LL == decision.LabelA) | (LL == decision.LabelB);
boundingBox = regionprops(localMask, 'BoundingBox');
if isempty(boundingBox)
    return;
end

box = boundingBox(1).BoundingBox;
rowStart = max(1, floor(box(2)) - 10);
colStart = max(1, floor(box(1)) - 10);
rowEnd = min(size(LL, 1), ceil(box(2) + box(4)) + 10);
colEnd = min(size(LL, 2), ceil(box(1) + box(3)) + 10);

cropLabels = LL(rowStart:rowEnd, colStart:colEnd);
figure('Name', sprintf('Ellipse merge - accepted pair %d/%d, iteration %d', ...
    decision.LabelA, decision.LabelB, iteration));

subplot(1, 2, 1);
imshow(label2rgb(cropLabels, 'jet', 'k', 'shuffle'));
title(sprintf('Before: %.3f + %.3f', decision.FitA.Error, decision.FitB.Error));
hold on;
plotEllipse(decision.FitA, colStart, rowStart, 'w-');
plotEllipse(decision.FitB, colStart, rowStart, 'c-');
hold off;

afterLabels = cropLabels;
afterLabels(afterLabels == decision.LabelA) = decision.LabelB;
subplot(1, 2, 2);
imshow(label2rgb(afterLabels, 'jet', 'k', 'shuffle'));
title(sprintf('After: %.3f', decision.FitMerged.Error));
hold on;
plotEllipse(decision.FitMerged, colStart, rowStart, 'w-');
hold off;
end

function plotEllipse(fit, colStart, rowStart, lineStyle)
if ~fit.IsValid
    return;
end

theta = linspace(0, 2 * pi, 200);
ellipsePoints = [fit.Axes(1) * cos(theta); fit.Axes(2) * sin(theta)];
rotationMatrix = [cos(fit.Rotation), -sin(fit.Rotation); sin(fit.Rotation), cos(fit.Rotation)];
rotatedPoints = rotationMatrix * ellipsePoints;
x = rotatedPoints(1, :) + fit.Center(1) - colStart + 1;
y = rotatedPoints(2, :) + fit.Center(2) - rowStart + 1;
plot(x, y, lineStyle, 'LineWidth', 1.5);
end
