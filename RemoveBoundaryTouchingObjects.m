function cleanedImage = RemoveBoundaryTouchingObjects(inputImage, boundaryMargin)
% RemoveBoundaryTouchingObjects - Remove objects touching the image border.
%
% Usage:
%   BWclean = RemoveBoundaryTouchingObjects(BW, 1)
%   Lclean  = RemoveBoundaryTouchingObjects(L, 5)
%
% Any connected object with at least one pixel inside the boundaryMargin-wide
% border band is removed. Use margin 1 to remove objects touching the four
% image edges; use a larger margin to remove objects near the edges.

if nargin < 2 || isempty(boundaryMargin)
    boundaryMargin = 1;
end

if boundaryMargin < 0
    error('RemoveBoundaryTouchingObjects:InvalidMargin', ...
        'boundaryMargin must be zero or a positive number.');
end

boundaryMargin = floor(boundaryMargin);
cleanedImage = inputImage;

if boundaryMargin == 0 || isempty(inputImage)
    return;
end

isBinaryInput = islogical(inputImage);
if isBinaryInput
    labelImage = bwlabel(inputImage);
else
    labelImage = inputImage;
end

[imageHeight, imageWidth] = size(labelImage);
rowMargin = min(boundaryMargin, imageHeight);
colMargin = min(boundaryMargin, imageWidth);

topLabels = labelImage(1:rowMargin, :);
bottomLabels = labelImage(imageHeight-rowMargin+1:imageHeight, :);
leftLabels = labelImage(:, 1:colMargin);
rightLabels = labelImage(:, imageWidth-colMargin+1:imageWidth);

labelsToRemove = unique([topLabels(:); bottomLabels(:); leftLabels(:); rightLabels(:)]);
labelsToRemove(labelsToRemove == 0) = [];

if isempty(labelsToRemove)
    return;
end

cleanedImage(ismember(labelImage, labelsToRemove)) = 0;

if isBinaryInput
    cleanedImage = logical(cleanedImage);
else
    cleanedImage = ReorderingLabels(cleanedImage);
end

end
