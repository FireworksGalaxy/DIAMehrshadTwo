
 function PostProcessing(FinalSegmentedObjects,NonSegOldobject,I_Formask,Object_information_last,config,TwoOverlapMask,MultipleOverlapMask)
% PostProcessing - Final visualization and reporting of the detected bubbles.
%
% Usage:
%   PostProcessing(FinalSegmentedObjects, NonSegOldobject, I_Formask, ...
%                  Object_information_last, config)
%
% Input:
%   FinalSegmentedObjects   - label image of the watershed-separated bubbles.
%   NonSegOldobject         - label image of the single (non-overlapped) bubbles.
%   I_Formask               - original RGB image, used as the display background.
%   Object_information_last  - per-object properties of the final combined image.
%   config                  - configuration struct.
%   TwoOverlapMask          - optional mask of original clusters with exactly
%                             two overlap/neck points.
%   MultipleOverlapMask     - optional mask of original clusters with three or
%                             more overlap/neck points.
%
% Produces an overlay of single bubbles, two-bubble overlaps, and multi-bubble
% overlaps using the colors configured in BubbleDetectionConfig, plus a
% Sauter-diameter plot and bubble-size-distribution histogram. Diameters are scaled by
% config.calibration_divisor.
% Segobject=load('FinalSegmentObject.mat');
% Segobject=Segobject.FinalSegmentedObjects;  
if nargin < 6 || isempty(TwoOverlapMask)
    TwoOverlapMask = false(size(FinalSegmentedObjects));
end
if nargin < 7 || isempty(MultipleOverlapMask)
    MultipleOverlapMask = false(size(FinalSegmentedObjects));
end
showObjectSeparationOverlay = config.show_figures;
if isfield(config, 'show_object_separation_overlay')
    showObjectSeparationOverlay = showObjectSeparationOverlay || config.show_object_separation_overlay;
end

if showObjectSeparationOverlay
    seg_obj_info=GetObjectProperties(FinalSegmentedObjects,config);
    % 
    % NonSegOldobject=load('LabelSeg.mat');
    % NonSegOldobject=NonSegOldobject.Labelseg;
    all_obj_info=GetObjectProperties(NonSegOldobject,config);
    % I_Formask = imread('Picture_3rd.jpg');
    % I_Formask = imcrop(I_Formask,[80 200 1460 1400]);
    figure('Name','PostProcessFunction - object seperation');
    imshow(rgb2gray(I_Formask));
    hold on;

    for k = 1:length(all_obj_info)
        if (all_obj_info(k).Roundness<config.roundness_threshold_display)
            boundaryWasPlotted = plotLabelBoundary(NonSegOldobject, k, all_obj_info(k).BoundingBox, ...
                config.single_bubble_display_color, config.object_separation_line_width);
            if boundaryWasPlotted && config.show_object_separation_numbers
                c = all_obj_info(k).Centroid;
                text(c(1), c(2), sprintf('%d', k), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle','Color','white');
            end
        end
    end

    segPixelStats = regionprops(FinalSegmentedObjects,'PixelIdxList');
    twoOverlapObjects = false(1, length(segPixelStats));
    multipleOverlapObjects = false(1, length(segPixelStats));
    for k = 1:length(segPixelStats)
        if isempty(segPixelStats(k).PixelIdxList)
            continue;
        end
        twoOverlapObjects(k) = any(TwoOverlapMask(segPixelStats(k).PixelIdxList));
        multipleOverlapObjects(k) = any(MultipleOverlapMask(segPixelStats(k).PixelIdxList));
    end

    redObjectCount = 0;
    greenObjectCount = 0;
    for k = 1:length(seg_obj_info)
        if (seg_obj_info(k).Roundness<config.roundness_threshold_display)
            boundaryColor = config.two_overlap_display_color;
            if multipleOverlapObjects(k)
                boundaryColor = config.multiple_overlap_display_color;
            end
            boundaryWasPlotted = plotLabelBoundary(FinalSegmentedObjects, k, seg_obj_info(k).BoundingBox, ...
                boundaryColor, config.object_separation_line_width);
            if boundaryWasPlotted
                if multipleOverlapObjects(k)
                    greenObjectCount = greenObjectCount + 1;
                elseif twoOverlapObjects(k)
                    redObjectCount = redObjectCount + 1;
                end
            end
            if boundaryWasPlotted && config.show_object_separation_numbers
                c = seg_obj_info(k).Centroid;
                text(c(1), c(2), sprintf('%d', k), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle','Color','white');
            end
        end
    end
    if config.verbose
        fprintf('Object separation overlay: %d segmented objects drawn as two-overlap color.\n', redObjectCount);
        fprintf('Object separation overlay: %d segmented objects drawn as multiple-overlap color.\n', greenObjectCount);
    end
    hold off;

end

showSummaryPlots = ~isfield(config, 'show_postprocessing_summary_plots') || config.show_postprocessing_summary_plots;
if showSummaryPlots
    figure('Name','PostProcessFunction - Sauter Diameter');
    diameters = [Object_information_last.SauterDia]/config.calibration_divisor;
    plot(diameters);
    pause(0.1);
    
    figure('Name','PostProcessFunction - BSD');
    histogram(diameters,50);
end

end

function boundaryWasPlotted = plotLabelBoundary(labelImage, labelIndex, boundingBox, boundaryColor, lineWidth)
boundaryWasPlotted = false;

rowStart = max(1, floor(boundingBox(2) + 0.5));
colStart = max(1, floor(boundingBox(1) + 0.5));
rowEnd = min(size(labelImage, 1), ceil(boundingBox(2) + boundingBox(4) - 0.5));
colEnd = min(size(labelImage, 2), ceil(boundingBox(1) + boundingBox(3) - 0.5));

if rowStart > rowEnd || colStart > colEnd
    return;
end

objectMask = labelImage(rowStart:rowEnd, colStart:colEnd) == labelIndex;
objectBoundaries = bwboundaries(objectMask, 'noholes');
for boundaryIndex = 1:length(objectBoundaries)
    boundary = objectBoundaries{boundaryIndex};
    plot(boundary(:,2) + colStart - 1, boundary(:,1) + rowStart - 1, boundaryColor, 'LineWidth', lineWidth)
    boundaryWasPlotted = true;
end
end
