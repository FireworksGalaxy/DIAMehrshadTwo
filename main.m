% main.m - Bubble detection pipeline (entry point)
% ------------------------------------------------------------------------
% High-level flow:
%   1. Load the image and pre-process it (sharpen + adaptive binarize).
%   2. Clean the binary mask (fill holes, morphological open, remove speckle).
%   3. Measure every candidate object (area, perimeter, Sauter diameter, roundness).
%   4. Split the pipeline into two branches:
%        a) OVERLAPPED bubbles  -> separated with a watershed transform.
%        b) SINGLE bubbles      -> reconstructed directly from their boundaries.
%   5. Combine both branches into one label image and report the size distribution.
%
% All tunable parameters live in BubbleDetectionConfig.m.
% ------------------------------------------------------------------------
clc;
clear;
close all;
% Load configuration (all thresholds / radii / toggles in one place)
config = BubbleDetectionConfig();

% --- Calibration switch ------------------------------------------------------
% Use 'manual' for normal runs after calibration is known.
% Use 'image' only when you want to measure pixels/mm from a ruler image.
imageFile = 'Image/SV_3.jpeg';
calibrationMode = 'manual';  % Options: 'manual' or 'image'
manualPixelsPerMm = 26.064;      % Replace this with your measured pixels/mm value
rulerLengthMm = 10;

switch lower(calibrationMode)
    case 'image'
        calibration = CalibrateImageFromRuler(imageFile, rulerLengthMm);
        config.calibration_divisor = calibration.pixelsPerMm;

    case 'manual'
        config.calibration_divisor = manualPixelsPerMm;

    otherwise
        error('main:InvalidCalibrationMode', 'calibrationMode must be ''manual'' or ''image''.');
end

if config.verbose
    fprintf('Using calibration: %.6f pixels/mm\n', config.calibration_divisor);
end

% --- Stage 1: Load the image -------------------------------------------------
% I2 is the grayscale image used for detection; I_Formask keeps the original
% RGB image, which is needed later when masking bubbles before the watershed.
I2 = rgb2gray(imread(imageFile));
if config.show_figures
    figure ('Name','The first show of picture after loading')
    imshow(I2)
end
I_Formask=imread(imageFile);


% --- Stage 2: Sharpen + binarize --------------------------------------------
% Sharpen to emphasise bubble edges, then adaptively threshold. The result is
% inverted (~) so that bubbles become foreground (true) and the liquid is
% background (false).
I3=FastUnsharpMask(I2,config.sharpen_radius_main,config.sharpen_amount_main); % Fast (FFT) unsharp mask; radius/amount from config
if config.show_figures
    figure('Name','after sharpening the picture')
    imshow(I3)
end
BW = ~imbinarize(I3,'adaptive','ForegroundPolarity','dark','Sensitivity',config.binarize_sensitivity_main);
if config.show_figures
    figure('Name','after sharpening the picture')
    imshow(BW)
end

%% The method of the University 
% --- Stage 3: Clean up the binary mask --------------------------------------
% Fill enclosed holes, then morphologically open to break thin connections and
% smooth edges (disk radius from config).
fillBW=imfill(BW,'holes');
% figure('Name','Filling BW - After Binarization it is filled');
% imshow(fillBW)
%clearBW=imclearborder(fillBW);
seBW=strel('disk',config.morphology_disk_size); % Disk size from config
openBW=imopen(fillBW,seBW);
% figure('Name','Strel after the filling');
% imshow(openBW)

%% Extracting Data 
% --- Stage 4: Extract objects and measure their properties ------------------
BW=openBW;
BW=bwareaopen(BW,config.min_object_area); % Remove small speckle noise before analysis
if config.remove_boundary_bubbles
    BW = RemoveBoundaryTouchingObjects(BW, config.boundary_margin_pixels);
end
% Compute object boundaries/labels once and reuse them downstream
% (avoids recomputing bwboundaries several times on the same image).
[B, L] = bwboundaries(BW,'noholes');
DrawLineandTextonImage(BW,I2,config,B,L)

% Per-object geometry (Area, Perimeter, Sauter diameter, Roundness, ...).
Object_information=GetObjectProperties(BW,config,B,L);
if config.show_convexity_defects
    convexityTimer = tic;
    [convexityDefects, convexityOverlay] = DetectConvexityDefects(BW, ...
        config.convexity_defect_depth_threshold, true, I_Formask, B, L);
    if config.verbose
        fprintf('Convexity defects above threshold: %d (%.2f s)\n', ...
            numel(convexityDefects), toc(convexityTimer));
    end
end
if config.show_curvature_overlaps
    curvatureTimer = tic;
    [curvatureOverlapPoints, curvatureOverlapOverlay, curvatureOverlapMask] = ...
        DetectCurvatureOverlapPoints(BW, config, I_Formask, B, L);
    if config.verbose
        fprintf('Curvature overlap points: %d (%.2f s)\n', ...
            numel(curvatureOverlapPoints), toc(curvatureTimer));
    end
    figure('Name', 'Curvature Overlap Points');
    imshow(curvatureOverlapOverlay);
end
% Candidate overlapped bubbles are the non-round (deformed) objects.
DeformedObject=DeformedObjectFinder(BW,Object_information,config,B);
% Build a label image containing only those deformed candidates.
Extracted_object=GetExtractedObject(DeformedObject,BW,L);
% Keep only the objects that match the selected overlap criterion.
switch lower(config.overlap_detection_method)
    case 'curvature'
        if ~exist('curvatureOverlapMask', 'var')
            curvatureOverlapPoints = DetectCurvatureOverlapPoints(BW, config, [], B, L);
            curvatureOverlapMask = ismember(L, unique([curvatureOverlapPoints.Label]));
        end
        RealOverlap = curvatureOverlapMask;

    case 'roundness_area'
        RealOverlap=FindRealOverlappedObjectsAmongExtracted(Extracted_object,DeformedObject,config.area_threshold_overlapped,config);

    case 'combined'
        if ~exist('curvatureOverlapMask', 'var')
            curvatureOverlapPoints = DetectCurvatureOverlapPoints(BW, config, [], B, L);
            curvatureOverlapMask = ismember(L, unique([curvatureOverlapPoints.Label]));
        end
        roundnessOverlap = FindRealOverlappedObjectsAmongExtracted(Extracted_object,DeformedObject,config.area_threshold_overlapped,config);
        RealOverlap = curvatureOverlapMask | roundnessOverlap;

    otherwise
        error('main:InvalidOverlapDetectionMethod', ...
            'overlap_detection_method must be ''curvature'', ''roundness_area'', or ''combined''.');
end
if config.show_figures
    figure('Name','Overlapped Objects');
    imshow(RealOverlap);
end

% --- Stage 5a: OVERLAPPED branch (strategy-based splitting) -----------------
% Use defect count to choose how each overlapped cluster is separated.
if ~exist('curvatureOverlapPoints', 'var')
    curvatureOverlapPoints = DetectCurvatureOverlapPoints(BW, config, [], B, L);
end
[FinalSegmentedObjects, overlapSeparationInfo] = SeparateOverlappedObjectsByStrategy(...
    BW, I_Formask, RealOverlap, L, curvatureOverlapPoints, config);
if isempty(overlapSeparationInfo)
    twoOverlapLabels = [];
    multipleOverlapLabels = [];
else
    twoOverlapLabels = [overlapSeparationInfo([overlapSeparationInfo.DefectCount] == config.two_bubble_defect_count).Label];
    multipleOverlapLabels = [overlapSeparationInfo([overlapSeparationInfo.DefectCount] >= config.multiple_bubble_min_defects).Label];
end
TwoOverlapMask = ismember(L, twoOverlapLabels);
MultipleOverlapMask = ismember(L, multipleOverlapLabels);
if config.verbose && ~isempty(overlapSeparationInfo)
    fprintf('Overlap separation decisions:\n');
    for infoIndex = 1:numel(overlapSeparationInfo)
        fprintf('  Label %d: %d defects -> %s\n', ...
            overlapSeparationInfo(infoIndex).Label, ...
            overlapSeparationInfo(infoIndex).DefectCount, ...
            overlapSeparationInfo(infoIndex).Method);
    end
    fprintf('Two-overlap source clusters: %d\n', numel(twoOverlapLabels));
    fprintf('Multiple-overlap source clusters: %d\n', numel(multipleOverlapLabels));
end

% --- Stage 5b: SINGLE-bubble branch (boundary reconstruction) ---------------
% Remove the overlapped objects from the full boundary set, then rebuild a
% label image for the remaining single (non-overlapped) bubbles.
struct_main=RemovalObjectFromMainBoundary(B,L,RealOverlap);
LabelledImage=ConstructLabelImageFromBoundary(BW,struct_main);
% Defensive guard: drop any single-bubble region that is implausibly large
% (e.g. a wall/border that survived as one giant "bubble"). The overlapped
% branch already filters big regions; this protects the single-bubble branch.
singleStats = regionprops(LabelledImage, 'Area');
if ~isempty(singleStats)
    oversized = find([singleStats.Area] > config.single_bubble_max_area);
    for lbl = oversized
        LabelledImage(LabelledImage == lbl) = 0;
    end
    LabelledImage = ReorderingLabels(LabelledImage);
end
Labelseg=LabelledImage;

% --- Stage 6: Combine both branches into one label image --------------------
% Combine watershed-segmented overlapped bubbles with the single bubbles.
% Place single-bubble labels (offset to stay unique) ONLY where the segmented
% image is empty, so spatially overlapping pixels cannot corrupt labels and
% create a spurious oversized region.
offset = max(FinalSegmentedObjects(:));
L_combined = FinalSegmentedObjects;
singleMask = (LabelledImage > 0) & (FinalSegmentedObjects == 0);
L_combined(singleMask) = LabelledImage(singleMask) + offset;

% --- Safety net: never lose a bubble ---------------------------------------
% A bubble can disappear when it was removed from the single-bubble branch as
% "overlapped" but the watershed branch failed to reproduce it (e.g. the
% cluster was not split and then dropped as an over-large region). Recover any
% part of the cleaned mask BW that ended up with no label so that every bubble
% present in BW is represented in the final result.
recoverMask = BW & (L_combined == 0);
% Drop thin watershed ridge slivers (1-2 px separation lines) so they do not
% become spurious sliver "bubbles"; keep genuine blob-shaped regions.
recoverMask = imopen(recoverMask, strel('disk', 3));
recoverMask = bwareaopen(recoverMask, config.min_object_area);
recoveredCC = bwconncomp(recoverMask);
recoverOffset = max(L_combined(:));
for recoverIndex = 1:recoveredCC.NumObjects
    recoverOffset = recoverOffset + 1;
    L_combined(recoveredCC.PixelIdxList{recoverIndex}) = recoverOffset;
end
if config.verbose && recoveredCC.NumObjects > 0
    fprintf('Safety net recovered %d bubble(s) dropped by both branches.\n', ...
        recoveredCC.NumObjects);
end

% Compact labels to 1-..-N so there are no empty gaps (empty labels would show
% up as zero-area phantom objects with NaN diameter).
L_combined = ReorderingLabels(L_combined);
if config.show_figures
    figure('Name','L_combined');
    imshow(L_combined)
end

% --- Stage 7: Final measurements, diagnostics and reporting -----------------
Object_information_last=GetObjectProperties(L_combined,config);
if config.verbose
    % Optional diagnostic: lists per-object sizes and flags statistical outliers.
    sizeReport = InspectBubbleSizes(L_combined,config);
end
PostProcessing(FinalSegmentedObjects,Labelseg,I_Formask,Object_information_last,config,TwoOverlapMask,MultipleOverlapMask);

%% Close all figures (cleanup)
% if ~config.show_figures
%     close all;
% end
% 
