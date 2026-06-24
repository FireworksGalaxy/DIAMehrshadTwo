function config = BubbleDetectionConfig()
% BubbleDetectionConfig - Centralized configuration for bubble detection algorithm
% Returns a struct with all tunable parameters
%
% Usage:
%   config = BubbleDetectionConfig();
%   sharpen_radius = config.sharpen_radius_main;
%
% All parameters can be easily modified here for different image resolutions/bubble sizes

    %% IMAGE PRE-PROCESSING PARAMETERS
    
    % Initial image sharpening.
    % NOTE: radius 200 / amount 10 is very aggressive and can create halos
    % around bubble edges. If boundaries look distorted, consider tuning
    % (typical radius ~2-50, amount ~0.5-3) against your own images.
    config.sharpen_radius_main = 200;        % Radius for initial sharpening (higher = more effect)
    config.sharpen_amount_main = 15;         % Amount intensity for sharpening
    
    % Adaptive binarization for main image
    config.binarize_sensitivity_main = 0.5;  % Sensitivity for foreground (0-1, higher = stricter)
    
    % Morphological filtering
    config.morphology_disk_size = 15;         % Disk structuring element radius for open operation
    
    % Minimum object area (noise removal). Objects smaller than this (px^2)
    % are removed before analysis to suppress speckle noise.
    config.min_object_area = 50;

    % Boundary bubble removal. When enabled, any detected object touching the
    % outer boundary band is removed before bubble analysis. Margin 1 removes
    % bubbles cut by the four image edges; increase it to remove bubbles near
    % the edges as well.
    config.remove_boundary_bubbles = true;
    config.boundary_margin_pixels = 1;
    
    %% OVERLAPPED REGION DETECTION PARAMETERS
    
    % Threshold for identifying overlapped bubbles
    config.area_threshold_overlapped = 10000; % Pixels² - regions larger than this are overlapped

    % Overlap detection method used before watershed. Options:
    % 'curvature'      - objects with sharp inward V-shaped boundary points
    % 'roundness_area' - original roundness + area heuristic
    % 'combined'       - either method can mark an object as overlapped
    config.overlap_detection_method = 'curvature';
    
    % Deformed object classification (roundness-based)
    config.roundness_threshold_deformed = 1.3; % Objects with roundness > 1.3 are classified as deformed

    % Experimental convexity-defect overlap detection. This can be slow on
    % large images, so it is disabled by default.
    config.convexity_defect_depth_threshold = 15; % Pixels
    config.show_convexity_defects = false;

    % Faster overlap detection based on sharp concave boundary curvature. A
    % V-shaped inward turn is treated as a candidate overlap/neck point.
    config.show_curvature_overlaps = true;
    config.curvature_window_size = 12;              % Boundary points on each side
    config.curvature_angle_threshold_degrees = 120; % Max V-angle; lower is sharper
    config.curvature_min_depth_pixels = 4;          % Local V-depth above chord
    config.curvature_min_point_separation = 20;     % Avoid duplicate nearby points
    config.curvature_overlay_point_radius = 6;

    % Separation strategy after overlapped objects are detected. Objects with
    % exactly two defect/neck points are treated as two-bubble overlaps and are
    % separated by watershed. Objects with more defect points can use either
    % watershed or an experimental implicit cubic curve split.
    config.separation_defect_count_source = 'convexity'; % Options: 'curvature' or 'convexity'
    config.two_bubble_defect_count = 2;
    config.multiple_bubble_min_defects = 3;
    config.multiple_overlap_separation_method = 'watershed'; % Options: 'watershed' or 'curve_fit'
    config.curve_fit_fallback_to_watershed = true;
    config.curve_fit_zero_tolerance = 0.04;
    config.curve_fit_max_boundary_points = 600;
    config.curve_fit_min_component_area = 50;
    
    %% WATERSHED SEGMENTATION PARAMETERS (For overlapped bubbles)
    
    % Sharpening before watershed
    config.sharpen_radius_watershed = 50;      % Radius for watershed pre-sharpening
    config.sharpen_amount_watershed = 10;      % Amount intensity for watershed sharpening
    
    % Binarization for watershed
    config.binarize_sensitivity_watershed = 0.5; % Sensitivity for watershed foreground
    
    % Extended minima for watershed markers
    config.extended_minima_value = 0.25;        % Threshold for creating watershed markers (sensitivity)
    
    % Adaptive watershed markers (optional). When true, the extended-minima
    % depth is scaled to a fraction of the per-cluster max distance transform
    % instead of the fixed extended_minima_value above. Default false keeps
    % the original behaviour.
    config.watershed_adaptive_markers = false;
    config.extended_minima_fraction = 0.25;     % Used only when adaptive markers are enabled
    
    %% POST-PROCESSING PARAMETERS
    
    % Removing large regions (likely walls/borders)
    config.area_threshold_big_regions = 8000;  % Pixels² - remove regions > this area

    % Defensive cap for the single-bubble branch: any non-overlapped region
    % larger than this (px^2) is treated as a wall/border artifact and dropped
    % so it cannot leak through as one giant spurious bubble. Set generously so
    % only clearly non-physical regions are removed; lower it if walls persist.
    config.single_bubble_max_area = 50000;
    
    % Calibration: pixels per mm, used to convert Sauter diameter from pixels
    % to mm. You can measure this with CalibrateImageFromRuler().
    config.calibration_divisor = 15;
    
    % Merging watershed regions that belong to the same bubble
    config.watershed_merge_strategy = 'ellipse'; % Options: 'ellipse' or 'centroid'
    config.merge_distance_threshold = 25;      % Pixels - used only by centroid merge strategy
    config.merge_max_iterations = 10;          % Max iterations for region merging

    % Ellipse-guided watershed merging. Neighboring watershed regions are
    % merged only when the merged boundary is a better ellipse fit than the
    % two separate boundaries.
    config.ellipse_merge_error_ratio = 0.90;       % merged RMSE must be below this fraction of separate RMSE
    config.ellipse_merge_min_improvement = 0.02;   % minimum absolute RMSE improvement
    config.ellipse_merge_max_merged_error = 0.35;  % reject merges with poor absolute ellipse fit
    config.ellipse_merge_neighbor_dilation_radius = 2; % pixels for finding neighbors across watershed ridges
    config.ellipse_merge_max_candidate_pairs = 150; % Max neighboring pairs evaluated per iteration
    config.ellipse_merge_max_boundary_points = 300; % Max boundary points used per ellipse fit
    config.show_ellipse_merge_figures = false;     % Show candidate pairs and accepted ellipse fits
    config.ellipse_merge_max_visual_pairs = 200;   % Limit candidate-pair lines drawn in diagnostics
    config.ellipse_merge_max_fit_figures = 8;      % Limit before/after ellipse fit diagnostic figures
    
    % Visualization filtering
    config.roundness_threshold_display = 1.8;  % Display only objects with roundness < 1.8
    
    % Size-outlier flagging (diagnostics). Objects whose diameter has a robust
    % z-score (MAD-based) above this are flagged by InspectBubbleSizes.
    config.size_outlier_mad_threshold = 3.5;
    
    %% FIGURE MANAGEMENT PARAMETERS
    
    % Control visualization output
    config.show_figures = false;                % Set to false to disable expensive/intermediate figure output
    config.show_object_separation_overlay = true; % Show final object-separation overlay figure
    config.show_postprocessing_summary_plots = true; % Show final Sauter-diameter and BSD figures
    config.verbose = true;                     % Set to false to suppress console output

    % Object-separation overlay appearance
    config.object_separation_line_width = 1;    % Boundary thickness in final overlay
    config.show_object_separation_numbers = false; % Show bubble index numbers in final overlay
    config.single_bubble_display_color = 'blue'; % Color for non-overlapped single bubbles
    config.two_overlap_display_color = 'red'; % Color for segmented objects from two-bubble overlap clusters
    config.multiple_overlap_display_color = 'green'; % Color for segmented objects from multi-overlap clusters
    
    %% DERIVED PARAMETERS (Auto-calculated, modify carefully)
    
    % These can be used to scale parameters based on image resolution
    % Leave empty or 0 to disable adaptive scaling
    config.reference_resolution = 0;           % Reference image height (0 = disabled)
    config.scale_factor = 1.0;                 % Current scale factor relative to reference
    
end
