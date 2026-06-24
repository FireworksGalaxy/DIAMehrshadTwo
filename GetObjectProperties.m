function [Object_information, B_Deformed] = GetObjectProperties(input, config, Bpre, Lpre)
% GetObjectProperties - Unified function to extract object properties
% Consolidates GetObjectInformation and GetObjectInformationfromLabelImage
%
% Usage:
%   Object_info = GetObjectProperties(BW)              % Binary image input (default config)
%   Object_info = GetObjectProperties(L, config)       % Label image input with custom config
%   Object_info = GetObjectProperties(BW, config, B, L)% Reuse precomputed boundaries/labels
%   [Object_info, B_Def] = GetObjectProperties(...)    % Also returns deformed objects
%
% Input:
%   input       - Either binary image (BW) or label image (L)
%   config      - (Optional) Configuration struct from BubbleDetectionConfig()
%   Bpre, Lpre  - (Optional) Precomputed boundaries/label image to avoid
%                 recomputing bwboundaries for the same input.
%
% Output:
%   Object_information - Struct array with properties (Area, BoundingBox, Eccentricity, Centroid, Perimeter, SauterDia, Roundness)
%   B_Deformed - Cell array of deformed object boundaries (optional, only if nargout > 1)

    % Use default config if not provided
    if nargin < 2
        config = BubbleDetectionConfig();
    end

    % Reuse precomputed boundaries/labels when supplied. Otherwise determine
    % the input type and extract boundaries and labels. Use the data type to
    % disambiguate: logical masks are labelled by connected components,
    % numeric label images are used directly. A numeric 0/1 mask (binary
    % stored as double) is also treated as a mask so it is labelled correctly
    % instead of collapsing into a single region.
    if nargin >= 4 && ~isempty(Bpre) && ~isempty(Lpre)
        B = Bpre;
        L = Lpre;
    elseif islogical(input)
        % Binary mask: label connected components
        [B, L] = bwboundaries(input, 'noholes');
    elseif max(input(:)) > 1
        % Numeric label image: use directly as the label matrix
        [B, ~] = bwboundaries(input, 'noholes');
        L = input;
    else
        % Numeric 0/1 mask: relabel connected components
        [B, L] = bwboundaries(logical(input), 'noholes');
    end
    
    % Extract region properties
    Object_information = regionprops(L, 'Area', 'BoundingBox', 'Eccentricity', 'Centroid', 'Perimeter');
    
    % Calculate derived properties (Sauter diameter and Roundness)
    counter = 0;
    B_Deformed = {};
    for i = 1:length(Object_information)
        A = Object_information(i).Area;
        P = Object_information(i).Perimeter;
        % Guard against degenerate regions (zero area/perimeter) so that
        % derived metrics never become NaN/Inf and poison downstream stats.
        if P > 0
            Object_information(i).SauterDia = (A * 4) / P;
        else
            Object_information(i).SauterDia = 0;
        end
        if A > 0
            Object_information(i).Roundness = (P^2) / (4 * pi * A);
        else
            Object_information(i).Roundness = 0;
        end
        
        % Collect deformed objects if second output requested (Roundness > threshold from config)
        if nargout > 1 && Object_information(i).Roundness > config.roundness_threshold_deformed
            counter = counter + 1;
            B_Deformed{counter, 1} = B{i, 1};
            B_Deformed{counter, 2} = i;
        end
    end
end
