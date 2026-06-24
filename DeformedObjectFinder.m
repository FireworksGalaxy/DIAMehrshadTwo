function DeformedObject=DeformedObjectFinder(BW,Object_information,config,B)
% DeformedObjectFinder - Identify non-round (deformed) objects, which are the
% candidates for being overlapped bubbles.
%
% Usage:
%   DeformedObject = DeformedObjectFinder(BW, Object_information, config, B)
%
% Input:
%   BW                 - binary mask (used only if B is not supplied).
%   Object_information - per-object struct array (must contain .Roundness),
%                        indexed the same way as the boundaries B.
%   config             - configuration struct (uses roundness_threshold_deformed).
%   B                  - (Optional) precomputed boundaries from bwboundaries.
%
% Output:
%   DeformedObject - Nx2 cell array: column 1 = boundary, column 2 = the
%                    object's index in the original label image. Empty when no
%                    object exceeds the roundness threshold.

% Use default config if not provided (backwards compatibility)
if nargin < 3 || isempty(config)
    config = BubbleDetectionConfig();
end

% Reuse precomputed boundaries when supplied; otherwise compute them.
if nargin < 4 || isempty(B)
    B = bwboundaries(BW,'noholes');
end

DeformedObject = {};   % Ensure defined even when no deformed objects exist
counter = 0;
for i = 1:numel(B)
    % In DeformedObjectFinder.m, change the condition to:
    if Object_information(i).Roundness > config.roundness_threshold_deformed 
        counter = counter + 1;
        DeformedObject{counter,1} = B{i,1};
        DeformedObject{counter,2} = i;
    end
end
end
