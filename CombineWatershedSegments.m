 function LL=CombineWatershedSegments(watershedLabels,config)
% CombineWatershedSegments - Post-process the raw watershed output into clean
% bubble labels.
%
% Usage:
%   LL = CombineWatershedSegments(watershedLabels, config)
%
% Steps:
%   1. Convert the watershed result to a label image.
%   2. Remove over-large regions (walls / merged blobs) via RemoveBigRegionandWall.
%   3. Iteratively merge regions whose centroids are very close together
%      (fixes small watershed fragments without ellipse-guided fitting).
%
% Output:
%   LL - final label image of the separated overlapped bubbles.
%% Initialization

[BB, LL] = bwboundaries(watershedLabels,'noholes');

LL=RemoveBigRegionandWall(LL,config.area_threshold_big_regions,length(BB));
LL=ReorderingLabels(LL);
ShowLabels(LL,'It is the fourth Image',config);

%% Preparation before merging
Struct_ERM=EstimateRegionsToMerge(LL,config);
Struct_BFM=Struct_ERM(1).Struct_BFMInside;
i=0;
while(~isempty(Struct_BFM(1).obj_one) && i<config.merge_max_iterations)
LL = MergeLabels(LL, Struct_BFM(1).obj_one, Struct_BFM(1).obj_two);
LL=ReorderingLabels(LL);
iterationText = sprintf('Iteration %d', i);
ShowLabels(LL,iterationText,config);
Struct_ERM=EstimateRegionsToMerge(LL,config);
Struct_BFM=Struct_ERM(1).Struct_BFMInside;
i=i+1;    
end

end
