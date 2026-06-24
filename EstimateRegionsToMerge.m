function Struct_ERM=EstimateRegionsToMerge(LL,config)
% EstimateRegionsToMerge - Find pairs of labelled regions close enough to merge.
%
% Usage:
%   Struct_ERM = EstimateRegionsToMerge(LL, config)
%
% Computes every region centroid and all pairwise distances, then returns the
% pairs whose centroid distance is below config.merge_distance_threshold.
%
% Output:
%   Struct_ERM(1).Struct_BFMInside - struct array of candidate merge pairs with
%       fields obj_one, obj_two and distance (empty if no pairs qualify).
%   Struct_ERM(1).length           - number of candidate pairs.
stats = regionprops(LL, 'Centroid'); % Centroids of current labelled regions

% Extract all centroids and compute pairwise Euclidean distances using a
% Gram-matrix formulation (fast BLAS matmul, no large 3-D intermediate and no
% Statistics Toolbox required).
centroids = cat(1, stats.Centroid);
n = size(centroids, 1);

sq = sum(centroids.^2, 2);            % squared norm of each centroid
D2 = sq + sq' - 2 * (centroids * centroids');
D2(D2 < 0) = 0;                       % clamp tiny negatives from round-off
D = sqrt(D2);

% Find all region pairs within merge threshold (upper triangle only, no redundancy)
[i, j] = find(D < config.merge_distance_threshold & triu(ones(size(D)), 1));

% Extract distances for the found pairs
distances = D(sub2ind(size(D), i, j));

% Create struct array from vectorized results
if ~isempty(i)
    Struct_BFM = struct('obj_one', num2cell(i), 'obj_two', num2cell(j), 'distance', num2cell(distances));
else
    Struct_BFM = struct('obj_one', [], 'obj_two', [], 'distance', []); % Empty if no pairs found
end

% Assign to Struct_ERM
Struct_ERM(1).Struct_BFMInside=Struct_BFM;
Struct_ERM(1).length=length(Struct_BFM);
end
