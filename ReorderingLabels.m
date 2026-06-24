function LL=ReorderingLabels(LL)
% Relabel an image so labels are consecutive 1..N (background 0 stays 0).
% Vectorized remap (no per-label loop).
[uniqueLabels, ~, newIdx] = unique(LL);
remap = zeros(numel(uniqueLabels), 1);
nonzero = uniqueLabels ~= 0;
remap(nonzero) = 1:nnz(nonzero);
LL = reshape(remap(newIdx), size(LL));
end
