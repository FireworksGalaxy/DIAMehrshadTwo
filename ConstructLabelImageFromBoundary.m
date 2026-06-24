function L_reconstructed=ConstructLabelImageFromBoundary(BW,struct_main)
% Reconstruct a label image for the given boundaries using the actual
% connected-component regions (exact and faster than polygon filling).
Llabel = bwlabel(BW);

% Build a remap from each connected-component label to its sequential index k,
% then apply it in a single vectorized lookup instead of scanning the whole
% image once per boundary (O(N*pixels) -> O(pixels)).
remap = zeros(max(Llabel(:)) + 1, 1);   % index = oldLabel + 1 (so label 0 stays 0)
for k = 1:numel(struct_main)
    boundary = struct_main{k};                  % k-th boundary ([row col] pairs)
    lbl = Llabel(boundary(1,1), boundary(1,2)); % Label at a boundary pixel
    if lbl > 0
        remap(lbl + 1) = k;                     % Last writer wins (matches loop behaviour)
    end
end
L_reconstructed = reshape(remap(Llabel + 1), size(Llabel));
end
