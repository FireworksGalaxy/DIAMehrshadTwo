function Extracted_object=GetExtractedObject(DeformedObject,BW,L)
% Build a label image of the deformed (candidate overlapped) objects using the
% original region labels. This is exact and faster than polygon filling.
if nargin < 3 || isempty(L)
    [~, L] = bwboundaries(BW,'noholes');
end
Extracted_object = zeros(size(L));  % Label matrix for the extracted objects

% Remap each original label to its sequential index k, then apply in a single
% vectorized lookup instead of scanning the image once per object.
remap = zeros(max(L(:)) + 1, 1);            % index = origLabel + 1 (label 0 stays 0)
for k = 1:size(DeformedObject,1)
    remap(DeformedObject{k,2} + 1) = k;     % Last writer wins (matches loop behaviour)
end
Extracted_object = reshape(remap(L + 1), size(L));
end
