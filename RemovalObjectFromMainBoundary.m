function struct_main=RemovalObjectFromMainBoundary(struct_main,L,RealOverlap)
% Remove boundaries whose region in the original label image L is (mostly)
% covered by the overlapped-object mask RealOverlap.
%
% struct_main : cell array of boundaries from [B,L]=bwboundaries(BW). Index i
%               corresponds to label i in L.
% L           : label image aligned with struct_main indices.
% RealOverlap : binary mask of detected overlapped objects.

overlapMask = RealOverlap > 0;

% Candidate labels: any label touching the overlap mask
candidateLabels = unique(L(overlapMask & L > 0))';

% Keep only labels that are mostly inside the overlap mask (>50% coverage),
% so neighbouring regions clipped by the mask are not removed by mistake.
labelsToRemove = [];
for lbl = candidateLabels
    region = (L == lbl);
    coverage = nnz(region & overlapMask) / nnz(region);
    if coverage > 0.5
        labelsToRemove(end+1) = lbl; %#ok<AGROW>
    end
end

% Remove in descending order so earlier indices stay valid
labelsToRemove = sort(labelsToRemove, 'descend');
for lbl = labelsToRemove
    if lbl >= 1 && lbl <= numel(struct_main)
        struct_main(lbl) = [];
    end
end
end
