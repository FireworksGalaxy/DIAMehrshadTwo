function  LL_new = MergeLabels(LL, label_one, label_two)
 
    % Copy the original labeled image
    LL_new = LL;

    % Create a binary mask for the first label
    label1Pixels = (LL == label_one);

    % Enlarge the mask (dilate one layer)
    label1Pixels = imdilate(label1Pixels, true(2)); % Structuring element: 2x2 ones

    % Assign label1 pixels to label2
    LL_new(label1Pixels) = label_two;

    % Note: label re-ordering is performed by the caller (ReorderingLabels)
    % immediately after merging, so it is intentionally not repeated here.
end
