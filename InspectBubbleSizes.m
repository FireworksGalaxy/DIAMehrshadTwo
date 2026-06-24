function outlierTable = InspectBubbleSizes(L_combined, config)
% InspectBubbleSizes - Report per-object size metrics and flag outliers.
%
% Helps locate erroneous data points (for example a single spuriously large
% "size 70" bubble) by listing every object's geometry and highlighting the
% ones whose diameter is statistically far from the rest.
%
% Usage:
%   sizeReport = InspectBubbleSizes(L_combined, config);
%
% Input:
%   L_combined - final label image of detected bubbles
%   config     - configuration struct from BubbleDetectionConfig()
%
% Output:
%   outlierTable - table with Label, Area, Perimeter, Roundness, SauterDia,
%                  Diameter, RobustZ and Outlier flag for every object.

    if nargin < 2
        config = BubbleDetectionConfig();
    end

    info = GetObjectProperties(L_combined, config);
    n = numel(info);
    if n == 0
        outlierTable = table();
        if config.verbose
            disp('InspectBubbleSizes: no objects found.');
        end
        return;
    end

    idx       = (1:n)';
    area      = [info.Area]';
    perimeter = [info.Perimeter]';
    roundness = [info.Roundness]';
    sauter    = [info.SauterDia]';
    diameter  = sauter / config.calibration_divisor;

    % Exclude degenerate/empty labels (zero area) from the statistics. These
    % arise from gaps in the label numbering and are not real objects.
    valid = area > 0 & isfinite(diameter);
    nDegenerate = nnz(~valid);

    % Robust outlier detection on diameter using the median absolute deviation
    % (MAD). This is insensitive to the very outliers we are trying to find.
    med  = median(diameter(valid));
    madv = median(abs(diameter(valid) - med));
    if isnan(madv) || madv == 0
        madv = eps;
    end
    robustZ = 0.6745 * (diameter - med) ./ madv;
    isOutlier = valid & (abs(robustZ) > config.size_outlier_mad_threshold);

    outlierTable = table(idx, area, perimeter, roundness, sauter, diameter, ...
        robustZ, isOutlier, 'VariableNames', ...
        {'Label','Area','Perimeter','Roundness','SauterDia','Diameter','RobustZ','Outlier'});

    if config.verbose
        fprintf('\nInspectBubbleSizes: %d labels (%d real, %d empty), median diameter = %.4g\n', ...
            n, nnz(valid), nDegenerate, med);
        flagged = outlierTable(isOutlier, :);
        if isempty(flagged)
            disp('No size outliers detected.');
        else
            fprintf('Flagged %d potential outlier(s) (largest first):\n', height(flagged));
            disp(sortrows(flagged, 'Diameter', 'descend'));
        end
    end
end
