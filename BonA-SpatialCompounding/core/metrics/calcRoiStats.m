function [stats, keepMask] = calcRoiStats(map, roiMask, baRef, useMad)

    % Calculate ROI statistics with optional three-scaled-MAD rejection
    % Use a scalar or map for baRef, or [] for descriptive statistics only

    % Select finite estimates within the geometric ROI
    assert(isequal(size(map),size(roiMask)), 'The B/A map and ROI mask must share one grid.');
    validMask = roiMask & isfinite(map);
    vals = map(validMask);
    assert(~isempty(vals), 'The ROI contains no finite B/A estimates.');
    keep = true(size(vals));

    % Apply the same MAD selection to all regional statistics
    if useMad
        medVal = median(vals);
        madScale = 1.4826 * median(abs(vals - medVal));
        if madScale > 0
            keep = abs(vals - medVal) <= 3 * madScale;
        else
            keep = vals == medVal;
        end
    end
    keepMask = false(size(roiMask));
    keepMask(validMask) = keep;
    kept = vals(keep);

    % Calculate the regional mean, spatial spread, and retained counts
    stats.mean = mean(kept);
    stats.spatialStd = std(kept,0);
    stats.cvPercent = 100 * stats.spatialStd / abs(stats.mean);
    stats.roiPixelCount = nnz(roiMask);
    stats.roiCount = numel(vals);
    stats.invalidRoiCount = nnz(roiMask & ~isfinite(map));
    stats.inlierCount = nnz(keep);
    stats.outlierCount = numel(vals) - stats.inlierCount;
    stats.outlierFraction = stats.outlierCount / stats.roiCount;

    % Calculate percentage mean error and pixel errors when a target is supplied
    if ~isempty(baRef)
        if isscalar(baRef)
            target = baRef;
        else
            target = baRef(keepMask);
        end
        stats.meanTruth = mean(target);
        stats.meanErrPct = 100 * abs(stats.mean - stats.meanTruth) / abs(stats.meanTruth);
        stats.rmse = sqrt(mean((kept - target).^2));
        stats.mapePercent = 100 * mean(abs(kept - target) ./ abs(target));
    end

end
