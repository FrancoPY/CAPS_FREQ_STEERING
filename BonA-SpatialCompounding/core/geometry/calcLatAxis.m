function xBA = calcLatAxis(xAxis, winSize, winStep)

    % Locate each lateral average at the mean coordinate of its input samples.
    nCols = floor((numel(xAxis) - winSize) / winStep) + 1;
    xBA = NaN(1, nCols);

    % Use the same windows as the MLR and DM envelope averages
    for iCol = 1:nCols
        colStart = 1 + (iCol - 1) * winStep;
        colEnd = colStart + winSize - 1;
        xBA(iCol) = mean(xAxis(colStart:colEnd));
    end
end
