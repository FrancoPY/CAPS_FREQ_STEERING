function prep = prepareBAdata(media, param, bmodeFlag, frameRange)

    % prepareBAdata  Load, filter, and preprocess RF data for B/A estimation
    %
    % Inputs
    %   media      : struct with input/output paths
    %   param      : struct with acoustic and processing parameters
    %   bmodeFlag  : true to compute B-mode previews
    %   frameRange : frames to include in the preprocessing
    %
    % Output
    %   prep       : struct with preprocessed data for EDM and DM

    % Params
    sample = param.sample;
    reference = param.reference;
    noiseLevel = param.noiseLevel;
    scaleFactor = param.pressureFactor;
    f0 = param.centralFreqMHz;
    c0 = param.speedOfSoundBf;
    order = param.filterOrder;

    % Load sample RF
    [lowSamFull, highSamFull, xAxis, zAxis, fs] = loadRfStack(media.samPath, sample, frameRange, noiseLevel);

    % Load reference RF
    [lowRefFull, highRefFull, ~, ~, fsRef] = loadRfStack(media.refPath, reference, frameRange, noiseLevel);

    % Filter
    fL = f0 - 1;
    fH = f0 + 1;
    nyq = fs / 2;
    hFilter = fir1(order, [fL * 1e6 / nyq, fH * 1e6 / nyq]);
    hFilter = cast(hFilter, 'like', lowSamFull);
    hKernel = reshape(hFilter, [], 1, 1);

    lowSamFull = convn(lowSamFull, hKernel, 'same');
    highSamFull = convn(highSamFull, hKernel, 'same');

    % Match the reference filter to its own sampling frequency
    hFilter = fir1(order, [fL * 1e6, fH * 1e6] / (fsRef / 2));
    hFilter = cast(hFilter, 'like', lowRefFull);
    hKernel = reshape(hFilter, [], 1, 1);
    lowRefFull = convn(lowRefFull, hKernel, 'same');
    highRefFull = convn(highRefFull, hKernel, 'same');

    % B-mode
    if bmodeFlag
        BmodeSAM = plotBmode(mean(lowSamFull, 3), xAxis, zAxis, false, [0.5e-2, 3.5e-2], false, 'Sample B-mode Image');
        BmodeREF = plotBmode(mean(lowRefFull, 3), xAxis, zAxis, false, [0.5e-2, 3.5e-2], false, 'Reference B-mode Image');
    else
        BmodeSAM = [];
        BmodeREF = [];
    end

    % Common envelopes
    envLowSamBase = mean(abs(hilbert(lowSamFull)), 3);
    envHighSamBase = mean(abs(hilbert(highSamFull)), 3);

    envLowRef = squeeze(mean(mean(abs(hilbert(lowRefFull)), 3), 2));
    envHighRef = squeeze(mean(mean(abs(hilbert(highRefFull)), 3), 2));

    envLowRef = envLowRef(:);
    envHighRef = envHighRef(:);

    % Lateral windows
    colLast = param.lastColIndex;
    overlap = param.overlap;
    width = param.smoothingWidth;
    colList = colLast:overlap:size(envLowSamBase, 2);
    xBA = xAxis(colList - overlap);

    envLowSamCols = movmean(envLowSamBase, [width - 1, 0], 2);
    envHighSamCols = movmean(envHighSamBase, [width - 1, 0], 2);

    envLowSamCols = envLowSamCols(:, colList);
    envHighSamCols = envHighSamCols(:, colList);

    % Common ratio term
    dz = zAxis(2) - zAxis(1);
    L = max(1, round(param.movingAvgWindow * (c0 / (f0 * 1e6)) / dz));

    dSam = scaleFactor * envLowSamCols - envHighSamCols;
    dRef = scaleFactor * envLowRef - envHighRef;

    pSam = sqrt(abs(dSam) .* envLowRef);
    pRef = sqrt(abs(dRef) .* envLowSamCols);

    pSamSmooth = movmean(pSam, L, 1);
    pRefSmooth = movmean(pRef, L, 1);

    % Output
    prep.xAxis = xAxis;
    prep.zAxis = zAxis;
    prep.xBA = xBA;
    prep.f1 = pSamSmooth ./ pRefSmooth;
    prep.BmodeSAM = BmodeSAM;
    prep.BmodeREF = BmodeREF;

end
