function [fileBA, prep] = computeEDMforHet(media, param, saveFlag, roi, frameRange, plotFlag)

    % computeEDM  Compute B/A map with the extended depletion model
    %
    % Inputs
    %   media      : struct with input/output paths
    %   param      : struct with acoustic and processing parameters
    %   saveFlag   : true to save .mat and figure outputs
    %   roi        : ROI [xMin xMax zMin zMax] in meters
    %   frameRange : frames to include in the preprocessing
    %   plotFlag   : true to compute B-mode previews
    %
    % Outputs
    %   fileBA     : output file path for saved B/A results
    %   prep       : struct with preprocessed data used in the computation

    if nargin < 6
        plotFlag = false;
    end

    % Prep
    needBmode = saveFlag || plotFlag;
    prep = prepareBAdata(media, param, needBmode, frameRange);

    % Params
    f0 = param.centralFreqMHz;
    probe = param.probe;
    freqStr = param.freqStr;
    sample = param.sample;
    reference = param.reference;

    % Axes
    z = prep.zAxis(:);
    xBA = prep.xBA(:).';
    zMat = z * ones(1, numel(xBA));
    f1 = prep.f1;

    % Acoustic defs
    if isfield(param, 'sampleDef')
        samDef = param.sampleDef;
    else
        samDef.ac = param.acSample;
        samDef.alphaPower = param.alphaPowerSample;
        samDef.speedOfSound = param.speedOfSoundSam;
        samDef.boa = 0;
    end

    if isfield(param, 'referenceDef')
        refDef = param.referenceDef;
    else
        refDef.ac = param.acReference;
        refDef.alphaPower = param.alphaPowerReference;
        refDef.speedOfSound = param.speedOfSoundRef;
        refDef.boa = param.boaReference;
    end

    % Cumulative effective maps
    samEff = buildPathEffectiveAcousticMaps(samDef, xBA, z, f0);
    refEff = buildPathEffectiveAcousticMaps(refDef, xBA, z, f0);

    alpha1S = samEff.cum.alpha1;
    alpha2S = samEff.cum.alpha2;
    alpha1R = refEff.cum.alpha1;
    alpha2R = refEff.cum.alpha2;
    cS = samEff.cum.speedOfSound;
    cR = refEff.cum.speedOfSound;
    boaR = refEff.cum.boa;

    % Model terms
    f2 = sqrt(edmKernel(alpha1R, alpha2R, zMat) ./ edmKernel(alpha1S, alpha2S, zMat));
    f3 = alpha1S ./ alpha1R;
    f4 = cS ./ cR;

    BAmap = (f1 .* f2 .* f3 .* f4.^3 .* (1 + 0.5 * boaR) - 1) * 2;

    % Output
    BAeff.image = BAmap;
    BAeff.lateral = xBA;
    BAeff.axial = z;

    outputDir = fullfile(media.resDir, 'edm');
    figDir = fullfile(media.figDir, 'mapsEDM', sprintf('%s_%s_%s', sample, probe, freqStr));

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end

    fileBA = fullfile(outputDir, sprintf('%s_ref_%s', sample, reference));

    % Save
    if saveFlag
        BmodeSAM = prep.BmodeSAM;
        BmodeREF = prep.BmodeREF;
        save(fileBA, 'BAeff', 'BmodeSAM', 'BmodeREF');

        saveName = sprintf('mapBA_%s_ref_%s_sam', sample, reference);
        figTitle = sprintf('%s - %s', sample, reference);

        plotBAmap(BAmap, xBA, z, [], [], [5, 12], 'turbo', figTitle, roi);
        exportgraphics(gcf, fullfile(figDir, [saveName '.png']), 'Resolution', 200);
        close(gcf);
    end

end

function K = edmKernel(alpha1, alpha2, zMat)

    % edmKernel  EDM attenuation kernel using alpha at f0 and 2f0

    tol = 1e-12;
    K = zeros(size(alpha1));

    idxY1 = abs(alpha2 - 2 * alpha1) < tol;

    if any(~idxY1, 'all')
        a1 = alpha1(~idxY1);
        a2 = alpha2(~idxY1);
        z = zMat(~idxY1);

        K(~idxY1) = a1 ./ (4 * (a2 - 2 * a1)) .* ...
            (0.5 * (1 - exp(-2 * a1 .* z)) - (a1 ./ a2) .* (1 - exp(-a2 .* z)));
    end

    if any(idxY1, 'all')
        a1 = alpha1(idxY1);
        z = zMat(idxY1);

        K(idxY1) = (1 / 16) * (1 - exp(-2 * a1 .* z) - 2 * a1 .* z .* exp(-2 * a1 .* z));
    end

    K(K <= 0) = eps;

end