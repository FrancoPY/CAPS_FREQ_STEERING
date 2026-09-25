function eff = buildPathEffectiveAcousticMaps(sampleDef, x, z, f0MHz)

    % buildPathEffectiveAcousticMaps  Local and cumulative acoustic maps
    %
    % Inputs
    %   sampleDef : struct with acoustic properties or profile definition
    %   x         : lateral axis [m]
    %   z         : axial axis [m]
    %   f0MHz     : central frequency [MHz]
    %
    % Output
    %   eff       : struct with local and cumulative maps

    x = x(:).';
    z = z(:);

    % Local maps
    [acMap, yMap, cMap, boaMap] = getLocalAcousticMaps(sampleDef, x, z);

    % Band attenuation
    alpha1Map = acMap .* f0MHz.^yMap * 100;
    alpha2Map = acMap .* (2 * f0MHz).^yMap * 100;

    % Output
    eff.local.ac = acMap;
    eff.local.alphaPower = yMap;
    eff.local.speedOfSound = cMap;
    eff.local.boa = boaMap;
    eff.local.alpha1 = alpha1Map;
    eff.local.alpha2 = alpha2Map;

    eff.cum.alpha1 = cumulativeMeanMap(alpha1Map, z);
    eff.cum.alpha2 = cumulativeMeanMap(alpha2Map, z);
    eff.cum.speedOfSound = cumulativeHarmonicMeanMap(cMap, z);
    eff.cum.boa = cumulativeMeanMap(boaMap, z);

end

function [acMap, yMap, cMap, boaMap] = getLocalAcousticMaps(sampleDef, x, z)

    % getLocalAcousticMaps  Build local property maps on the BA grid

    Nx = numel(x);
    Nz = numel(z);
    [X, Z] = meshgrid(x, z);

    % Homogeneous fallback
    acMap = sampleDef.ac * ones(Nz, Nx);
    yMap = sampleDef.alphaPower * ones(Nz, Nx);
    cMap = sampleDef.speedOfSound * ones(Nz, Nx);
    boaMap = sampleDef.boa * ones(Nz, Nx);

    % Full maps
    if isfield(sampleDef, 'x') && isfield(sampleDef, 'z') && ...
            isfield(sampleDef, 'acMap') && isfield(sampleDef, 'alphaPowerMap') && ...
            isfield(sampleDef, 'speedOfSoundMap') && isfield(sampleDef, 'boaMap')

        acMap = interp2(sampleDef.x(:).', sampleDef.z(:), sampleDef.acMap, X, Z, 'nearest');
        yMap = interp2(sampleDef.x(:).', sampleDef.z(:), sampleDef.alphaPowerMap, X, Z, 'nearest');
        cMap = interp2(sampleDef.x(:).', sampleDef.z(:), sampleDef.speedOfSoundMap, X, Z, 'nearest');
        boaMap = interp2(sampleDef.x(:).', sampleDef.z(:), sampleDef.boaMap, X, Z, 'nearest');

        acMap(isnan(acMap)) = sampleDef.ac;
        yMap(isnan(yMap)) = sampleDef.alphaPower;
        cMap(isnan(cMap)) = sampleDef.speedOfSound;
        boaMap(isnan(boaMap)) = sampleDef.boa;
        return
    end

    % Profile type
    if ~isfield(sampleDef, 'profileType')
        return
    end

    switch lower(sampleDef.profileType)

        case 'layered'
            if ~isfield(sampleDef, 'zBreaks') || ...
                    ~isfield(sampleDef, 'acLayers') || ...
                    ~isfield(sampleDef, 'alphaPowerLayers') || ...
                    ~isfield(sampleDef, 'speedOfSoundLayers') || ...
                    ~isfield(sampleDef, 'boaLayers')
                return
            end

            zEdges = [0, sampleDef.zBreaks(:).', inf];
            nLayers = numel(zEdges) - 1;

            acMap = zeros(Nz, Nx);
            yMap = zeros(Nz, Nx);
            cMap = zeros(Nz, Nx);
            boaMap = zeros(Nz, Nx);

            for iLayer = 1:nLayers
                mask = Z >= zEdges(iLayer) & Z < zEdges(iLayer + 1);
                acMap(mask) = sampleDef.acLayers(iLayer);
                yMap(mask) = sampleDef.alphaPowerLayers(iLayer);
                cMap(mask) = sampleDef.speedOfSoundLayers(iLayer);
                boaMap(mask) = sampleDef.boaLayers(iLayer);
            end

        case 'inclusion'
            if ~isfield(sampleDef, 'acInc') || ~isfield(sampleDef, 'acBgnd') || ...
                    ~isfield(sampleDef, 'alphaPowerInc') || ~isfield(sampleDef, 'alphaPowerBgnd') || ...
                    ~isfield(sampleDef, 'speedOfSoundInc') || ~isfield(sampleDef, 'speedOfSoundBgnd') || ...
                    ~isfield(sampleDef, 'boaInc') || ~isfield(sampleDef, 'boaBgnd')
                return
            end

            if isfield(sampleDef, 'maskFcn')
                incMask = logical(sampleDef.maskFcn(X, Z));
            elseif isfield(sampleDef, 'inclusionCenter') && isfield(sampleDef, 'inclusionRadius')
                x0 = sampleDef.inclusionCenter(1);
                z0 = sampleDef.inclusionCenter(2);
                r = sampleDef.inclusionRadius;
                incMask = (X - x0).^2 + (Z - z0).^2 <= r^2;
            elseif isfield(sampleDef, 'xLimInc') && isfield(sampleDef, 'zLimInc')
                incMask = X >= sampleDef.xLimInc(1) & X <= sampleDef.xLimInc(2) & ...
                          Z >= sampleDef.zLimInc(1) & Z <= sampleDef.zLimInc(2);
            else
                return
            end

            acMap = sampleDef.acBgnd * ones(Nz, Nx);
            yMap = sampleDef.alphaPowerBgnd * ones(Nz, Nx);
            cMap = sampleDef.speedOfSoundBgnd * ones(Nz, Nx);
            boaMap = sampleDef.boaBgnd * ones(Nz, Nx);

            acMap(incMask) = sampleDef.acInc;
            yMap(incMask) = sampleDef.alphaPowerInc;
            cMap(incMask) = sampleDef.speedOfSoundInc;
            boaMap(incMask) = sampleDef.boaInc;

    end

end

function meanMap = cumulativeMeanMap(map, z)

    % cumulativeMeanMap  Cumulative mean along depth

    z = z(:);

    if z(1) > 0
        zPad = [0; z];
        mapPad = [map(1, :); map];
    else
        zPad = z;
        mapPad = map;
    end

    intMap = cumtrapz(zPad, mapPad, 1);
    meanPad = zeros(size(mapPad));

    idx = zPad > 0;
    meanPad(idx, :) = intMap(idx, :) ./ zPad(idx);
    meanPad(~idx, :) = mapPad(~idx, :);

    if z(1) > 0
        meanMap = meanPad(2:end, :);
    else
        meanMap = meanPad;
    end

end

function hMap = cumulativeHarmonicMeanMap(map, z)

    % cumulativeHarmonicMeanMap  Cumulative harmonic mean along depth

    invMean = cumulativeMeanMap(1 ./ map, z);
    hMap = 1 ./ invMean;

end