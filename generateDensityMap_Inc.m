sample = 'Phantom_inc_BA11_FrecSteering_v1';

Nx = 2000;
Ny = 2000;
dx = 0.02e-3;

sd_bg = 0.02;          % desviación estándar de fondo
scatterBoost = 3;      % factor multiplicador dentro de la inclusión
sd_inc = sd_bg * scatterBoost;

nFrames = 6;

% Geometría de la inclusión (debe coincidir con simulateInc9_copia.m)
centerDepth = 22.5e-3;  % [m]
radius = 7.5e-3;        % [m]
inclusion = makeDisc(Nx, Ny, round(centerDepth/dx), Ny/2, round(radius/dx));

outputDir = fullfile(pwd, 'densityMaps', sample);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Frames
for ii = 1:nFrames

    outName = fullfile(outputDir, sprintf('frame%d', ii));
    if exist([outName '.mat'], 'file')
        fprintf('[frame %d] skipped\n', ii);
        continue
    end

    % Density de fondo
    rng('shuffle');

    density = 1000 * ones(Nx, Ny);
    density = density + 1000 * sd_bg * randn(size(density));

    % Refuerzo de dispersores dentro de la inclusión
    extraNoise = 1000 * sd_inc * randn(size(density));
    density(inclusion > 0) = 1000 + extraNoise(inclusion > 0);

    density = single(density);

    % Save
    save(outName, 'density', '-v7.3');
    fprintf('[frame %d] saved\n', ii);

end
