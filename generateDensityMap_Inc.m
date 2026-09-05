sample = 'Phantom_bgnd_BA6_CAPS_FreqSteering_v1';

Nx = 2000;
Ny = 2000;

sd = 0.02;
nFrames = 1;

outputDir = fullfile(pwd, 'densityMaps', sample);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Frames
for ii = 1:nFrames

    outName = fullfile(outputDir, sprintf('frame%d', ii));
    if exist([outName '.mat'], 'file')
        fprintf('[frame %d] skipped\n', ii);
        continue
    end

    % Density
    rng('shuffle');

    density = 1000 * ones(Nx, Ny);
    density = density + 1000 * sd * randn(size(density));
    density = single(density);

    % Save
    save(outName, 'density', '-v7.3');
    fprintf('[frame %d] saved\n', ii);

end
