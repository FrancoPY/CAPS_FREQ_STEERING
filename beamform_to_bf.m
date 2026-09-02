baseDir = fullfile(pwd,'FRECUENCIA-STEERING');
probe = 'L14-5u';

freqs = [5e6, 7e6, 9e6];
steeringAngles = [-5, 0, 5, 10, 15];

nFrames = 6;
fNumber = 3;
c0_bf   = 1500;

% Lista de phantoms a procesar (nombres exactos de tus scripts de barrido)
phantomNames = { ...
    'Phantom_bgnd_BA6_CAPS_FreqSteering_v1', ...
    'Phantom_inc_BA11_CAPS_FreqSteering_v1', ...
};

for p = 1:numel(phantomNames)
    sample = phantomNames{p};

    for f0 = freqs
        freqStr = sprintf('%dMHz', round(f0/1e6));

        for steeringAngle = steeringAngles

            rawDir = fullfile(baseDir, sample, probe, freqStr, 'rf', sprintf('Angle_%d', steeringAngle));
            bfDir  = fullfile(baseDir, sample, probe, freqStr, 'bf', sprintf('Angle_%d', steeringAngle));
            if ~exist(bfDir, 'dir'); mkdir(bfDir); end

            if ~exist(rawDir, 'dir')
                fprintf('--- SKIP %s | %s | Angle %d (no existe rawDir) ---\n', sample, freqStr, steeringAngle);
                continue
            end

            fprintf('--- Procesando %s | %s | Angle %d ---\n', sample, freqStr, steeringAngle);

            for i = 1:nFrames

                candidates = dir(fullfile(rawDir, sprintf('%s_f%d_sub1_of_m*_*kPa.mat', sample, i)));
                if isempty(candidates)
                    error('No se encontraron archivos de subapertura para %s, %s, angulo %d, frame %d', ...
                        sample, freqStr, steeringAngle, i);
                end
                tok = regexp(candidates(1).name, 'sub1_of_m(\d+)_', 'tokens');
                m = str2double(tok{1}{1});

                % Cargar y sumar coherentemente las m subaperturas 
                rfSumPrebf = [];
                for j = 1:m
                    fileSub = dir(fullfile(rawDir, sprintf('%s_f%d_sub%d_of_m%d_*kPa.mat', sample, i, j, m)));
                    if isempty(fileSub)
                        error('Falta la subapertura %d/%d para %s, %s, angulo %d, frame %d', ...
                            j, m, sample, freqStr, steeringAngle, i);
                    end
                    sJ = load(fullfile(rawDir, fileSub(1).name));
                    if j == 1
                        rfSumPrebf = sJ.rfPrebf;
                        fs = sJ.fs;
                        offsetBase = sJ.offset_base;
                    else
                        rfSumPrebf = rfSumPrebf + sJ.rfPrebf;
                    end
                end

                % Cargar la transmisión de apertura completa 
                fileFull = dir(fullfile(rawDir, sprintf('%s_f%d_full_m%d_*kPa.mat', sample, i, m)));
                if isempty(fileFull)
                    error('Falta la apertura completa para %s, %s, angulo %d, frame %d', ...
                        sample, freqStr, steeringAngle, i);
                end
                sFull = load(fullfile(rawDir, fileFull(1).name));

                rfLp = bfPlaneWaveSimu(rfSumPrebf', fs, fNumber, steeringAngle, c0_bf, offsetBase);
                rfHp = bfPlaneWaveSimu(sFull.rfPrebf', fs, fNumber, steeringAngle, c0_bf, offsetBase);

                nSamples = size(rfLp, 1);
                nLines   = size(rfLp, 2);
                zAxis = (1:nSamples)' * (c0_bf / fs) / 2;
                xAxis = (0:nLines-1) * 0.3e-3;
                xAxis = xAxis - mean(xAxis);

                save(fullfile(bfDir, sprintf('%s_f%d_LP.mat', sample, i)), 'rfLp', 'xAxis', 'zAxis', 'fs', 'm');
                save(fullfile(bfDir, sprintf('%s_f%d_HP.mat', sample, i)), 'rfHp', 'xAxis', 'zAxis', 'fs', 'm');

                fprintf('  Frame %d (m=%d) beamformeado y guardado en bf/\n', i, m);
            end
        end
    end
end

fprintf('\nListo. Carpetas "bf" generadas para cada phantom, frecuencia y ángulo (CAPS).\n');