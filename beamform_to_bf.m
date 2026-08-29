%% beamform_to_bf.m (versión CAPS)
% Suma coherentemente las m señales RF de subapertura de cada frame
% (reconstruyendo la señal equivalente a "apertura completa a presión
% m veces menor", ecuación 10 del paper de Avilés et al.), la usa como
% señal de baja presión (LP), y beamformea también la transmisión de
% apertura completa como señal de alta presión (HP).
%
% m se lee directamente de los archivos .mat generados por
% simulateBgnd6_copia.m / simulateInc9_copia.m (versión CAPS), así que
% no hace falta declararlo aquí: cada phantom puede haberse simulado con
% un m distinto.

baseDir = 'C:\Users\FRANCO PERALTA\Documents\MATLAB\LIM\CAPS\TAREA_CAPS';
probe   = 'L14-5u';
freqStr = '7MHz';

nFrames = 6;
fNumber = 3;
steeringAngle = 15;
c0_bf = 1500;   % velocidad de sonido

% Lista de phantoms a procesar
phantomNames = { ...
    'Phantom_homo_BA6_Steering15_CAPS', ...
    'Phantom_inc_BA11_Steering15_CAPS', ...
};

for p = 1:numel(phantomNames)

    sample = phantomNames{p};

    rawDir = fullfile(baseDir, sample, probe, freqStr, 'rf');
    bfDir  = fullfile(baseDir, sample, probe, freqStr, 'bf');
    if ~exist(bfDir, 'dir'); mkdir(bfDir); end

    fprintf('--- Procesando %s ---\n', sample);

    for i = 1:nFrames

        % Detectar m a partir de los archivos de subapertura existentes
        % (busca sub1_of_mK, K = 2,4,8,...,numElem)
        mDetected = [];
        candidates = dir(fullfile(rawDir, sprintf('%s_f%d_sub1_of_m*_*kPa.mat', sample, i)));
        if isempty(candidates)
            error('No se encontraron archivos de subapertura para %s, frame %d', sample, i);
        end
        tok = regexp(candidates(1).name, 'sub1_of_m(\d+)_', 'tokens');
        mDetected = str2double(tok{1}{1});
        m = mDetected;

        % --- Cargar y sumar coherentemente las m subaperturas ---
        rfSumPrebf = [];
        for j = 1:m
            fileSub = dir(fullfile(rawDir, sprintf('%s_f%d_sub%d_of_m%d_*kPa.mat', sample, i, j, m)));
            if isempty(fileSub)
                error('Falta la subapertura %d/%d para %s, frame %d', j, m, sample, i);
            end
            sJ = load(fullfile(rawDir, fileSub(1).name));

            if j == 1
                rfSumPrebf = sJ.rfPrebf;
                fs = sJ.fs;
                c0 = sJ.c0;
                offsetBase = sJ.offset_base;
            else
                rfSumPrebf = rfSumPrebf + sJ.rfPrebf;
            end
        end

        % --- Cargar la transmisión de apertura completa (alta presión) ---
        fileFull = dir(fullfile(rawDir, sprintf('%s_f%d_full_m%d_*kPa.mat', sample, i, m)));
        if isempty(fileFull)
            error('Falta la transmisión de apertura completa para %s, frame %d', sample, i);
        end
        sFull = load(fullfile(rawDir, fileFull(1).name));

        % offset_base debe ser el mismo para subaperturas y apertura
        % completa (mismo steering, misma física), se usa uno solo
        rfLp = bfPlaneWaveSimu(rfSumPrebf', fs, fNumber, steeringAngle, c0_bf, offsetBase);
        rfHp = bfPlaneWaveSimu(sFull.rfPrebf', fs, fNumber, steeringAngle, c0_bf, offsetBase);

        nSamples = size(rfLp, 1);
        nLines = size(rfLp, 2);
        zAxis = (1:nSamples)' * (c0_bf / fs) / 2;
        xAxis = (0:nLines-1) * 0.3e-3;
        xAxis = xAxis - mean(xAxis);

        save(fullfile(bfDir, sprintf('%s_f%d_LP.mat', sample, i)), 'rfLp', 'xAxis', 'zAxis', 'fs', 'm');
        save(fullfile(bfDir, sprintf('%s_f%d_HP.mat', sample, i)), 'rfHp', 'xAxis', 'zAxis', 'fs', 'm');

        fprintf('  Frame %d (m=%d) beamformeado y guardado en bf/\n', i, m);

    end

end

fprintf('\nListo. Carpetas "bf" generadas para cada phantom (CAPS).\n');
