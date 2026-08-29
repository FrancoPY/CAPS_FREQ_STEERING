sample = 'Phantom_inc_BA11_FrecSteering_v1';
probe = 'L14-5u';

% Frecuencias [Hz]
freqs = [5e6, 7e6, 9e6];

% Ángulos de steering [grados]
steeringAngles = [-5, 0, 5, 10, 15];

% m subaperturas
m = 2;

densityDir = fullfile(pwd, 'densityMaps', sample);

nFrames = 6;
dataCast = 'gpuArray-single';

% Grid
pmlXSize = 25;
pmlYSize = 25;

nx = 2000;
ny = 2000;

dx = 0.02e-3; % tamaño de un punto, [m]
dy = dx;

kgrid = kWaveGrid(nx, dx, ny, dy);

% Inclusion parameters
centerDepth = 22.5e-3;  % [m]
radius = 7.5e-3;        % [m]
inclusion = makeDisc(nx, ny, round(centerDepth/dx), ny/2, round(radius/dx)); % radius o centerDepth / dx, es en puntos grilla

% Medium
c0   = 1500.0;
rho0 = 1000;

a_bg  = 0.10;   y_bg  = 2.00;
a_inc = 0.10;   y_inc = 2.00;

f_ref_MHz = 7;  % frecuencia de referencia para el ajuste de atenuación
y_ref = 2.00;

a_bg_eff  = a_bg  * f_ref_MHz^(y_bg  - y_ref);
a_inc_eff = a_inc * f_ref_MHz^(y_inc - y_ref);

medium.BonA            = 6.00 * ones(nx, ny);
medium.sound_speed     = c0   * ones(nx, ny);
medium.sound_speed_ref = c0;

medium.alpha_power     = y_ref;
medium.alpha_coeff     = a_bg_eff * ones(nx, ny);
medium.alpha_coeff(inclusion > 0) = a_inc_eff;

medium.BonA        = medium.BonA + (11.00 - 6.00) * inclusion;
medium.sound_speed = medium.sound_speed + (1500.0 - c0) * inclusion;

medium.alpha_mode  = 'no_dispersion';

% Time grid
tEnd = 2.3 * (nx * dx) / c0;
kgrid.makeTime(c0, [], tEnd);

toneBurstCycles = 10;

% Probe 
numElem   = 128;
elemWidth = round(0.3e-3 / dx);
kerf      = round(0 / dx);
pitch     = elemWidth + kerf;
apWidth   = numElem * elemWidth + (numElem - 1) * kerf;

y0 = floor((ny - apWidth) / 2) + 1; %empieza en el 41

x_ap  = 5;
x_rcv = 5;

condition = numElem - round(numElem/m)*m==0; %para saber si se divide exacto, sino esta mal y manda error
if condition
    element_index = -(numElem-1)/2 : (numElem-1)/2;
    elementSpacing = pitch*dx;
    
    sourceStrength = 400e3; % una sola presión para todas las transmisiones
    
    nPtsPerElem = elemWidth;
    
    % Máscara para apertura completa
    src_mask_full = false(kgrid.Nx, kgrid.Ny);

    for e = 0:numElem-1
        ys = y0 + e * pitch;
        ye = ys + elemWidth - 1;
        src_mask_full(x_ap, ys:ye) = true;
    end
    
    % Máscara para las m subaperturas
    src_mask_sub = cell(1, m);
    for j=1:m
        active_Elems = j:m:numElem;

        maskM = false(kgrid.Nx, kgrid.Ny); % máscara de cada subapertura
        for e = active_Elems
            ys = y0 + (e-1) * pitch; %porque j empieza en 1, no en cero
            ye = ys + elemWidth - 1;
            maskM(x_ap, ys:ye) = true;
        end
        src_mask_sub{j} = maskM;
    end

    % Máscara para la recepción, siempre está activa, sin importar las
    % subaperturas
    rcv_mask = false(kgrid.Nx, kgrid.Ny);
    elem_label = zeros(kgrid.Nx, kgrid.Ny);
    for e = 0:numElem-1
        ys = y0 + e * pitch;
        ye = ys + elemWidth - 1;
        rcv_mask(x_rcv, ys:ye) = true;
        elem_label(x_rcv, ys:ye) = e + 1;
    end
    
    % Sensor
    sensor.mask = rcv_mask;
    sensor.directivity_angle = zeros(size(sensor.mask));
    sensor.directivity_size = 10 * kgrid.dx;
    sensor.record = {'p'};
    
    % k-Wave
    inputArgs = { ...
        'PMLInside',  false, ...
        'PMLSize',    [pmlXSize, pmlYSize], ...
        'DataCast',   dataCast, ...
        'DataRecast', true, ...
        'PlotSim',    false ...
    };
    
    labels = elem_label(sensor.mask);
    fs = 1 / kgrid.dt;
    
    % Barrido de frecuencia y ángulo
    for f0 = freqs
    
        freqStr = sprintf('%dMHz', round(f0/1e6));
        toneBurstFreq = f0;
    
        for steering_angle = steeringAngles
    
            fprintf('=== f0 = %d MHz, angulo = %d grados ===\n', round(f0/1e6), steering_angle);
    
            % Offset base y retardo por elemento (dependen de frecuencia y ángulo)
            max_offset = elementSpacing * max(abs(element_index)) * ...
                         abs(sin(steering_angle * pi/180)) / (c0 * kgrid.dt);
            offset_base = ceil(max_offset) + 20;
    
            tone_burst_offset = offset_base + elementSpacing * element_index * ...
                sin(steering_angle * pi/180) / (c0 * kgrid.dt);
    
            pulseNorm = toneBurst(1/kgrid.dt, toneBurstFreq, toneBurstCycles, ...
                'SignalOffset', tone_burst_offset);
    
            nPtsPerElem = elemWidth;
            pulseNormExpanded_full = zeros(numElem * nPtsPerElem, size(pulseNorm, 2));
            for e = 1:numElem
                rows = (e-1)*nPtsPerElem + 1 : e*nPtsPerElem;
                pulseNormExpanded_full(rows, :) = repmat(pulseNorm(e,:), nPtsPerElem, 1);
            end

            pulseNormExpanded_sub = cell(1, m);
            for j=1:m
                active_Elems = j:m:numElem;
                pulseM = zeros(numel(active_Elems) * nPtsPerElem, size(pulseNorm, 2));
                for idx = 1:numel(active_Elems)
                    e = active_Elems(idx);
                    rows = (idx-1)*nPtsPerElem + 1 : idx*nPtsPerElem;
                    pulseM(rows, :) = repmat(pulseNorm(e,:), nPtsPerElem, 1);
                end
                pulseNormExpanded_sub{j} = pulseM;
            end
    
            outputDir = fullfile(pwd, 'FRECUENCIA-STEERING', sample, probe, freqStr, 'rf', ...
                sprintf('Angle_%d', steering_angle));
            if ~exist(outputDir, 'dir'); mkdir(outputDir); end
    
            % Frames
            for frame = 1:nFrames
    
                densityFile = sprintf('frame%d', frame);
                load(fullfile(densityDir, densityFile));
    
                medium.density = single(density);
    
                % m transmisiones, no la completa
                for j=1:m
                    fprintf('[f0=%dMHz][angulo %d][frame %d] subapertura%d/%d\n', ...
                        round(f0/1e6), steering_angle, frame, j, m);
    
                    vdrive = single((sourceStrength / (c0 * rho0)) * pulseNormExpanded_sub{j});
    
                    source.ux = vdrive;
                    source.u_mask = src_mask_sub{j};
    
                    sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputArgs{:});
    
                    rf_allpts = gather(sensorData.p);
    
                    rfElem = zeros(numElem, size(rf_allpts, 2), 'like', rf_allpts);
                    for e = 1:numElem
                        rfElem(e, :) = mean(rf_allpts(labels == e, :), 1);
                    end
    
                    rfPrebf = rfElem.';
                    outName = fullfile(outputDir, ...
                        sprintf('%s_f%d_sub%d_of_m%d_%dkPa.mat', sample, frame, j, m, round(sourceStrength/1e3)));
    
                    save(outName, 'rfPrebf', 'fs', 'c0', 'offset_base', 'f0', 'steering_angle','m','j');
    
                    source = struct();
                    clear sensorData rf_allpts rfElem rfPrebf
    
                end
    
                % Transmisión de apertura completa
                fprintf('  [f0=%dMHz][angulo %d][frame %d] apertura completa\n', ...
                    round(f0/1e6), steering_angle, frame);

                 vdrive = single((sourceStrength / (c0 * rho0)) * pulseNormExpanded_full);

                source.ux = vdrive;
                source.u_mask = src_mask_full;
            
                sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputArgs{:});
            
                rf_allpts = gather(sensorData.p);
            
                rfElem = zeros(numElem, size(rf_allpts, 2), 'like', rf_allpts);
                for e = 1:numElem
                    rfElem(e, :) = mean(rf_allpts(labels == e, :), 1);
                end
            
                rfPrebf = rfElem.';
                outName = fullfile(outputDir, ...
                    sprintf('%s_f%d_full_m%d_%dkPa.mat', sample, frame, m, round(sourceStrength / 1e3)));
            
                save(outName, 'rfPrebf', 'fs', 'c0', 'offset_base', 'f0', 'steering_angle', 'm');
            
                source = struct();
                clear sensorData rf_allpts rfElem rfPrebf
            
                clear density
            end
    
        end
    
    end
else
    fprintf('m debe dividir de manera exacta el número de elementos');
end
