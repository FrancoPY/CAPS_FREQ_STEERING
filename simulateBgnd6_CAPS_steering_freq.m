sample = 'Phantom_bgnd_BA6_CAPS_FreqSteering_v2';
probe = 'L14-5u';

% Frecuencias [Hz]
freqs = [4e6, 5e6, 6e6];

% Ángulos de steering [grados]
steeringAngles = [-5, 0, 5, 10, 15];

m = 2;

densityDir = fullfile(pwd, 'densityMaps', sample);

nFrames = 1;
dataCast = 'gpuArray-single';

% Grid
pmlXSize = 25;
pmlYSize = 25;

nx = 2000;
ny = 2000;

dx = 0.02e-3;
dy = dx;

kgrid = kWaveGrid(nx, dx, ny, dy);

% Medium
c0   = 1500.0;
rho0 = 1000;

medium.BonA            = 6.00;
medium.sound_speed     = c0;
medium.sound_speed_ref = c0;
medium.alpha_coeff     = 0.1;
medium.alpha_power     = 2.00;
medium.alpha_mode      = 'no_dispersion';

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

y0 = floor((ny - apWidth) / 2) + 1;

condition = numElem - round(numElem/m)*m==0;

if condition
    x_ap  = 5;
    x_rcv = 5;
    
    element_index = -(numElem-1)/2 : (numElem-1)/2;
    elementSpacing = 0.3e-3;
    
    sourceStrength = 400e3;
    
    % Máscaras de apertura completa
    src_mask_full = false(kgrid.Nx, kgrid.Ny);
    for e = 0:numElem-1
        ys = y0 + e * pitch;
        ye = ys + elemWidth - 1;
        src_mask_full(x_ap, ys:ye) = true;
    end
    
    %Máscara de m subaperturas
    src_mask_sub = cell(1,m);
    for j=1:m
        active_Elems = j:m:numElem;
        src_mask_M = false(kgrid.Nx, kgrid.Ny);
        for e=active_Elems
            ys = y0 + (e-1)* elemWidth;
            ye= ys + elemWidth-1;
            src_mask_M(x_ap, ys:ye) = true;
        end
        src_mask_sub{j} = src_mask_M;
    end
    
    % Máscara de recepción, debe estar completa simepre
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
    
            fprintf('f0 - %d MHz, angulo - %d grados\n', round(f0/1e6), steering_angle);
    
            % Offset base y retardo por elemento (dependen de frecuencia y ángulo)
            max_offset = elementSpacing * max(abs(element_index)) * ...
                         abs(sin(steering_angle * pi/180)) / (c0 * kgrid.dt);
            offset_base = ceil(max_offset) + 20;
    
            tone_burst_offset = offset_base + elementSpacing * element_index * ... %que es un vector de 128 valores
                sin(steering_angle * pi/180) / (c0 * kgrid.dt);
    
            pulseNorm = toneBurst(1/kgrid.dt, toneBurstFreq, toneBurstCycles, ...
                'SignalOffset', tone_burst_offset); %(128, 15333), necesitamos expandirlo
    
            nPtsPerElem = elemWidth;
            pulseNormExpanded_full = zeros(numElem * nPtsPerElem, size(pulseNorm, 2)); % (1920,15333)
            for e = 1:numElem
                rows = (e-1)*nPtsPerElem + 1 : e*nPtsPerElem;
                pulseNormExpanded_full(rows, :) = repmat(pulseNorm(e,:), nPtsPerElem, 1);
            end
            
            % Pulso para las subaperturas
            pulseNormExpanded_sub = cell(1,m); % (1, 2)
            for j = 1:m
                active_Elems = j:m:numElem;
                pulseNormExpanded_M = zeros(numel(active_Elems) * nPtsPerElem, size(pulseNorm, 2)); %(960,15333)
                for idx = 1:numel(active_Elems) 
                    e = active_Elems(idx);
                    rows = (idx-1)*nPtsPerElem + 1 : idx*nPtsPerElem;
                    pulseNormExpanded_M(rows, :) = repmat(pulseNorm(e,:), nPtsPerElem, 1);
                end
                pulseNormExpanded_sub{j} = pulseNormExpanded_M;
            end
    
            outputDir = fullfile(pwd, 'CAPS-FRECUENCIA-STEERING', sample, probe, freqStr, 'rf', ...
                sprintf('Angle_%d', steering_angle));
            if ~exist(outputDir, 'dir'); mkdir(outputDir); end
    
            % Frames
            for frame = 1:nFrames
    
                densityFile = sprintf('frame%d', frame);
                load(fullfile(densityDir, densityFile));
    
                medium.density = single(density);
                for j=1:m
                    fprintf('  [f0=%dMHz][angulo %d][frame %d] subapertura%d/%d\n', ...
                        round(f0/1e6), steering_angle, frame, j, m);
    
                    vdrive = single((sourceStrength / (c0 * rho0)) * pulseNormExpanded_sub{j});
    
                    source.ux = vdrive;
                    source.u_mask = src_mask_sub{j};
        
                    sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputArgs{:});
        
                    rf_allpts = gather(sensorData.p);
    
                    % Esto promedia los 1920 puntos a 128
                    rfElem = zeros(numElem, size(rf_allpts, 2), 'like', rf_allpts);
                    for e = 1:numElem
                        rfElem(e, :) = mean(rf_allpts(labels == e, :), 1);
                    end
    
                    rfPrebf = rfElem.';
                    outName = fullfile(outputDir, ...
                        sprintf('%s_f%d_sub%d_of_m%d_%dkPa.mat', sample, frame, j, m, round(sourceStrength / 1e3)));
        
                    save(outName, 'rfPrebf', 'fs', 'c0', 'offset_base', 'f0', 'steering_angle', 'm', 'j');
        
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
    
                % Esto promedia los 1920 puntos a 128
                rfElem = zeros(numElem, size(rf_allpts, 2), 'like', rf_allpts);
                for e = 1:numElem
                    rfElem(e, :) = mean(rf_allpts(labels == e, :), 1);
                end
    
                rfPrebf = rfElem.';
                outName = fullfile(outputDir, ...
                    sprintf('%s_f%d_full_m%d_%dkPa.mat', sample, frame, m, round(sourceStrength / 1e3)));
    
                save(outName, 'rfPrebf', 'fs', 'c0', 'offset_base', 'f0', 'steering_angle','m');
    
                source = struct();
                clear sensorData rf_allpts rfElem rfPrebf
    
                clear density
    
            end
    
        end
    
    end
else
    fprintf('m debe dividir de manera exacta el número de elementos');
end