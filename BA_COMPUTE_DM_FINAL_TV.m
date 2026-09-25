refSample = 'Simulation_Bgnd_CAPS_FreqSteering_v3';
samSample = 'Simulation_Inc_CAPS_FreqSteering_v3';
freq_vect = [4e6, 5e6, 6e6];


probe = 'L14-5u';
basedir = fullfile(pwd,'CAPS-FRECUENCIA-STEERING');
angles = {'Angle_-5', 'Angle_0', 'Angle_5', 'Angle_10', 'Angle_15'};

nFrames = 1;
    
c0 = 1500;
B_r = 4;     % beta de referencia (calibrado con B/A_ref = 6)

BA_mean = 0;
P_ref_L_acumu = 0;
P_ref_H_acumu = 0;
P_sam_L_acumu = 0;
P_sam_H_acumu = 0;
N = 18325;

totalCombinaciones = length(freq_vect) * length(angles);   

for freq = freq_vect
    f_fund = freq;
    f_fund_r = freq;
    
    alpha = 0.10*((f_fund/1e6)^2)*100/8.686;
    alpha_r = 0.10*((f_fund_r/1e6)^2)*100/8.686;
    freqStr = sprintf('%dMHz',freq/1e6);

    for a=1:length(angles)
        angleStr = angles{a};
        
        refDir = fullfile(basedir,refSample, probe, freqStr, 'bf', angleStr);
        samDir = fullfile(basedir,samSample, probe, freqStr, 'bf', angleStr);
    
        % Cargar la referencia
        ref_L0 = load(fullfile(refDir, sprintf('%s_f1_LP.mat',refSample)));
        fs = ref_L0.fs;
        zAxis = ref_L0.zAxis;
        xAxis = ref_L0.xAxis;
        
        bw = 0.8e6;
        f_low = f_fund - bw/2;
        f_high = f_fund + bw/2;
        Wn = [f_low f_high] / (fs/2);
        order = 200;
        hFilter = fir1(order, Wn);
        
        z = (1:N)' * (c0 / fs) / 2;
        x = xAxis;
        term = ((1-exp(-2*alpha_r*z))./(1-exp(-2*alpha*z)))*(alpha/alpha_r);
        
        P_ref_L_sum = 0;
        P_ref_H_sum = 0;
        P_sam_L_sum = 0;
        P_sam_H_sum = 0;
        
        for i = 1:nFrames
            % Cargar
            rL = load(fullfile(refDir, sprintf('%s_f%d_LP.mat', refSample, i)));
            rH = load(fullfile(refDir, sprintf('%s_f%d_HP.mat', refSample, i)));
            sL = load(fullfile(samDir, sprintf('%s_f%d_LP.mat', samSample, i)));
            sH = load(fullfile(samDir, sprintf('%s_f%d_HP.mat', samSample, i)));
        
            ref_L    = rL.rfLp;   % suma coherente de subaperturas (ν·P_L)
            ref_H    = rH.rfHp;   % apertura completa
            sample_L = sL.rfLp;
            sample_H = sH.rfHp;
        
            % Filtrado
            ref_L_f = filtfilt(hFilter, 1, ref_L(1:N,:));
            ref_H_f = filtfilt(hFilter, 1, ref_H(1:N,:));
            sample_L_f = filtfilt(hFilter, 1, sample_L(1:N,:));
            sample_H_f = filtfilt(hFilter, 1, sample_H(1:N,:));
        
            % Envolventes
            P_ref_L_sum = P_ref_L_sum + abs(hilbert(ref_L_f));
            P_ref_H_sum = P_ref_H_sum + abs(hilbert(ref_H_f));
            P_sam_L_sum = P_sam_L_sum + abs(hilbert(sample_L_f));
            P_sam_H_sum = P_sam_H_sum + abs(hilbert(sample_H_f));
        
            if i == 1
                P_ref_L_frame1 = abs(hilbert(ref_L_f));
                P_ref_H_frame1 = abs(hilbert(ref_H_f));
                P_sam_L_frame1 = abs(hilbert(sample_L_f));
                P_sam_H_frame1 = abs(hilbert(sample_H_f));
            end
        
        end
        
        P_ref_L_mean = P_ref_L_sum / nFrames;
        P_ref_H_mean = P_ref_H_sum / nFrames;
        P_sam_L_mean = P_sam_L_sum / nFrames;
        P_sam_H_mean = P_sam_H_sum / nFrames;
        
        % Calcular B/A (CAPS: sin factor v, ya está implícito en la suma coherente)
        numerador = (P_sam_L_mean - P_sam_H_mean).*(P_ref_L_mean);
        denominador = (P_ref_L_mean - P_ref_H_mean).*(P_sam_L_mean);
        ratio = numerador ./ denominador;
        B = B_r * sqrt(abs(ratio)) .* term;
        B_A = 2*(B - 1);
        
        BA_mean = BA_mean + B_A;
    
        P_ref_L_acumu = P_ref_L_acumu + P_ref_L_frame1;
        P_ref_H_acumu = P_ref_H_acumu + P_ref_H_frame1;
        P_sam_L_acumu = P_sam_L_acumu + P_sam_L_frame1;
        P_sam_H_acumu = P_sam_H_acumu + P_sam_H_frame1;
    end
end

P_ref_L_meanAngles = P_ref_L_acumu / totalCombinaciones;
P_ref_H_meanAngles = P_ref_H_acumu / totalCombinaciones;
P_sam_L_meanAngles = P_sam_L_acumu / totalCombinaciones;
P_sam_H_meanAngles = P_sam_H_acumu / totalCombinaciones;

BA_mean = BA_mean / totalCombinaciones;

% ---- Inversión local con TV (reemplaza el recorte manual del B/A) ----
BAeff.image   = BA_mean;   % B/AC completo
BAeff.lateral = x;         % [m]
BAeff.axial   = z;         % [m], eje completo (1:N)

invParam.zCrop    = [3.5e-3, 50e-3];   % <-- aquí defines tu recorte 3.5-50mm
invParam.zInv     = [3.5e-3, 50e-3];
invParam.gridSize = 0.3e-3;
invParam.plotFlag = false;

mu = 1;
BAloc = invertBAeffTV(BAeff, mu, invParam);

BA_final = BAloc.image;   % mapa LOCAL con TV, ya recortado 3.5-50mm
z_crop   = BAloc.axial;   % eje z que corresponde a BA_final (nueva grilla)
x_crop   = BAloc.lateral; % eje x que corresponde a BA_final
    
% Recorte de las imágenes
idxStart = find(z >= 3.5e-3, 1, 'first');
idxEnd   = find(z <= 50e-3, 1, 'last'); 
z_crop_bmode = z(idxStart:idxEnd);

P_ref_L_meanAngles = P_ref_L_meanAngles(idxStart:idxEnd, :);
P_ref_H_meanAngles = P_ref_H_meanAngles(idxStart:idxEnd, :);
P_sam_L_meanAngles = P_sam_L_meanAngles(idxStart:idxEnd, :);
P_sam_H_meanAngles = P_sam_H_meanAngles(idxStart:idxEnd, :);

refMax = max(P_ref_L_meanAngles(:));
P_ref_L_Bmode = 20*log10(P_ref_L_meanAngles / refMax);
P_ref_H_Bmode = 20*log10(P_ref_H_meanAngles / refMax);
P_sam_L_Bmode = 20*log10(P_sam_L_meanAngles / refMax);
P_sam_H_Bmode = 20*log10(P_sam_H_meanAngles / refMax);

% Crear carpeta para guardar las figuras
figDir = fullfile(basedir, 'figuras_v3_TV');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% Visualización B/A
fig1 = figure('Visible', 'off');   % <-- 'off' porque no hay pantalla en el cluster
imagesc(x_crop*1000, z_crop*1000, BA_final);
axis image; colormap(turbo); colorbar;
title(sprintf('Mapa B/A local (TV, \\mu=%.2g), todo promediado', mu));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
clim([5 12]);

% Guardar como PNG
outNamePNG = fullfile(figDir, sprintf('BA_TV_all mean_mu1.png'));
saveas(fig1, outNamePNG);
close(fig1);
fprintf('Figura B/A guardada en: %s\n', outNamePNG);


% Visualización B-mode
dynamicRange = 60;                                      

fig2 = figure('Visible', 'off');
subplot(2,2,1); imagesc(x*1000, z_crop_bmode*1000, P_ref_L_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title(sprintf('B-mode Referencia Low'));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,2); imagesc(x*1000, z_crop_bmode*1000, P_ref_H_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title(sprintf('B-mode Referencia High'));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,3); imagesc(x*1000, z_crop_bmode*1000, P_sam_L_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title(sprintf('B-mode Muestra Low'));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,4); imagesc(x*1000, z_crop_bmode*1000, P_sam_H_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title(sprintf('B-mode Muestra High'));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

% Guardar como PNG
outNamePNG2 = fullfile(figDir, sprintf('Bmode_all_angles_all_freq.png_mu1'));
saveas(fig2, outNamePNG2);
close(fig2);
fprintf('Figura B-mode guardada en: %s\n', outNamePNG2);