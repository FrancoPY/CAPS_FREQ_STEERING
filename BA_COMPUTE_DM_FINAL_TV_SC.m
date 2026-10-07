refSample = 'Simulation_Bgnd_CAPS_FreqSteering_v4';
samSample = 'Simulation_Inc_CAPS_FreqSteering_v4';
freq_vect = 5e6;

probe = 'L14-5u';
basedir = fullfile(pwd,'CAPS-FRECUENCIA-STEERING');
angles = {'Angle_-15'};

nFrames = 1;
    
c0 = 1500;
B_r = 4;     % beta de referencia (calibrado con B/A_ref = 6)

zMaxProc = 55e-3;
mu = 0.1;
invParam.zCrop = [10e-3 51e-3];  invParam.zInv = [10e-3 51e-3];   % se invierte hasta 51 mm...
invParam.zOut  = [10e-3 50e-3];   
invParam.gridSize = 0.3e-3;      invParam.plotFlag = false;
BAL_sum = [];  BAL_count = [];

P_ref_L_acumu = 0;
P_ref_H_acumu = 0;
P_sam_L_acumu = 0;
P_sam_H_acumu = 0;

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
        thetaDeg = angleStrToDeg(angleStr);
        fs = ref_L0.fs;
        N = min(round(2*zMaxProc*fs/c0), size(ref_L0.rfLp,1));
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
        
        B_A(B_A < -3 | B_A > 25) = NaN;      % descarta valores fuera de rango físico razonable
        B_A = fillmissing(B_A, 'linear', 1); % rellena esos NaN interpolando en la dirección axial
        B_A = fillmissing(B_A, 'linear', 2); % rellena lo que quede, interpolando lateralmente

        % Suavizado axial y lateral
        lambda = c0/f_fund;
        muestras_ventana = round((2*10*lambda*fs)/c0);
        B_suave_axial = movmean(B_A, muestras_ventana, 1);
        B_final = movmean(B_suave_axial, 4, 2);
        
        BAeff_theta.image = B_final;  BAeff_theta.lateral = x;  BAeff_theta.axial = z;
        BAL_theta = invertLocalBA_Steered(BAeff_theta, mu, invParam, thetaDeg);
        
        if isempty(BAL_sum)
            BAL_sum = zeros(size(BAL_theta.image));  BAL_count = BAL_sum;
            x_crop = BAL_theta.lateral;  z_crop = BAL_theta.axial;
        end
        ok = ~isnan(BAL_theta.image);  tmp = BAL_theta.image;  tmp(~ok) = 0;
        BAL_sum = BAL_sum + tmp;  BAL_count = BAL_count + ok;
    
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

BA_final = BAL_sum ./ BAL_count;
fprintf('B/A compounded: min=%.2f max=%.2f, #NaN=%d\n', ...
    min(BA_final(:),[],'omitnan'), max(BA_final(:),[],'omitnan'), sum(isnan(BA_final(:))));
muStr = sprintf('%g', mu);

if numel(angles) == 1
    angStr = angles{1};   
else
    angStr = sprintf('SC%dang', numel(angles));
end
    
% Recorte de las imágenes
idxStart = find(z >= 10e-3, 1, 'first');
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
figDir = fullfile(basedir, 'figuras_v4_TV_m2');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% Visualización B/A
fig1 = figure('Visible', 'off');   % <-- 'off' porque no hay pantalla en el cluster
imagesc(x_crop*1000, z_crop*1000, BA_final);
axis image; colormap("pink"); colorbar;
title(sprintf('Mapa B/A local (TV, \\mu=%.2g), Freq %s %s', mu,freqStr,angStr));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
clim([5 12]);

% Guardar como PNG
outNamePNG = fullfile(figDir, sprintf('BA_TV_%s_%s_mu%s_V4_Pink1.png',freqStr,angStr,muStr));
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
outNamePNG2 = fullfile(figDir, sprintf('B-mode_TV_%s_%s_mu%s_V4_Pink1.png',freqStr,angStr,muStr));
saveas(fig2, outNamePNG2);
close(fig2);
fprintf('Figura B-mode guardada en: %s\n', outNamePNG2);


function thetaDeg = angleStrToDeg(s)
    thetaDeg = str2double(strrep(s,'Angle_',''));
end