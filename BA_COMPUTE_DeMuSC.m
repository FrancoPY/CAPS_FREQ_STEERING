refSample = 'Simulation_Bgnd_CAPS_FreqSteering_v5_m2';
samSample = 'Simulation_Inc_CAPS_FreqSteering_v5_m2';
freq_vect = [5e6];
 
probe = 'L14-5u';
basedir = fullfile(pwd,'CAPS-FRECUENCIA-STEERING');
angles = {'Angle_0', 'Angle_-15', 'Angle_-10', 'Angle_-5', 'Angle_5', 'Angle_10', 'Angle_15'};
 
nFrames = 1;
 
c0 = 1500;
B_r = 4;     % beta de referencia (calibrado con B/A_ref = 6)
zMaxProc = 55e-3;
 
% Parametros del DeMuSC (funciones prepDeMuSc, calcScRayOp, invertDeMuSc, solveDeMuScTv)
mu = 10^-0.5;
param.gridSize = 0.6e-3;
param.zMin     = 10e-3;     % desde aqui se usan observaciones
param.zMax     = 50e-3;     % hasta aqui se usan observaciones
param.mu       = mu;
param.tau      = 0.1;
param.maxIter  = 7;
 
muStr = sprintf('%.2g', mu);
nAng = numel(angles);
nFreq = numel(freq_vect);
angNum = zeros(1, nAng);
for a = 1:nAng
    angNum(a) = angleStrToDeg(angles{a});
end
 
% Crear carpeta para guardar las figuras
figDir = fullfile(basedir, 'figuras_v5_DeMuSC_m2');
if ~exist(figDir, 'dir'); mkdir(figDir); end
 
BAfreq = [];           % aqui se guarda el B/A_L de cada frecuencia
P_ref_L_acumu = 0;
P_ref_H_acumu = 0;
P_sam_L_acumu = 0;
P_sam_H_acumu = 0;
iF = 0;
 
for freq = freq_vect
    iF = iF + 1;
    f_fund = freq;
    f_fund_r = freq;
    alpha = 0.10*((f_fund/1e6)^2)*100/8.686;
    alpha_r = 0.10*((f_fund_r/1e6)^2)*100/8.686;
    freqStr = sprintf('%dMHz',freq/1e6);
 
    % Datos que no dependen del angulo (se leen del primer angulo)
    refDir0 = fullfile(basedir, refSample, probe, freqStr, 'bf', angles{1});
    ref0 = load(fullfile(refDir0, sprintf('%s_f1_LP.mat', refSample)));
    fs = ref0.fs;
    x = ref0.xAxis;
    N = min(round(2*zMaxProc*fs/c0), size(ref0.rfLp,1));
    clear ref0
 
    z = (1:N)' * (c0/fs) / 2;
    term = ((1-exp(-2*alpha_r*z))./(1-exp(-2*alpha*z)))*(alpha/alpha_r);
 
    bw = 0.8e6;
    Wn = [f_fund-bw/2, f_fund+bw/2] / (fs/2);
    hFilter = fir1(200, Wn);
 
    lambda = c0/f_fund;
    muestras_ventana = round((2*10*lambda*fs)/c0);
 
    baStack = zeros(N, numel(x), nAng);   % B/A_C de cada angulo para esta frecuencia
 
    for a=1:nAng
        angleStr = angles{a};
 
        refDir = fullfile(basedir,refSample, probe, freqStr, 'bf', angleStr);
        samDir = fullfile(basedir,samSample, probe, freqStr, 'bf', angleStr);
 
        P_ref_L_sum = 0;
        P_ref_H_sum = 0;
        P_sam_L_sum = 0;
        P_sam_H_sum = 0;
 
        for i = 1:nFrames
            % Cargar, filtrar y sacar envolvente
            refL = leerEnvolvente(fullfile(refDir, sprintf('%s_f%d_LP.mat', refSample, i)), 'rfLp', N, hFilter);
            refH = leerEnvolvente(fullfile(refDir, sprintf('%s_f%d_HP.mat', refSample, i)), 'rfHp', N, hFilter);
            samL = leerEnvolvente(fullfile(samDir, sprintf('%s_f%d_LP.mat', samSample, i)), 'rfLp', N, hFilter);
            samH = leerEnvolvente(fullfile(samDir, sprintf('%s_f%d_HP.mat', samSample, i)), 'rfHp', N, hFilter);
 
            P_ref_L_sum = P_ref_L_sum + refL;
            P_ref_H_sum = P_ref_H_sum + refH;
            P_sam_L_sum = P_sam_L_sum + samL;
            P_sam_H_sum = P_sam_H_sum + samH;
 
            % Para el B-mode solo se usa el frame 1
            if i == 1
                P_ref_L_acumu = P_ref_L_acumu + refL;
                P_ref_H_acumu = P_ref_H_acumu + refH;
                P_sam_L_acumu = P_sam_L_acumu + samL;
                P_sam_H_acumu = P_sam_H_acumu + samH;
            end
        end
 
        P_ref_L_mean = P_ref_L_sum / nFrames;
        P_ref_H_mean = P_ref_H_sum / nFrames;
        P_sam_L_mean = P_sam_L_sum / nFrames;
        P_sam_H_mean = P_sam_H_sum / nFrames;
 
        % B/A acumulado (B/A_C)
        numerador = (P_sam_L_mean - P_sam_H_mean).*(P_ref_L_mean);
        denominador = (P_ref_L_mean - P_ref_H_mean).*(P_sam_L_mean);
        ratio = numerador ./ denominador;
        B = B_r * sqrt(abs(ratio)) .* term;
        B_A = 2*(B - 1);
 
        B_A(B_A < -3 | B_A > 25) = NaN;
        B_A = fillmissing(B_A, 'linear', 1);
        B_A = fillmissing(B_A, 'linear', 2);
 
        % Suavizado axial y lateral
        B_suave_axial = movmean(B_A, muestras_ventana, 1);
        B_final = movmean(B_suave_axial, 4, 2);
 
        baStack(:,:,a) = B_final;
    end
 
    % DeMuSC para esta frecuencia: un solo sistema con todos los angulos
    prob = prepDeMuSc(baStack, x, z, angNum, param);
    res = invertDeMuSc(prob, param);
    fprintf('%s | DeMuSC convergio: %d (iter %d), init flag %d\n', freqStr, res.diag.converged, res.diag.nIter, res.diag.init.flag);
 
    % Quedarse con 10 a 50 mm
    filas = res.z >= 10e-3 & res.z <= 50e-3;
    z_crop = res.z(filas);
    x_crop = res.x;
    BAfreq(:,:,iF) = res.image(filas, :);
end
 
% Promedio de los angulos y de las frecuencias
totalCombinaciones = nFreq * nAng;
P_ref_L_meanAngles = P_ref_L_acumu / totalCombinaciones;
P_ref_H_meanAngles = P_ref_H_acumu / totalCombinaciones;
P_sam_L_meanAngles = P_sam_L_acumu / totalCombinaciones;
P_sam_H_meanAngles = P_sam_H_acumu / totalCombinaciones;
 
BA_final = mean(BAfreq, 3, 'omitnan');
 
fprintf('B/A compounded: min=%.2f max=%.2f, #NaN=%d\n', ...
    min(BA_final(:),[],'omitnan'), max(BA_final(:),[],'omitnan'), sum(isnan(BA_final(:))));
 
% ROIs (en mm): cuadrado dentro de la inclusión y el mismo cuadrado partido en fondo izq y der
zc = 22.5;           % centro de la inclusión en profundidad
xc = 0;              % centro lateral de la inclusión
R = 7.5;             % radio de la inclusión
L = 8;               % lado del cuadrado (debe caber dentro de la inclusión)
gap = 3;             % separación entre el borde de la inclusión y el fondo
 
xi = [xc-L/2, xc+L/2];                    % inclusión
zi = [zc-L/2, zc+L/2];
xbL = [xc-R-gap-L/2, xc-R-gap];           % fondo izquierdo (ancho L/2)
xbR = [xc+R+gap, xc+R+gap+L/2];           % fondo derecho (ancho L/2)
zb = zi;
 
X = x_crop(:).'*1000;   Z = z_crop(:)*1000;
[XX, ZZ] = meshgrid(X, Z);
maskInc = XX>=xi(1)  & XX<=xi(2)  & ZZ>=zi(1) & ZZ<=zi(2);
maskBg  = ((XX>=xbL(1) & XX<=xbL(2)) | (XX>=xbR(1) & XX<=xbR(2))) & ZZ>=zb(1) & ZZ<=zb(2);
 
vInc = BA_final(maskInc);  vInc = vInc(~isnan(vInc));
vBg  = BA_final(maskBg);   vBg  = vBg(~isnan(vBg));
 
muInc = mean(vInc);  sdInc = std(vInc);
muBg  = mean(vBg);   sdBg  = std(vBg);
CNR   = abs(muInc - muBg) / sqrt(sdInc^2 + sdBg^2);
 
fprintf('Inclusión: %.2f +- %.2f (n=%d)\n', muInc, sdInc, numel(vInc));
fprintf('Fondo:     %.2f +- %.2f (n=%d)\n', muBg, sdBg, numel(vBg));
fprintf('Diferencia: %.2f | CNR: %.2f\n', muInc-muBg, CNR);
 
% Nombres para los archivos
if nAng == 1
    angStr = angles{1};
else
    angStr = sprintf('SC%dang', nAng);
end
if nFreq == 1
    freqName = sprintf('%dMHz', freq_vect(1)/1e6);
else
    freqName = sprintf('%dfreq', nFreq);
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
 
% Visualización B/A
fig1 = figure('Visible', 'off');
imagesc(x_crop*1000, z_crop*1000, BA_final, 'AlphaData', ~isnan(BA_final));
axis image; colormap("pink"); colorbar;
title(sprintf('Mapa B/A DeMuSC (TV, \\mu=%.2g), Freq %s %s', mu, freqName, angStr));
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
clim([5 12]);
 
hold on;
rectangle('Position',[xi(1) zi(1) L L], 'EdgeColor','c','LineWidth',1.5);
rectangle('Position',[xbL(1) zb(1) L/2 L], 'EdgeColor','w','LineWidth',1.5);
rectangle('Position',[xbR(1) zb(1) L/2 L], 'EdgeColor','w','LineWidth',1.5);
hold off;
 
% Guardar como PNG
outNamePNG = fullfile(figDir, sprintf('BA_DeMuSC_cuadrado_%s_%s_mu%s_V6_Pink.png',freqName,angStr,muStr));
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
outNamePNG2 = fullfile(figDir, sprintf('B-mode_DeMuSC_cuadrado_%s_%s_mu%s_V6_Pink.png',freqName,angStr,muStr));
saveas(fig2, outNamePNG2);
close(fig2);
fprintf('Figura B-mode guardada en: %s\n', outNamePNG2);
 
 
function thetaDeg = angleStrToDeg(s)
    thetaDeg = str2double(strrep(s,'Angle_',''));
end
 
function env = leerEnvolvente(archivo, campo, N, hFilter)
    % Carga el RF, lo filtra y devuelve la envolvente (primeras N muestras)
    datos = load(archivo, campo);
    rf = datos.(campo);
    rf = filtfilt(hFilter, 1, rf(1:N,:));
    env = abs(hilbert(rf));
end
