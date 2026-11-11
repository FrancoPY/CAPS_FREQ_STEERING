refSample = 'Simulation_Bgnd_CAPS_FreqSteering_v6_m2';
samSample = 'Simulation_Inc_CAPS_FreqSteering_v6_m2';
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
figDir = fullfile(basedir, 'figuras_v6.2_DeMuSC_m2');
if ~exist(figDir, 'dir'); mkdir(figDir); end
 
BAfreq = [];           % aqui se guarda el B/A_L de cada frecuencia
BAfreqNo = [];         % Para los mapas sin TVR
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
            if i == 1 && angNum(a) == 0
                P_ref_L_acumu = refL;
                P_ref_H_acumu = refH;
                P_sam_L_acumu = samL;
                P_sam_H_acumu = samH;
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

    % Mapa de entrada (B/A_C antes de invertir)
    filasIn = z >= 10e-3 & z <= 50e-3;
    baMean  = mean(baStack, 3);          % promedio de ángulos, solo para mirar
    a0 = find(angNum == 0);
    
    figure;
    subplot(1,2,1); imagesc(x*1000, z(filasIn)*1000, baStack(filasIn,:,a0));
    axis image; colormap("pink"); colorbar; clim([5 12]); title('B/A_C ángulo 0°');
    subplot(1,2,2); imagesc(x*1000, z(filasIn)*1000, baMean(filasIn,:));
    axis image; colormap("pink"); colorbar; clim([5 12]); title('B/A_C promedio ángulos');
    
    % Perfil en profundidad: centro de la inclusión vs. fondo
    figure; hold on;
    plot(z*1000, mean(baMean(:, abs(x) < 2e-3), 2), 'DisplayName', 'centro (inclusión)');
    plot(z*1000, mean(baMean(:, abs(x) > 10e-3), 2), 'DisplayName', 'fondo');
    xlabel('Profundidad (mm)'); ylabel('B/A_C'); legend; grid on;
 
    % DeMuSC para esta frecuencia: un solo sistema con todos los angulos
    prob = prepDeMuSc(baStack, x, z, angNum, param);
    res = invertDeMuSc(prob, param);
    fprintf('%s | DeMuSC convergio: %d (iter %d), init flag %d\n', freqStr, res.diag.converged, res.diag.nIter, res.diag.init.flag);
 
    % Quedarse con 10 a 50 mm
    filas = res.z >= 10e-3 & res.z <= 50e-3;
    z_crop = res.z(filas);
    x_crop = res.x;
    BAfreq(:,:,iF) = res.image(filas, :);
    BAfreqNo(:,:,iF) = res.noReg(filas, :);   % sin TVR
end
 
BA_final = mean(BAfreq, 3, 'omitnan');
BA_finalNo = mean(BAfreqNo, 3, 'omitnan'); %SIN TVR
 
fprintf('B/A compounded: min=%.2f max=%.2f, #NaN=%d\n', ...
    min(BA_final(:),[],'omitnan'), max(BA_final(:),[],'omitnan'), sum(isnan(BA_final(:))));
 
% ROIs (en mm)
zc = 22.5;  xc = 0;  R = 7.5;        % inclusión (círculo)
xR = [-18 18];  zR = [12 47];        % cuadrado rojo (región total)
margen = 0;                          % mm de separación alrededor de la inclusión

X = x_crop(:).'*1000;   Z = z_crop(:)*1000;
[XX, ZZ] = meshgrid(X, Z);
dist = sqrt((XX-xc).^2 + (ZZ-zc).^2);
maskInc = dist <= R;
maskBg  = XX>=xR(1) & XX<=xR(2) & ZZ>=zR(1) & ZZ<=zR(2) & dist > R + margen;

mapas    = {BA_final, BA_finalNo};
etiqueta = {'TV', 'sinTVR'};
 
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

for k = 1:2
    BA = mapas{k};

    vInc = BA(maskInc);  vInc = vInc(~isnan(vInc));
    vBg  = BA(maskBg);   vBg  = vBg(~isnan(vBg));
    
    muInc = mean(vInc);  sdInc = std(vInc);
    muBg  = mean(vBg);   sdBg  = std(vBg);
    CNR   = abs(muInc - muBg) / sqrt(sdInc^2 + sdBg^2);
    
    fprintf('Inclusión: %.2f +- %.2f (n=%d)\n', muInc, sdInc, numel(vInc));
    fprintf('Fondo:     %.2f +- %.2f (n=%d)\n', muBg, sdBg, numel(vBg));
    fprintf('Diferencia: %.2f | CNR: %.2f\n', muInc-muBg, CNR);

    % Visualización B/A
    fig1 = figure('Visible','off');
    imagesc(x_crop*1000, z_crop*1000, BA, 'AlphaData', ~isnan(BA));
    axis image; colormap("pink"); colorbar; clim([5 12]);
    title(sprintf('DeMuSC %s (\\mu=%.2g), %s %s', etiqueta{k}, mu, freqName, angStr));
    xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
    hold on;
    rectangle('Position',[xR(1) zR(1) diff(xR) diff(zR)], 'EdgeColor','r','LineWidth',1.5);
    rectangle('Position',[xc-R zc-R 2*R 2*R], 'Curvature',[1 1], 'EdgeColor','k','LineWidth',1.5);
    hold off;
    saveas(fig1, fullfile(figDir, sprintf('BA_DeMuSC_%s_%s_%s_mu%s.png', etiqueta{k}, freqName, angStr, muStr)));
    close(fig1);
end
    
% Recorte de las imágenes
idxStart = find(z >= 10e-3, 1, 'first');
idxEnd   = find(z <= 50e-3, 1, 'last');
z_crop_bmode = z(idxStart:idxEnd);

P_ref_L_acumu = P_ref_L_acumu(idxStart:idxEnd, :);
P_ref_H_acumu = P_ref_H_acumu(idxStart:idxEnd, :);
P_sam_L_acumu = P_sam_L_acumu(idxStart:idxEnd, :);
P_sam_H_acumu = P_sam_H_acumu(idxStart:idxEnd, :);
 
refMax = max(P_ref_L_acumu(:));
P_ref_L_Bmode = 20*log10(P_ref_L_acumu / refMax);
P_ref_H_Bmode = 20*log10(P_ref_H_acumu / refMax);
P_sam_L_Bmode = 20*log10(P_sam_L_acumu / refMax);
P_sam_H_Bmode = 20*log10(P_sam_H_acumu / refMax);
 
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

% Prueba: cuántas iteraciones usar en el sin TVR
AtA = prob.A' * prob.A;
Aty = prob.A' * prob.y;
filas = prob.z >= 10e-3 & prob.z <= 50e-3;
iters = [7 15 30 50 100 200];

for k = 1:numel(iters)
    [raw,~,relRes] = cgs(AtA, Aty, 1e-5, iters(k), [], [], 6*ones(size(prob.A,2),1));
    mapa = reshape(raw, prob.size);
    mapa = mapa(filas,:);

    vInc = mapa(maskInc);  vBg = mapa(maskBg);
    CNR = abs(mean(vInc)-mean(vBg)) / sqrt(std(vInc)^2 + std(vBg)^2);
    fprintf('iter %3d | Inc %.2f +- %.2f | Fondo %.2f +- %.2f | CNR %.2f | relRes %.3f\n', ...
        iters(k), mean(vInc), std(vInc), mean(vBg), std(vBg), CNR, relRes);

    fig = figure('Visible','off');
    imagesc(x_crop*1000, z_crop*1000, mapa); axis image; colormap("pink"); colorbar; clim([5 12]);
    title(sprintf('DeMuSC sinTVR, %d iter', iters(k)));
    saveas(fig, fullfile(figDir, sprintf('BA_DeMuSC_sinTVR_iter%d.png', iters(k))));
    close(fig);
end
 
 
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
