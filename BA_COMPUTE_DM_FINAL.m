% Cómputo de B/A con CAPS (Complementary-Aperture Pulse Sequence).
%
% Diferencia clave respecto a la versión DM clásica: la señal "L"
% (rfLp) ya NO es una adquisición física de baja presión, sino la suma
% coherente de las m subaperturas (ver beamform_to_bf.m), que reconstruye
% una señal equivalente a "apertura completa a presión ν≈m veces menor"
% por construcción geométrica. Por eso aquí NO se multiplica por ningún
% factor de escala v: la ecuación de depleción usa directamente
% (P_L,Σ - P_H) en vez de (v*P_L - P_H), como indica la ecuación 12 del
% paper de Avilés et al. El resto del método (calibración con phantom de
% referencia, corrección por atenuación) es idéntico al DM clásico.

refSample = 'Phantom_homo_BA6_Steering15_CAPS';
samSample = 'Phantom_inc_BA11_Steering15_CAPS';

refDir = 'C:\Users\FRANCO PERALTA\Documents\MATLAB\LIM\CAPS\TAREA_CAPS\Phantom_homo_BA6_Steering15_CAPS\L14-5u\7MHz\bf';
samDir = 'C:\Users\FRANCO PERALTA\Documents\MATLAB\LIM\CAPS\TAREA_CAPS\Phantom_inc_BA11_Steering15_CAPS\L14-5u\7MHz\bf';

% Cargar la referencia
ref_L0 = load(fullfile(refDir, 'Phantom_homo_BA6_Steering15_CAPS_f1_LP.mat'));
fs = ref_L0.fs;
zAxis = ref_L0.zAxis;
xAxis = ref_L0.xAxis;
c0 = 1500;

fprintf('m usado en la referencia: %d\n', ref_L0.m);

B_r = 4;     % beta de referencia (calibrado con B/A_ref = 6)

f_fund = 7e6;
f_fund_r = 7e6;

alpha = 0.10*((f_fund/1e6)^2)*100/8.686;
alpha_r = 0.10*((f_fund_r/1e6)^2)*100/8.686;

bw = 0.8e6;
f_low = f_fund - bw/2;
f_high = f_fund + bw/2;
Wn = [f_low f_high] / (fs/2);
order = 200;
hFilter = fir1(order, Wn);

z = zAxis(:);
x = xAxis;
term = ((1-exp(-2*alpha_r*z))./(1-exp(-2*alpha*z)))*(alpha/alpha_r);

P_ref_L_sum = 0;
P_ref_H_sum = 0;
P_sam_L_sum = 0;
P_sam_H_sum = 0;

for i = 1:6
    % Cargar
    rL = load(fullfile(refDir, sprintf('%s_f%d_LP.mat', refSample, i)));
    rH = load(fullfile(refDir, sprintf('%s_f%d_HP.mat', refSample, i)));
    sL = load(fullfile(samDir, sprintf('%s_f%d_LP.mat', samSample, i)));
    sH = load(fullfile(samDir, sprintf('%s_f%d_HP.mat', samSample, i)));

    ref_L    = rL.rfLp;   % suma coherente de subaperturas (≈ ν·P_L)
    ref_H    = rH.rfHp;   % apertura completa
    sample_L = sL.rfLp;
    sample_H = sH.rfHp;

    % Filtrado
    ref_L_f = filtfilt(hFilter, 1, ref_L);
    ref_H_f = filtfilt(hFilter, 1, ref_H);
    sample_L_f = filtfilt(hFilter, 1, sample_L);
    sample_H_f = filtfilt(hFilter, 1, sample_H);

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

P_ref_L_mean = P_ref_L_sum / 6;
P_ref_H_mean = P_ref_H_sum / 6;
P_sam_L_mean = P_sam_L_sum / 6;
P_sam_H_mean = P_sam_H_sum / 6;

% Calcular B/A (CAPS: sin factor v, ya está implícito en la suma coherente)
numerador = (P_sam_L_mean - P_sam_H_mean).*(P_ref_L_mean);
denominador = (P_ref_L_mean - P_ref_H_mean).*(P_sam_L_mean);
ratio = numerador ./ denominador;
B = B_r * sqrt(abs(ratio)) .* term;
B_A = 2*(B - 1);

% Suavizado axial y lateral
lambda = c0/f_fund;
muestras_ventana = round((2*10*lambda*fs)/c0);
B_suave_axial = movmean(B_A, muestras_ventana, 1);
B_final = movmean(B_suave_axial, 4, 2);

% Visualización
figure;
imagesc(x*1000, z*1000, B_final);
axis image; colormap(turbo); colorbar;
title('Mapa B/A (CAPS, con steering)'); xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
clim([5 12]);
fila_inicio = find(z>=0.015, 1, 'first');
fila_fin = find(z<=0.030, 1, 'last');

media = mean(B_final(fila_inicio:fila_fin, :), 'all');
fprintf('La media es: %.2f\n', media);

% Visualización B-mode

P_ref_L_Bmode = 20*log10(P_ref_L_frame1 / max(P_ref_L_frame1(:)));
P_ref_H_Bmode = 20*log10(P_ref_H_frame1 / max(P_ref_H_frame1(:)));
P_sam_L_Bmode = 20*log10(P_sam_L_frame1 / max(P_sam_L_frame1(:)));
P_sam_H_Bmode = 20*log10(P_sam_H_frame1 / max(P_sam_H_frame1(:)));

dynamicRange = 60;

figure;
subplot(2,2,1); imagesc(x*1000, z*1000, P_ref_L_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title('B-mode Referencia \Sigma-subaperturas (Homo BA6)');
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,2); imagesc(x*1000, z*1000, P_ref_H_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title('B-mode Referencia Apertura completa (Homo BA6)');
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,3); imagesc(x*1000, z*1000, P_sam_L_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title('B-mode Muestra \Sigma-subaperturas (Inc BA11)');
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');

subplot(2,2,4); imagesc(x*1000, z*1000, P_sam_H_Bmode);
axis image; colormap gray; colorbar; clim([-dynamicRange 0]);
title('B-mode Muestra Apertura completa (Inc BA11)');
xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
