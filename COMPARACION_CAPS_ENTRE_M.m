% Diferencia entre m=2 y m=4 (mismo mapa de densidad)
basedir = fullfile(pwd,'CAPS-FRECUENCIA-STEERING');
probe = 'L14-5u';
tipo = 'Inc';             % 'Inc' o 'Bgnd'
freqStr = '5MHz';
angleStr = 'Angle_0';
i = 1;                    % frame

s2 = ['Simulation_' tipo '_CAPS_FreqSteering_v6_m2'];
s4 = ['Simulation_' tipo '_CAPS_FreqSteering_v6_m4'];

%% SIN beamforming (carpeta rf)
rf2 = fullfile(basedir, s2, probe, freqStr, 'rf', angleStr);
rf4 = fullfile(basedir, s4, probe, freqStr, 'rf', angleStr);

% LP = suma de las subaperturas
LP2 = 0;
for j = 1:2
    f = dir(fullfile(rf2, sprintf('%s_f%d_sub%d_of_m2_*kPa.mat', s2, i, j)));
    d = load(fullfile(rf2, f(1).name));
    LP2 = LP2 + d.rfPrebf;
end
LP4 = 0;
for j = 1:4
    f = dir(fullfile(rf4, sprintf('%s_f%d_sub%d_of_m4_*kPa.mat', s4, i, j)));
    d = load(fullfile(rf4, f(1).name));
    LP4 = LP4 + d.rfPrebf;
end

% HP = apertura completa
f = dir(fullfile(rf2, sprintf('%s_f%d_full_m2_*kPa.mat', s2, i)));
d = load(fullfile(rf2, f(1).name));
HP2 = d.rfPrebf;
f = dir(fullfile(rf4, sprintf('%s_f%d_full_m4_*kPa.mat', s4, i)));
d = load(fullfile(rf4, f(1).name));
HP4 = d.rfPrebf;

fprintf('--- SIN beamforming ---\n');
mostrarDif('LP', LP2, LP4);
mostrarDif('HP', HP2, HP4);

%% CON beamforming (carpeta bf)
bf2 = fullfile(basedir, s2, probe, freqStr, 'bf', angleStr);
bf4 = fullfile(basedir, s4, probe, freqStr, 'bf', angleStr);

a = load(fullfile(bf2, sprintf('%s_f%d_LP.mat', s2, i)));
b = load(fullfile(bf4, sprintf('%s_f%d_LP.mat', s4, i)));
bfLP2 = a.rfLp;  bfLP4 = b.rfLp;

a = load(fullfile(bf2, sprintf('%s_f%d_HP.mat', s2, i)));
b = load(fullfile(bf4, sprintf('%s_f%d_HP.mat', s4, i)));
bfHP2 = a.rfHp;  bfHP4 = b.rfHp;

fprintf('--- CON beamforming ---\n');
mostrarDif('LP', bfLP2, bfLP4);
mostrarDif('HP', bfHP2, bfHP4);

% Mapa de la diferencia del LP beamformeado (en dB respecto al maximo)
D = bfLP2 - bfLP4;
figure;
imagesc(20*log10(abs(D)/max(abs(bfLP4(:)))));
colorbar; clim([-60 0]); title('Diferencia LP m2 - m4 (dB)');

function mostrarDif(nombre, A, B)
    A = double(A);  B = double(B);
    D = A - B;
    nrmse = norm(D,'fro') / norm(B,'fro') * 100;
    c = corr(A(:), B(:));
    fprintf('%s: diferencia = %.2f %% | correlacion = %.4f\n', nombre, nrmse, c);
end