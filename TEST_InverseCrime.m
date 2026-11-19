% TEST_InverseCrime.m
% Pasos 1 a 3: mapa verdadero -> operador directo (A) -> inversión sin ruido

param.gridSize = 0.3e-3;
param.zMin     = 10e-3;
param.zMax     = 50e-3;

paso   = 2.5;                               % 5, 2.5, 1 o 0.5
angNum = [0, setdiff(-15:paso:15, 0)];
nAng   = numel(angNum);

figDir = fullfile(pwd, 'figuras_InverseCrime');
if ~exist(figDir, 'dir'); mkdir(figDir); end

%% Paso 1: ejes y mapa verdadero
x = (0:127)*0.3e-3;  x = x - mean(x);       % 128 elementos, pitch 0.3 mm
z = (0:0.1e-3:65e-3)';

% B/A de relleno, solo para armar la matriz A (y se reemplaza más abajo)
ba   = ones(numel(z), numel(x), nAng);
prob = prepDeMuSc(ba, x, z, angNum, param);

[XX, ZZ] = meshgrid(prob.x*1000, prob.z*1000);   % en mm
zc = 22.5;  xc = 0;  R = 7.5;
dist  = sqrt((XX-xc).^2 + (ZZ-zc).^2);
uVerd = 6*ones(prob.size);                       % fondo
uVerd(dist <= R) = 11;                           % inclusión

%% Paso 2: operador directo (inverse crime: es el mismo A de la inversión)
y = prob.A * uVerd(:);

%% ROIs (iguales a los de los BA_COMPUTE)
filas = prob.z >= 10e-3 & prob.z <= 50e-3;
Xc = XX(filas,:);  Zc = ZZ(filas,:);
distC = sqrt((Xc-xc).^2 + (Zc-zc).^2);
xR = [-18 18];  zR = [12 47];
maskInc = distC <= R;
maskBg  = Xc>=xR(1) & Xc<=xR(2) & Zc>=zR(1) & Zc<=zR(2) & distC > R;
uC = uVerd(filas,:);
xmm = prob.x*1000;
zmm = prob.z(filas)*1000;

%% Figura: mapa verdadero y B/A_C que sale del operador (ángulo 0 y último)
ini = 0;  yAng = cell(1, nAng);
for a = 1:nAng
    n = nnz(prob.mask{a});
    m = nan(prob.size);
    m(prob.mask{a}) = y(ini+1:ini+n);
    yAng{a} = m(filas,:);
    ini = ini + n;
end

fig = figure('Visible','off');
subplot(1,3,1); imagesc(xmm, zmm, uC); axis image; colormap("pink"); colorbar; clim([5 12]);
title('B/A local verdadero');
subplot(1,3,2); imagesc(xmm, zmm, yAng{1}, 'AlphaData', ~isnan(yAng{1})); axis image; colorbar; clim([5 12]);
title(sprintf('B/A_C ángulo %g°', angNum(1)));
subplot(1,3,3); imagesc(xmm, zmm, yAng{end}, 'AlphaData', ~isnan(yAng{end})); axis image; colorbar; clim([5 12]);
title(sprintf('B/A_C ángulo %g°', angNum(end)));
saveas(fig, fullfile(figDir, sprintf('Verdadero_y_acumulado_%dang.png', nAng)));
close(fig);

%% Paso 3: inversión sin ruido (el mismo cgs del "sin TVR")
AtA = prob.A' * prob.A;
Aty = prob.A' * y;
iters = [2 5 7 15 50 100 200];
mapas = cell(1, numel(iters));

fprintf('\n===== InverseCrime sin ruido | %d ángulos | grid %.1f mm =====\n', nAng, param.gridSize*1e3);
for k = 1:numel(iters)
    [raw, ~] = cgs(AtA, Aty, 1e-5, iters(k), [], [], 6*ones(size(prob.A,2),1));
    mapa = reshape(raw, prob.size);
    mapa = mapa(filas,:);
    mapas{k} = mapa;

    imprimirStats(sprintf('iter %d', iters(k)), mapa, uC, maskInc, maskBg);

    fig = figure('Visible','off');
    imagesc(xmm, zmm, mapa); axis image; colormap("pink"); colorbar; clim([5 12]);
    title(sprintf('Sin ruido, %d iter, %d ángulos', iters(k), nAng));
    xlabel('Posición Lateral (mm)'); ylabel('Profundidad (mm)');
    saveas(fig, fullfile(figDir, sprintf('SinRuido_%dang_iter%d.png', nAng, iters(k))));
    close(fig);
end

% Perfil en profundidad por el centro de la inclusión
fig = figure('Visible','off'); hold on;
centro = abs(xmm) < 2;
plot(zmm, mean(uC(:,centro), 2), 'k', 'LineWidth', 1.5, 'DisplayName', 'verdadero');
for k = [1 3 numel(iters)]
    plot(zmm, mean(mapas{k}(:,centro), 2), 'DisplayName', sprintf('iter %d', iters(k)));
end
xlabel('Profundidad (mm)'); ylabel('B/A'); legend; grid on;
saveas(fig, fullfile(figDir, sprintf('Perfil_SinRuido_%dang.png', nAng)));
close(fig);

function imprimirStats(nombre, mapa, uVerd, maskInc, maskBg)
    vInc = mapa(maskInc);  vBg = mapa(maskBg);
    muInc = mean(vInc);  sdInc = std(vInc);
    muBg  = mean(vBg);   sdBg  = std(vBg);
    CNR   = abs(muInc - muBg) / sqrt(sdInc^2 + sdBg^2);
    roi   = maskInc | maskBg;
    rmse  = sqrt(mean((mapa(roi) - uVerd(roi)).^2));
    fprintf('%-10s | Inc %.2f +- %.2f | Fondo %.2f +- %.2f | Dif %.2f | CNR %.2f | RMSE %.2f\n', ...
        nombre, muInc, sdInc, muBg, sdBg, muInc-muBg, CNR, rmse);
end