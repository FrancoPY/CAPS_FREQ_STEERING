% TEST_InverseCrime_Ruido.m
% Mapa verdadero -> A -> ruido gaussiano -> inversión con invertDeMuSc (sin TVR y TV)
% Caso A: ruido en el acumulado (y + ruido)
% Caso B: ruido en el mapa local antes de A (A*(u + ruido))

rng(1);                                     % misma semilla siempre

g = 0.6e-3;
param.gridSize = g;
param.zMin     = 10e-3;
param.zMax     = 50e-3;

paso   = 2.5;                               % 5, 2.5, 1 o 0.5
angNum = [0, setdiff(-15:paso:15, 0)];
nAng   = numel(angNum);

sigmas = [0.1 0.5 1.0];
iters  = [2 5 7 15 50];
param.mu  = 10^-0.5;           % los mismos de tu BA_COMPUTE_DeMuSC
param.tau = 0.1;

figDir = fullfile(pwd, 'figuras_InverseCrime_Ruido');
if ~exist(figDir, 'dir'); mkdir(figDir); end

%% Mapa verdadero y operador A
x = (0:127)*g;  x = x - mean(x);
z = (0:0.1e-3:65e-3)';

ba   = ones(numel(z), numel(x), nAng);      % relleno, solo para armar A
prob = prepDeMuSc(ba, x, z, angNum, param);

[XX, ZZ] = meshgrid(prob.x*1000, prob.z*1000);
zc = 22.5;  xc = 0;  R = 7.5;
dist  = sqrt((XX-xc).^2 + (ZZ-zc).^2);
uVerd = 6*ones(prob.size);
uVerd(dist <= R) = 11;

y = prob.A * uVerd(:);                      % acumulado sin ruido

%% ROIs
filas = prob.z >= 10e-3 & prob.z <= 50e-3;
Xc = XX(filas,:);  Zc = ZZ(filas,:);
distC = sqrt((Xc-xc).^2 + (Zc-zc).^2);
xR = [-18 18];  zR = [12 47];
maskInc = distC <= R;
maskBg  = Xc>=xR(1) & Xc<=xR(2) & Zc>=zR(1) & Zc<=zR(2) & distC > R;
uC  = uVerd(filas,:);
xmm = prob.x*1000;
zmm = prob.z(filas)*1000;

%% Ruido y inversión (con invertDeMuSc)
casos = {'A', 'B'};
for s = 1:numel(sigmas)
    sig = sigmas(s);

    yA = y + sig*randn(size(y));                    % A: ruido en el acumulado
    uB = uVerd + sig*randn(size(uVerd));            % B: ruido en el local
    yB = prob.A * uB(:);
    yCasos = {yA, yB};

    for c = 1:2
        fprintf('\n===== Ruido %s | sigma %.1f | %d ángulos | grid %.1f mm =====\n', ...
            casos{c}, sig, nAng, g*1e3);
        prob.y = yCasos{c};                         % DeMuSC invierte este y

        fig = figure('Visible','off','Position',[100 100 1500 300]);
        for k = 1:numel(iters)
            param.maxIter = iters(k);
            res = invertDeMuSc(prob, param);

            mapa = res.noReg(filas,:);              % sin TVR
            imprimirStats(sprintf('sinTVR %d', iters(k)), mapa, uC, maskInc, maskBg);

            if iters(k) == 7                        % TV solo para iter 7 (la de tus corridas)
                imprimirStats('TV (iter 7)', res.image(filas,:), uC, maskInc, maskBg);
            end

            subplot(1, numel(iters), k);
            imagesc(xmm, zmm, mapa); axis image; colormap("pink"); clim([5 12]);
            title(sprintf('iter %d', iters(k)));
        end
        sgtitle(sprintf('Ruido %s, sigma %.1f, %d ángulos', casos{c}, sig, nAng));
        saveas(fig, fullfile(figDir, sprintf('Ruido%s_sigma%.1f_%dang.png', casos{c}, sig, nAng)));
        close(fig);
    end
end

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