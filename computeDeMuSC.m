function [BA, info] = computeDeMuSC(BAeffCell, anglesDeg, zGrid, xGrid, mu, opts)
% BA local DeMuSC a partir de los mapas de B/A efectivo de cada angulo
if nargin == 0; selfTest(); return; end
if nargin < 6 || isempty(opts); opts = struct; end
if isfield(opts, 'weightMode'); wm = opts.weightMode; else; wm = 'surface'; end

zGrid = zGrid(:); xGrid = xGrid(:).';
m = numel(zGrid); n = numel(xGrid);

Wc = cell(1, numel(anglesDeg));
Yc = Wc;
for k = 1:numel(anglesDeg)
    [Wk, rinfo] = buildRayOperator(zGrid, xGrid, anglesDeg(k), wm);
    img = BAeffCell{k};
    y = img(rinfo.rowPix);
    yb = interp1(xGrid, img(1,:), rinfo.xTop, 'linear', 'extrap');
    Y = y - rinfo.baseRatio .* yb(:);
    good = ~isnan(Y);
    Wc{k} = Wk(good,:);
    Yc{k} = Y(good);
end

[BA, info] = solveDeMuSC(vertcat(Wc{:}), vertcat(Yc{:}), m, n, mu, opts);
info.nRows = sum(cellfun(@numel, Yc));
end

function [W, info] = buildRayOperator(zGrid, xGrid, thetaDeg, weightMode)
% matriz de rayos W para un angulo
if nargin < 4 || isempty(weightMode); weightMode = 'surface'; end
zGrid = zGrid(:); xGrid = xGrid(:).';
m = numel(zGrid); n = numel(xGrid);
dz = zGrid(2) - zGrid(1);
dx = xGrid(2) - xGrid(1);
th = thetaDeg*pi/180;
sT = sin(th);
cT = cos(th);
useSurface = ~strcmpi(weightMode, 'domain');

zTop = zGrid(1) - dz/2;
xMin = xGrid(1) - dx/2;
xMax = xGrid(end) + dx/2;
tol = 1e-9*dz;

rowsC = cell(m,1);
colsC = cell(m,1);
valsC = cell(m,1);
xTopFull = nan(m*n,1);
baseFull = zeros(m*n,1);
validFull = false(m*n,1);

for i = 1:m
    zi = zGrid(i);
    tEnd = (zi - zTop)/cT;

    % cruces del rayo con las lineas de la grilla
    tH = ((0:i-2) + 0.5)*dz/cT;
    if abs(sT) > 1e-12
        K = floor(tEnd*abs(sT)/dx - 0.5);
        tV = ((0:max(K,-1)) + 0.5)*dx/abs(sT);
        tV = tV(tV < tEnd - tol);
    else
        tV = [];
    end
    tAll = sort([0, tH, tV, tEnd]);
    L = diff(tAll);
    tm = 0.5*(tAll(1:end-1) + tAll(2:end));
    keep = L > tol;
    L = L(keep);
    tm = tm(keep);

    rIdx = floor((zi - tm*cT - zTop)/dz) + 1;
    dk = round(-tm*sT/dx);

    if useSurface; w = L/(zi/cT); else; w = L/tEnd; end

    % rayos que no se salen por los costados
    xEnd = xGrid - tEnd*sT;
    okj = find(xEnd >= xMin - tol & xEnd <= xMax + tol);
    if isempty(okj); continue; end

    S = numel(L);
    pix = i + m*(okj - 1);
    colMat = okj(:) + dk(:).';
    cellId = rIdx(:).' + m*(colMat - 1);
    rowM = repmat(pix(:), 1, S);
    valM = repmat(w(:).', numel(okj), 1);
    rowsC{i} = rowM(:);
    colsC{i} = cellId(:);
    valsC{i} = valM(:);

    validFull(pix) = true;
    xTopFull(pix) = xEnd(okj);
    if useSurface; baseFull(pix) = zTop/zi; end
end

Wfull = sparse(vertcat(rowsC{:}), vertcat(colsC{:}), vertcat(valsC{:}), m*n, m*n);
sel = find(validFull);
W = Wfull(sel, :);

info.rowPix = sel;
info.xTop = xTopFull(sel);
info.zTop = zTop;
info.baseRatio = baseFull(sel);
info.theta = thetaDeg;
info.nRays = numel(sel);
info.nTotal = m*n;
end

function [xMap, info] = solveDeMuSC(W, Y, m, n, mu, opts)
% minimiza 0.5*||W x - Y||^2 + mu*TV(x)
if nargin < 6 || isempty(opts); opts = struct; end
maxOuter = getOpt(opts, 'maxOuter', 25);
tolOuter = getOpt(opts, 'tolOuter', 1e-4);
maxPCG = getOpt(opts, 'maxPCG', 300);
tolPCG = getOpt(opts, 'tolPCG', 1e-6);
epsTV = getOpt(opts, 'epsTV', 1e-2);
verbose = getOpt(opts, 'verbose', true);

N = m*n;
Y = Y(:);

% derivadas hacia adelante
Dm = spdiags([-ones(m,1) ones(m,1)], [0 1], m, m);
Dm(m,:) = 0;
Dn = spdiags([-ones(n,1) ones(n,1)], [0 1], n, n);
Dn(n,:) = 0;
Dz = kron(speye(n), Dm);
Dx = kron(Dn, speye(m));
Dz2 = Dz.^2;
Dx2 = Dx.^2;

Wt = W';
b = Wt*Y;
dW = full(sum(W.^2, 1))';

rowSum = full(sum(W, 2));
x = ones(N,1) * median(Y ./ max(rowSum, eps));

cost = zeros(1, maxOuter);
for it = 1:maxOuter
    gz = Dz*x;
    gx = Dx*x;
    phi = 1 ./ sqrt(gz.^2 + gx.^2 + epsTV^2);

    Afun = @(v) Wt*(W*v) + mu*(Dz'*(phi.*(Dz*v)) + Dx'*(phi.*(Dx*v)));
    Mp = spdiags(dW + mu*(Dz2'*phi + Dx2'*phi) + 1e-8, 0, N, N);

    [xn, ~] = pcg(Afun, b, tolPCG, maxPCG, Mp, [], x);

    rel = norm(xn - x) / max(norm(x), eps);
    x = xn;
    res = W*x - Y;
    cost(it) = 0.5*(res'*res) + mu*sum(sqrt((Dz*x).^2 + (Dx*x).^2));
    if verbose
        fprintf('iter %d  costo=%.4e  cambio=%.2e\n', it, cost(it), rel);
    end
    if rel < tolOuter; break; end
end

xMap = reshape(x, m, n);
info.cost = cost(1:it);
info.iters = it;
end

function v = getOpt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function selfTest()
rng(2);
dz = 0.3e-3;
zGrid = (10e-3+dz/2 : dz : 45e-3-dz/2)';
xGrid = (-63.5:63.5)*0.3e-3;
m = numel(zGrid);
n = numel(xGrid);
[XG, ZG] = meshgrid(xGrid, zGrid);
xT = 6*ones(m,n);
xT((XG.^2 + (ZG-22.5e-3).^2) <= (7.5e-3)^2) = 11;
ang = -15:5:15;
E = cell(1,7);
for k = 1:7
    [W, rinfo] = buildRayOperator(zGrid, xGrid, ang(k), 'surface');
    img = nan(m,n);
    img(rinfo.rowPix) = W*xT(:) + rinfo.baseRatio*6 + 0.3*randn(rinfo.nRays,1);
    E{k} = img;
end
o.maxOuter = 8;
o.verbose = false;
BA = computeDeMuSC(E, ang, zGrid, xGrid, 0.05, o);
inc = (XG.^2 + (ZG-22.5e-3).^2) <= (5e-3)^2;
bg = (abs(XG) > 12e-3) & (ZG > 15e-3) & (ZG < 40e-3);
fprintf('inclusion %.2f | fondo %.2f\n', mean(BA(inc)), mean(BA(bg)));
end
