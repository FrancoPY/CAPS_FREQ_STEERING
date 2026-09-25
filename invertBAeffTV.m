function BAloc = invertBAeffTV(BAeff, mu, param)

    % invertBAeffTV  Invert cumulative B/A map into a local B/A map
    %
    % Inputs
    %   BAeff : struct with fields image, lateral, axial
    %   mu    : TV regularization weight
    %   param : struct with inversion settings
    %
    % Output
    %   BAloc : struct with local B/A result and diagnostics

    if nargin < 3
        param = struct;
    end

    % Params
    tau = getOption(param, 'tau', 0.01);
    gridSize = getOption(param, 'gridSize', 0.3e-3);
    zCrop = getOption(param, 'zCrop', [BAeff.axial(1), BAeff.axial(end)]);
    zInv = getOption(param, 'zInv', zCrop);
    maxIter = getOption(param, 'maxIter', 2000);
    tol = getOption(param, 'tol', 1e-7);
    stableIter = getOption(param, 'stableIter', 5);
    plotFlag = getOption(param, 'plotFlag', false);

    % Input
    x = BAeff.lateral(:).';
    z = BAeff.axial(:);
    BAimage = BAeff.image;

    % Axial crop
    zMask = z >= zCrop(1) & z <= zCrop(2);
    z = z(zMask);
    BAimage = BAimage(zMask, :);

    % Uniform grid
    [X, Z] = meshgrid(x, z);
    xNew = x(1):gridSize:x(end);
    zNew = z(1):gridSize:z(end);
    [Xq, Zq] = meshgrid(xNew, zNew);

    BAnew = interp2(X, Z, BAimage, Xq, Zq, 'linear');

    nanMask = isnan(BAnew);
    if any(nanMask, 'all')
        BAnear = interp2(X, Z, BAimage, Xq, Zq, 'nearest');
        BAnew(nanMask) = BAnear(nanMask);
    end

    % Inversion window
    [~, id0] = min(abs(zNew - zInv(1)));
    [~, idH] = min(abs(zNew - zInv(2)));

    if idH <= id0 + 1
        error('zInv must contain at least two depth samples after the baseline.');
    end

    z0 = zNew(id0);
    y0 = BAnew(id0, :);

    zLoc = zNew(id0+1:idH).';
    BAeffCrop = BAnew(id0+1:idH, :);

    % Residual effective map after removing baseline contribution
    Y = BAeffCrop - (z0 ./ zLoc) * y0;

    % Forward operator
    Az = buildResidualCumulativeOperator(zLoc, z0);

    % TV inversion
    [xTV, cost, errorCost, fide, regul] = totalVariationBA( ...
        Y, Az, mu, tau, maxIter, tol, stableIter);

    % Reference inversions
    xNaive = Az \ Y;

    AzT = Az';
    AzTAz = AzT * Az;
    AzTY = AzT * Y;
    xCGS = zeros(size(Y));
    xCGS1 = zeros(size(Y));
    xCGS2 = zeros(size(Y));
    
    for iCol = 1:size(Y, 2)
        xCGS1(:, iCol) = cgs(AzTAz, AzTY(:, iCol), 1e-6, 1);
        xCGS2(:, iCol) = cgs(AzTAz, AzTY(:, iCol), 1e-6, 2);
        xCGS(:, iCol) = cgs(AzTAz, AzTY(:, iCol), 1e-6, 200);
    end

    % Output
    BAloc.image = xTV;
    BAloc.noisyimage = xCGS;
    BAloc.naiveimage = xNaive;
    BAloc.effimage = BAeffCrop;
    BAloc.residEffimage = Y;
    BAloc.lateral = xNew;
    BAloc.axial = zLoc;
    BAloc.baselineDepth = z0;
    BAloc.baselineValue = y0;
    BAloc.cost = cost;
    BAloc.error = errorCost;
    BAloc.fide = fide;
    BAloc.regul = regul;
    BAloc.cgsIter1image = xCGS1;
    BAloc.cgsIter2image = xCGS2;
    BAloc.cgsDelta12image = abs(xCGS2 - xCGS1);
    BAloc.cgsWeight12image = buildCGSWeightMap(xCGS1, xCGS2);

    if plotFlag
        showInversionPlots(BAloc, mu, tau);
    end

end

function Az = buildResidualCumulativeOperator(zLoc, z0)

    % buildResidualCumulativeOperator  Residual cumulative operator

    dz = diff([z0; zLoc]);
    n = numel(zLoc);

    Az = tril(ones(n, n));
    Az = Az .* (ones(n, 1) * dz.');
    Az = Az ./ (zLoc * ones(1, n));

end

function [x, cost, errorCost, fide, regul] = totalVariationBA(Y, Az, lambda, tau, maxIter, tol, stableIter)

    % totalVariationBA  TV inversion for local B/A image

    rho = 1.99;
    sigma = 1 / tau / 8;

    [H, W] = size(Y);
    At = Az';
    AtA = At * Az;
    AtY = At * Y;

    R = chol(speye(H) + tau * AtA);

    x2 = Y;
    u2 = zeros(H, W, 2);

    cost = zeros(1, maxIter);
    errorCost = zeros(1, maxIter);
    fide = zeros(1, maxIter);
    regul = zeros(1, maxIter);

    cost(1) = inf;
    errorCost(1) = 1;

    ee = 1;
    iter = 1;

    while iter < maxIter && ee > tol

        x = proxTauF(x2 - tau * opDadj(u2), R, AtY, tau);
        u = proxSigmaGConj(u2 + sigma * opD(2 * x - x2), lambda);

        x2 = x2 + rho * (x - x2);
        u2 = u2 + rho * (u - u2);

        aux = Az * x - Y;
        fide(iter + 1) = 0.5 * sum(aux(:).^2);
        regul(iter + 1) = lambda * TVV(opD(x));
        cost(iter + 1) = fide(iter + 1) + regul(iter + 1);

        if isinf(cost(iter))
            ee = 1;
        else
            ee = abs(cost(iter + 1) - cost(iter)) / max(cost(iter + 1), eps);
        end

        if iter < stableIter
            ee = 1;
        end

        errorCost(iter + 1) = ee;
        iter = iter + 1;
    end

    cost = cost(1:iter);
    errorCost = errorCost(1:iter);
    fide = fide(1:iter);
    regul = regul(1:iter);

end

function U = proxTauF(V, R, AtY, tau)

    rhs = V + tau * AtY;
    U = R \ (R' \ rhs);

end

function U = proxSigmaGConj(V, lambda)

    mag = sqrt(sum(V.^2, 3));
    U = V ./ max(mag / lambda, 1);

end

function U = opD(V)

    [H, W] = size(V);
    U = cat(3, [diff(V, 1, 1); zeros(1, W)], [diff(V, 1, 2), zeros(H, 1)]);

end

function U = opDadj(V)

    U = -[V(1, :, 1); diff(V(:, :, 1), 1, 1)] ...
        -[V(:, 1, 2), diff(V(:, :, 2), 1, 2)];

end

function val = TVV(DV)

    val = sum(sum(sqrt(sum(DV.^2, 3))));

end

function showInversionPlots(BAloc, mu, tau)

    font = 18;
    lw = 2;
    bAxis = [0, 12];

    figure;
    subplot(1, 2, 1);
    plot(BAloc.cost, 'LineWidth', lw);
    grid on;
    xlabel('Iteration');
    ylabel('Cost');
    title(sprintf('\\mu = %.3g', mu));
    set(gca, 'FontSize', font);

    subplot(1, 2, 2);
    plot(BAloc.error, 'LineWidth', lw);
    grid on;
    xlabel('Iteration');
    ylabel('Relative cost change');
    title(sprintf('\\tau = %.3g', tau));
    set(gca, 'FontSize', font);

    figure;
    subplot(2, 2, 1);
    imagesc(1e3 * BAloc.lateral, 1e3 * BAloc.axial, BAloc.residEffimage, bAxis);
    axis image;
    colormap turbo;
    colorbar;
    clim([5, 12]);
    title('Residual effective');
    set(gca, 'FontSize', font);

    subplot(2, 2, 2);
    imagesc(1e3 * BAloc.lateral, 1e3 * BAloc.axial, BAloc.naiveimage, bAxis);
    axis image;
    colormap turbo;
    colorbar;
    clim([5, 12]);
    title('Local Naive');
    set(gca, 'FontSize', font);

    subplot(2, 2, 3);
    imagesc(1e3 * BAloc.lateral, 1e3 * BAloc.axial, BAloc.noisyimage, bAxis);
    axis image;
    colormap turbo;
    colorbar;
    clim([5, 12]);
    title('Local CGS');
    set(gca, 'FontSize', font);

    subplot(2, 2, 4);
    imagesc(1e3 * BAloc.lateral, 1e3 * BAloc.axial, BAloc.image, bAxis);
    axis image;
    colormap turbo;
    colorbar;
    clim([5, 12]);
    title('Local TV');
    set(gca, 'FontSize', font);

end

function val = getOption(s, fieldName, defaultVal)

    if isfield(s, fieldName) && ~isempty(s.(fieldName))
        val = s.(fieldName);
    else
        val = defaultVal;
    end

end

function W = buildCGSWeightMap(x1, x2)

    % buildCGSWeightMap  Normalized spatial weights from early CGS change

    D = abs(x2 - x1);
    D = D / prctile(D(:), 95);
    D = min(D, 1);

    W = 1 - D;
    W = max(W, 0);

end