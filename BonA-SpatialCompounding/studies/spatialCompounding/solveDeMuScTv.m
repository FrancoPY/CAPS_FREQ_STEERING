function [x,diag] = solveDeMuScTv(inputImage,systemOperator,mu,tau,maxIter,nRow,nCol,initMode)

    % Set the inherited primal-dual relaxation and step sizes
    rho = 1.99;
    sigma = 0.99 / (8 * tau);
    H = nRow;
    W = nCol;
    maxIterations = 1000;
    tolerance = 1e-6;
    stableIterations = 5;

    % Prepare normal equations and the proximal system
    adjointOp   = systemOperator';
    gramMatrix  = adjointOp * systemOperator;
    projectedY  = adjointOp * inputImage;
    nPixels     = size(systemOperator, 2);
    systemMatrix   = speye(nPixels) + tau * gramMatrix;
    pcgTol      = 1e-6;
    pcgMaxIter  = 50;

    % Initialize from six using the selected historical Krylov method
    if strcmp(initMode,'pcg')
        [xInitVec,initFlag,initRes,initIter,initHist] = pcg(gramMatrix, projectedY(:), [], maxIter, [], [], 6 * ones(nPixels,1));
    else
        [xInitVec,initFlag,initRes,initIter,initHist] = cgs(gramMatrix, projectedY(:), [], maxIter, [], [], 6 * ones(nPixels,1));
    end
    xPrev = reshape(xInitVec, H, W);
    uPrev = zeros(H, W, 2);

    % Initialize objective and relative-change diagnostics
    fidelityHistory       = zeros(1, maxIterations);
    regularizationHistory = zeros(1, maxIterations);
    costHistory           = zeros(1, maxIterations);
    errorHistory          = zeros(1, maxIterations);
    proxFlag = zeros(1,maxIterations-1);
    proxRes = zeros(1,maxIterations-1);
    proxIter = zeros(1,maxIterations-1);
    fidelityHistory(1)       = 0.5 * sum((systemOperator * xPrev(:) - inputImage(:)).^2);
    regularizationHistory(1) = mu * computeTVNorm(opD(xPrev));
    costHistory(1)           = fidelityHistory(1) + regularizationHistory(1);
    errorHistory(1)          = 1;
    iteration = 1;
    costDiff  = 1;

    % Update the primal and dual variables with the original stopping rule
    while iteration < maxIterations && costDiff > tolerance
        [xCurr,proxFlag(iteration),proxRes(iteration),proxIter(iteration)] = proxTauFGeneric(xPrev - tau * opDAdj(uPrev), inputImage, tau,  adjointOp, systemMatrix, pcgTol, pcgMaxIter,H,W);
        uCurr = proxSigmaGConj(uPrev + sigma * opD(2 * xCurr - xPrev), mu);
        xPrev = xPrev + rho * (xCurr - xPrev);
        uPrev = uPrev + rho * (uCurr - uPrev);
        fidelityHistory(iteration+1)       = 0.5 * sum((systemOperator * xPrev(:) - inputImage(:)).^2);
        regularizationHistory(iteration+1) = mu * computeTVNorm(opD(xPrev));
        costHistory(iteration+1)           = fidelityHistory(iteration+1) + regularizationHistory(iteration+1);
        costDiff = abs(costHistory(iteration+1) - costHistory(iteration)) / costHistory(iteration+1);
        if iteration < stableIterations
            costDiff = 1;
        end
        errorHistory(iteration+1) = costDiff;
        iteration = iteration + 1;
    end

    % Return the numerical map and convergence history
    x = xPrev;
    diag.cost = costHistory(1:iteration);
    diag.relChange = errorHistory(1:iteration);
    diag.fidelity = fidelityHistory(1:iteration);
    diag.tvCost = regularizationHistory(1:iteration);
    diag.nIter = iteration;
    diag.converged = costDiff <= tolerance;
    diag.init = struct('method',initMode,'image',reshape(xInitVec,H,W),'flag',initFlag,'relRes',initRes,'nIter',initIter,'resHist',initHist);
    diag.prox = struct('flag',proxFlag(1:iteration-1),'relRes',proxRes(1:iteration-1),'nIter',proxIter(1:iteration-1),'tol',pcgTol,'maxIter',pcgMaxIter);
    diag.dual = uPrev;
    diag.rho = rho;
    diag.sigma = sigma;
    diag.stopTol = tolerance;
    diag.stopRule = 'relative cost change';
    diag.nUpdates = iteration-1;
    diag.hitCap = iteration==maxIterations && ~diag.converged;
end
function [x2d,flag,relRes,nIter] = proxTauFGeneric(v,inputImage,tau,adjointOp,systemMatrix,pcgTol,pcgMaxIter,H,W)
    b = v(:) + tau * adjointOp * inputImage(:);
    [xSol,flag,relRes,nIter] = pcg(systemMatrix, b, pcgTol, pcgMaxIter);
    x2d = reshape(xSol, [H, W]);
end
function u = proxSigmaGConj(v, lambda)
    normV = sqrt(sum(v.^2, 3));
    u = v ./ max(normV / lambda, 1);
end
function Dv = opD(v)
    [H, W] = size(v);
    u1 = [diff(v,1,1); zeros(1, W)];
    u2 = [diff(v,1,2), zeros(H, 1)];
    Dv = cat(3, u1, u2);
end
function v = opDAdj(u)
    div1 = -[u(1,:,1); diff(u(:,:,1),1,1)];
    div2 = -[u(:,1,2), diff(u(:,:,2),1,2)];
    v = div1 + div2;
end
function tv = computeTVNorm(Dv)
    tv = sum(sqrt(sum(Dv.^2,3)),'all');
end
