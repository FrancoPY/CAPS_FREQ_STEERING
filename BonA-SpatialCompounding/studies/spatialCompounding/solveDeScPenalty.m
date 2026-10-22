function [u,diag,uInit] = solveDeScPenalty(A,b,mask,mu,cfg,uInit)
    % Solve matching DeSC penalties; preserve gradient-coupled padding.

    % Assemble the inherited pixel differences and fixed data terms
    [m,n] = size(mask);
    dRow = spdiags([-ones(m,1),ones(m,1)],[0,1],m,m+1);
    dRow(:,end) = [];
    dRow(m,m) = 0;
    dCol = spdiags([-ones(n,1),ones(n,1)],[0,1],n,n+1);
    dCol(:,end) = [];
    dCol(n,n) = 0;
    D = [kron(speye(n),dRow);kron(dCol,speye(m))];
    gradMask = [mask(:);mask(:)];
    ata = A'*A;
    rhs = A'*b;
    quad = makeQuad(32);
    dataIdx = find(sum(abs(A),1)' > 0);
    gradOp = spdiags(gradMask,0,2*m*n,2*m*n)*D;
    solveIdx = find(sum(abs(A),1)'+sum(abs(gradOp),1)' > 0);
    nullIdx = setdiff((1:m*n)',solveIdx);

    % Use a data-only initialization shared by all penalties
    initDiag = struct('flag',NaN,'relres',NaN,'nIter',NaN,'provided',~isempty(uInit));
    if isempty(uInit)
        u = zeros(m*n,1);
        if strcmp(cfg.innerSolve,'direct')
            u(dataIdx) = A(:,dataIdx)\b;
            initDiag.flag = 0;
        else
            [u(dataIdx),initDiag.flag,initDiag.relres,initDiag.nIter] = pcg(ata(dataIdx,dataIdx),rhs(dataIdx),cfg.pcgTol,cfg.maxPcg);
        end
    else
        u = uInit(:);
    end
    uInit = u;
    nullInit = u(nullIdx);
    [own,report,g] = evalCost(u,A,b,D,mask,mu,cfg,quad);
    normRhs = max(norm(rhs(solveIdx),inf),1);
    hist.own = own;
    hist.report = report;
    hist.stationarity = norm(g(solveIdx),inf)/normRhs;
    hist.ownDelta = NaN;
    hist.reportDelta = NaN;
    hist.dataResidual = norm(A*u-b);
    hist.innerFlag = NaN;
    hist.innerRelres = NaN;
    hist.innerIter = NaN;
    hist.systemResidual = NaN;
    snap = struct('report',[],'own',[],'station',[]);
    iIter = 0;
    done = false;

    % Follow one trajectory and retain the first crossing of each stopping rule
    while iIter < cfg.maxOuter && ~done
        iIter = iIter+1;
        r = hypot(D(1:m*n,:)*u,D(m*n+1:end,:)*u);
        [~,w] = evalPenalty(r,cfg,quad);
        wVec = [w;w].*gradMask;
        H = ata+mu*(D'*spdiags(wVec,0,2*m*n,2*m*n)*D);
        if strcmp(cfg.innerSolve,'direct')
            u(solveIdx) = H(solveIdx,solveIdx)\rhs(solveIdx);
            flag = 0;
            relres = norm(H(solveIdx,:)*u-rhs(solveIdx))/max(norm(rhs(solveIdx)),eps);
            nInner = 0;
        else
            [u(solveIdx),flag,relres,nInner] = pcg(H(solveIdx,solveIdx),rhs(solveIdx),cfg.pcgTol,cfg.maxPcg,[],[],u(solveIdx));
        end
        [own,report,g] = evalCost(u,A,b,D,mask,mu,cfg,quad);
        ownDelta = abs(own-hist.own(end));
        reportDelta = abs(report-hist.report(end));
        station = norm(g(solveIdx),inf)/normRhs;
        hist.own(end+1,1) = own;
        hist.report(end+1,1) = report;
        hist.stationarity(end+1,1) = station;
        hist.ownDelta(end+1,1) = ownDelta;
        hist.reportDelta(end+1,1) = reportDelta;
        hist.dataResidual(end+1,1) = norm(A*u-b);
        hist.innerFlag(end+1,1) = flag;
        hist.innerRelres(end+1,1) = relres;
        hist.innerIter(end+1,1) = nInner;
        hist.systemResidual(end+1,1) = norm(H(solveIdx,:)*u-rhs(solveIdx),inf)/normRhs;
        point = struct('u',u,'nIter',iIter,'own',own,'report',report,'stationarity',station,'ownDelta',ownDelta,'reportDelta',reportDelta);
        if isempty(snap.report) && reportDelta <= cfg.objTol
            snap.report = point;
        end
        if isempty(snap.own) && ownDelta <= cfg.objTol
            snap.own = point;
        end
        if isempty(snap.station) && station <= cfg.stationTol
            snap.station = point;
        end
        done = ~isempty(snap.station);
        if strcmp(cfg.kind,'inherit')
            done = done && ~isempty(snap.report) && ~isempty(snap.own);
        end
    end
    diag = struct('kind',cfg.kind,'innerSolve',cfg.innerSolve,'nIter',iIter,'converged',~isempty(snap.station),'capHit',~done,'init',initDiag,'history',hist,'snap',snap,'nNull',numel(nullIdx),'nullChange',norm(u(nullIdx)-nullInit,inf),'nInnerFail',nnz(hist.innerFlag(2:end)~=0),'maxInnerIter',max(hist.innerIter(2:end)),'nOwnIncrease',nnz(diff(hist.own)>1e-10));
    if ~isempty(snap.station)
        u = snap.station.u;
    end

end

function [own,report,g] = evalCost(u,A,b,D,mask,mu,cfg,quad)

    % Evaluate the chosen penalty and the inherited reported TV on the same field
    nPix = numel(mask);
    du = D*u;
    r = hypot(du(1:nPix),du(nPix+1:end));
    [rho,w] = evalPenalty(r,cfg,quad);
    res = A*u-b;
    fidelity = 0.5*(res'*res);
    own = fidelity+mu*sum(mask(:).*rho);
    report = fidelity+mu*sum(mask(:).*r);
    wVec = [w;w].*[mask(:);mask(:)];
    g = A'*res+mu*D'*(wVec.*du);

end

function [rho,w] = evalPenalty(r,cfg,quad)

    % Match each IRLS weight exactly to its radial penalty
    switch cfg.kind
        case 'inherit'
            vMax = asinh(r.^2/cfg.epsVal);
            rho = sqrt(cfg.epsVal)/(2*sqrt(2))*vMax.*sum(sqrt(cosh(vMax*quad.nodes)).*quad.weights,2);
            w = (r.^4+cfg.epsVal^2).^(-0.25)/sqrt(2);
        case 'pseudoHuber'
            rho = (sqrt(r.^2+cfg.epsVal)-sqrt(cfg.epsVal))/sqrt(2);
            w = 1./sqrt(2*(r.^2+cfg.epsVal));
        case 'smoothedTV'
            rho = sqrt(r.^2+cfg.epsVal)-sqrt(cfg.epsVal);
            w = 1./sqrt(r.^2+cfg.epsVal);
        otherwise
            error('Unknown trial penalty.');
    end

end

function quad = makeQuad(nQuad)

    % Integrate the inherited penalty with the verified sinh-Gauss rule
    k = (1:nQuad-1)';
    beta = k./sqrt(4*k.^2-1);
    [vec,val] = eig(diag(beta,1)+diag(beta,-1));
    [nodes,idx] = sort(diag(val));
    quad.nodes = (nodes+1)'/2;
    quad.weights = vec(1,idx).^2;

end
