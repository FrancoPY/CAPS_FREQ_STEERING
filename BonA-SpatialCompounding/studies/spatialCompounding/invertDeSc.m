function result = invertDeSc(ba,x,z,angles,param,doTv)
    % Invert cumulative B/A per angle on source-relative paths; axes are meters.
    if nargin < 6
        doTv = true;
    end
    assert(angles(1)==0,'DeSC requires broadside first for the selected study convention.');

    % Crop and resample cumulative observations onto the declared physical grid
    [~,idL] = min(abs(z));
    [~,idH] = min(abs(z-80e-3));
    zCrop = z(idL:idH);
    xNew = x(1):param.gridSize:x(end);
    zNew = zCrop(1):param.gridSize:zCrop(end);
    [X,Z] = meshgrid(xNew,zNew);
    cfg = struct('ds',param.gridSize,'zMin',param.zMin,'weightTol',1e-12,'centre',true);
    solver = struct('innerSolve','pcg','pcgTol',param.pcgTol,'maxPcg',param.maxIter,'maxOuter',param.maxOuter,'epsVal',param.epsVal,'objTol',param.gradTol,'stationTol',param.gradTol,'kind','pseudoHuber');
    nAng = numel(angles);
    regStack = NaN([size(X),nAng]);
    rawStack = NaN(size(regStack));
    effStack = NaN(size(regStack));
    support = false(size(regStack));
    diag = cell(1,nAng);

    % Solve each physical column problem with the selected geometry and penalty
    for iAng = 1:nAng
        cum = interp2(x,zCrop,ba(idL:idH,:,iAng),X,Z);
        model = prepDeSc(cum,xNew,zNew,angles(iAng),cfg,true);
        if doTv
            [u,fit,uInit] = solveDeScPenalty(model.A,model.b,model.mask,param.mu,solver,[]);
            if ~fit.converged || fit.nInnerFail>0 || fit.init.flag~=0
                warning('bona:DeScSolve','Angle %g: outer convergence %d, initialization flag %d, inner failures %d.',angles(iAng),fit.converged,fit.init.flag,fit.nInnerFail)
            end
        else
            % Keep an accurate data-only reference and expose deliberate truncation
            dataIdx = find(sum(abs(model.A),1)' > 0);
            ata = model.A'*model.A;
            rhs = model.A'*model.b;
            uInit = zeros(size(model.A,2),1);
            [uInit(dataIdx),flag,relres,nIter] = pcg(ata(dataIdx,dataIdx),rhs(dataIdx),param.pcgTol,param.maxIter);
            u = uInit;
            init = struct('flag',flag,'relres',relres,'nIter',nIter,'provided',false);
            fit = struct('kind','none','nIter',0,'converged',flag==0,'capHit',flag~=0,'init',init,'nInnerFail',0);
        end

        % Register cell-centre estimates and reject partial interpolation support
        field = reshape(u,size(model.mask));
        raw = reshape(uInit,size(model.mask));
        field(~model.mask) = NaN;
        raw(~model.mask) = NaN;
        [regStack(:,:,iAng),support(:,:,iAng)] = registerDeSc(field,model,xNew,zNew,cfg,true,true);
        rawStack(:,:,iAng) = registerDeSc(raw,model,xNew,zNew,cfg,true,true);
        effStack(:,:,iAng) = cum;
        geom = rmfield(model,{'A','b','raw'});
        diag{iAng} = struct('angle',angles(iAng),'nObs',nnz(model.mask),'nCol',size(model.A,2),'solver',fit,'geom',geom);
    end

    % Average supported angle estimates on the same physical grid
    result.image = mean(regStack,3,'omitnan');
    result.noReg = mean(rawStack,3,'omitnan');
    result.eff = mean(effStack,3,'omitnan');
    result.angleMaps = regStack;
    result.angleSupport = support;
    result.nAngle = sum(isfinite(regStack),3);
    result.x = xNew;
    result.z = zNew;
    result.diag = diag;
    result.param = param;
    result.method = 'Physical DeSC with matched pseudo-Huber';
    result.converged = all(cellfun(@(d) d.solver.converged && d.solver.init.flag==0 && d.solver.nInnerFail==0,diag));

end
