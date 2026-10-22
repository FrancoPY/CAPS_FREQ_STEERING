function result = invertDeMuSc(prob,param)

    % Retain the original CGS unregularized estimate with uniform B/A six
    A = prob.A;
    y = prob.y;
    [raw,flag,relRes,nIter] = cgs(A'*A,A'*y,1e-5,param.maxIter,[],[],6*ones(size(A,2),1));
    result.noReg = reshape(raw,prob.size);
    result.cgs = struct('flag',flag,'relRes',relRes,'nIter',nIter);

    % Retain the separate PCG initialization and primal-dual TV solver
    [result.image,result.diag] = solveDeMuScTv(y,A,param.mu,param.tau,param.maxIter,prob.size(1),prob.size(2),'pcg');
    result.x = prob.x;
    result.z = prob.z;
    result.srcZ = prob.srcZ;
    result.obsMask = prob.mask;
    result.param = param;
    result.fitRmse = norm(A*result.image(:)-y)/sqrt(numel(y));

end
