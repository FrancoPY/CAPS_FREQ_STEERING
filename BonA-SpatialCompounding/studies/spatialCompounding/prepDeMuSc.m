function prob = prepDeMuSc(ba,x,z,angles,param)

    % Interpolate effective B/A and assemble the full-depth multi-view system
    prob.x = x(1):param.gridSize:x(end);
    prob.z = z(1):param.gridSize:z(end);
    [Xq,Zq] = meshgrid(prob.x,prob.z);
    prob.srcZ = 0;
    [nRow,nCol] = size(Xq);
    aCell = cell(1,numel(angles));
    yCell = cell(1,numel(angles));
    maskCell = cell(1,numel(angles));
    for iAng = 1:numel(angles)
        baGrid = interp2(x,z,ba(:,:,iAng),Xq,Zq);
        ell = (Zq-prob.srcZ)/cosd(angles(iAng));
        obs = ell>=param.zMin-1e-12 & Zq<=param.zMax+1e-12;
        model = calcScRayOp(prob.x,prob.z,angles(iAng),prob.srcZ,obs);
        assert(all(model.prefix(model.mask)<1e-12),'The full-depth joint system cannot omit a shallow contribution.');
        aCell{iAng} = model.A(model.mask(:),:);
        yCell{iAng} = baGrid(model.mask);
        maskCell{iAng} = model.mask;
    end
    prob.A = vertcat(aCell{:});
    prob.y = vertcat(yCell{:});
    assert(all(isfinite(prob.y)),'The selected inversion observations contain undefined measured B/A.');
    prob.size = [nRow,nCol];
    prob.mask = maskCell;
    prob.param = param;

end
