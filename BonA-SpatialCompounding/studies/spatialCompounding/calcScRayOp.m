function model = calcScRayOp(x,z,theta,srcZ,obs)
    % Calculate full-path ray averages on physical Cartesian cells (meters).

    % Define physical cell faces and retain any unmodeled shallow prefix
    dx = x(2)-x(1);
    dz = z(2)-z(1);
    xEdge = [x(1)-dx/2,(x(1:end-1)+x(2:end))/2,x(end)+dx/2];
    zEdge = [max(srcZ,z(1)-dz/2),(z(1:end-1)+z(2:end))/2,z(end)+dz/2];
    [X,Z] = meshgrid(x,z);
    [nRow,nCol] = size(X);
    nPix = numel(X);
    sx = sind(theta);
    cz = cosd(theta);
    ell = (Z-srcZ)/cz;
    xSrc = X-ell*sx;
    tol = 1e-12;
    valid = obs & ell>tol & xSrc>=xEdge(1)-tol & xSrc<=xEdge(end)+tol;
    ids = find(valid);
    rowCell = cell(numel(ids),1);
    colCell = cell(size(rowCell));
    valCell = cell(size(rowCell));
    prefix = zeros(size(X));

    % Split each ray at grid faces and identify cells from segment midpoints
    for k = 1:numel(ids)
        id = ids(k);
        tEnd = ell(id);
        tZ = (Z(id)-zEdge)/cz;
        if sx==0
            tX = [];
        else
            tX = (X(id)-xEdge)/sx;
        end
        t = sort([0,tEnd,tX(tX>tol & tX<tEnd-tol),tZ(tZ>tol & tZ<tEnd-tol)]);
        t = t([true,diff(t)>tol]);
        dt = diff(t);
        mid = (t(1:end-1)+t(2:end))/2;
        iCol = discretize(X(id)-sx*mid,xEdge);
        iRow = discretize(Z(id)-cz*mid,zEdge);
        keep = ~isnan(iCol) & ~isnan(iRow);
        prefix(id) = sum(dt(~keep))/tEnd;
        cols = iRow(keep)+nRow*(iCol(keep)-1);
        rowCell{k} = repmat(id,numel(cols),1);
        colCell{k} = cols(:);
        valCell{k} = dt(keep)'/tEnd;
    end
    A = sparse(vertcat(rowCell{:}),vertcat(colCell{:}),vertcat(valCell{:}),nPix,nPix);
    model = struct('A',A,'mask',valid,'ell',ell,'prefix',prefix,'xEdge',xEdge,'zEdge',zEdge,'srcZ',srcZ,'theta',theta,'size',[nRow,nCol]);

end
