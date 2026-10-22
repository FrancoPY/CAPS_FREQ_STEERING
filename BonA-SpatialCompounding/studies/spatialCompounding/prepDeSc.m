function model = prepDeSc(cum,x,z,theta,cfg,fixPad)
    % Prepare uniform-grid cumulative B/A; x, z and cfg.ds are in meters.

    % Preserve source observation selection while fixing only the padding switch
    [nInput,nCol] = size(cum);
    nPad = abs(round(nInput*tand(theta)));
    pad = NaN(nInput,nCol+nPad);
    good = false(size(pad));
    if theta >= 0
        pad(:,1:nCol) = cum;
        good(:,1:nCol) = true;
        if theta > 0 && ~fixPad
            pad(:,nCol) = NaN;
            good(:,nCol) = false;
        end
        nLost = nnz(~good(:,1:nCol));
    else
        pad(:,nPad+(1:nCol)) = cum;
        good(:,nPad+(1:nCol)) = true;
        nLost = nnz(~good(:,nPad+(1:nCol)));
    end
    image = imrotate(pad,-theta,'bilinear','loose');
    weight = imrotate(double(good),-theta,'bilinear','loose');
    outside = ~imrotate(true(size(pad)),-theta,'bilinear','loose');
    image(outside & ~imclearborder(outside)) = NaN;
    [nRow,nRotCol] = size(image);
    [dc,dr] = meshgrid((1:nRotCol)-(nRotCol+1)/2,(1:nRow)-(nRow+1)/2);
    padMid = mean(x([1,end]))+sign(theta)*nPad*cfg.ds/2;
    zMid = mean(z([1,end]));
    ell = zMid/cosd(theta)+cfg.ds*(dr-tand(theta)*dc);
    u = padMid*cosd(theta)-zMid*sind(theta)+cfg.ds*dc;
    image(ell <= 0) = NaN;
    if theta < 0
        idLx = find(any(isfinite(image),1),1,'first');
        idHx = find(isfinite(image(2,:)),1,'last');
    elseif theta == 0
        idLx = 1;
        idHx = nRotCol;
    else
        idLx = find(isfinite(image(2,:)),1,'first');
        idHx = find(any(isfinite(image),1),1,'last');
    end
    cols = idLx:idHx;
    image = image(:,cols);
    ell = ell(:,cols);
    u = u(:,cols);
    full = weight(:,cols) > 1-cfg.weightTol;
    nCol = numel(cols);
    nEnd = zeros(1,nCol);
    for col = 1:nCol
        last = find(isfinite(image(:,col)),1,'last');
        if isempty(last)
            nEnd(col) = nRow;
        else
            nEnd(col) = nRow-last;
        end
    end
    [minEnd,minIdx] = min(nEnd);
    idH = round(cfg.zMin/cfg.ds);
    target = minEnd+idH;
    favorRight = theta < 0;
    blocks = cell(1,nCol);
    b = zeros(size(image));
    mask = false(size(image));
    raw = NaN(size(image));
    path = repmat(struct('anchor',NaN,'last',NaN,'nRow',0,'trim',0),1,nCol);
    startNan = nRow;
    endNan = nRow;
    for col = 1:nCol
        blocks{col} = sparse(nRow,nRow);
        first = find(isfinite(image(:,col)),1,'first');
        if isempty(first)
            startNan = nRow;
            endNan = nRow;
            continue
        end
        startNan = first-1;
        endNan = nEnd(col);
        onFavor = (favorRight && col > minIdx) || (~favorRight && col < minIdx);
        onOther = (favorRight && col < minIdx) || (~favorRight && col > minIdx);
        if onFavor
            trim = idH;
        elseif onOther && endNan > target
            trim = 0;
        else
            trim = target-endNan;
        end
        last = nRow-endNan-trim;
        anchor = find(isfinite(image(:,col)) & ell(:,col)>0 & ell(:,col)<=cfg.zMin+32*eps(cfg.zMin),1,'last');
        if isempty(anchor) || last <= anchor
            continue
        end
        rows = (anchor+1:last)';
        ellA = ell(anchor,col);
        ellObs = ell(rows,col);
        step = diff([ellA;ellObs]);
        block = (tril(ones(numel(rows),numel(rows))).*step')./ellObs;
        blocks{col}(rows,rows) = sparse(block);
        b(rows,col) = image(rows,col)-(ellA./ellObs)*image(anchor,col);
        raw(rows,col) = diff([ellA;ellObs].*image([anchor;rows],col))./step;
        mask(rows,col) = true;
        path(col) = struct('anchor',anchor,'last',last,'nRow',numel(rows),'trim',trim);
        assert(all(isfinite(image([anchor;rows],col))),'The candidate observation window contains missing cumulative data.');
    end
    A = blkdiag(blocks{:});
    xNode = padMid+cfg.ds*(dc(:,cols)*cosd(theta)+dr(:,cols)*sind(theta));
    zNode = zMid+cfg.ds*(-dc(:,cols)*sind(theta)+dr(:,cols)*cosd(theta));
    model = struct('A',A,'b',b(:),'mask',mask,'raw',raw,'xCell',xNode-0.5*cfg.ds*sind(theta),'zCell',zNode-0.5*cfg.ds*cosd(theta),'ell',ell,'u',u,'full',full,'path',path,'theta',theta,'nPad',nPad,'nLost',nLost,'nFalseObs',nnz(mask & ~full),'padMid',padMid,'zMid',zMid,'rotSize',[nRow,nRotCol],'idLx',idLx,'idHx',idHx,'startNan',startNan,'endNan',endNan);

end

