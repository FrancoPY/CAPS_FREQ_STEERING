function [map,fullMap] = registerDeSc(field,model,x,z,cfg,fixSupport,fixCoords)

    % Keep half-cell resampling fixed before either registration calculation
    [nRow,nCol] = size(field);
    if cfg.centre
        field = interp2(1:nCol,(1:nRow)',field,1:nCol,(1:nRow)'+0.5,'linear',NaN);
    end
    valid = isfinite(field);
    if fixCoords
        % Query the broadside physical grid directly in the rotated canvas
        [X,Z] = meshgrid(x,z);
        colQ = ((X-model.padMid)*cosd(model.theta)-(Z-model.zMid)*sind(model.theta))/cfg.ds+(model.rotSize(2)+1)/2-model.idLx+1;
        rowQ = ((X-model.padMid)*sind(model.theta)+(Z-model.zMid)*cosd(model.theta))/cfg.ds+(model.rotSize(1)+1)/2;
        padded = zeros(nRow+2,nCol+2);
        padded(2:end-1,2:end-1) = field;
        good = zeros(size(padded));
        good(2:end-1,2:end-1) = valid;
        map = interp2(0:nCol+1,0:nRow+1,padded,colQ,rowQ,'linear',0);
        weight = interp2(0:nCol+1,0:nRow+1,good,colQ,rowQ,'linear',0);
    else
        % Preserve the source crop so single-switch comparisons isolate changes
        back = imrotate(field,model.theta,'bilinear','loose');
        weight = imrotate(double(valid),model.theta,'bilinear','loose');
        outside = ~imrotate(true(nRow,nCol),model.theta,'bilinear','loose');
        back(outside & ~imclearborder(outside)) = NaN;
        if model.theta < 0
            range = [max(1,model.endNan),size(back,1)-model.endNan,max(1,model.nPad),size(back,2)];
        elseif model.theta == 0
            range = [1,size(back,1),1,size(back,2)];
        else
            range = [max(1,model.startNan),size(back,1)-model.startNan,1,max(size(back,2)-model.nPad,1)];
        end
        map = fitSource(back(range(1):range(2),range(3):range(4)),[numel(z),numel(x)]);
        map = fitSource(map,[numel(z),numel(x)]);
        weight = fitSource(weight(range(1):range(2),range(3):range(4)),[numel(z),numel(x)]);
        weight = fitSource(weight,[numel(z),numel(x)]);
    end
    fullMap = weight > 1-cfg.weightTol;
    map(map==0) = NaN;
    if fixSupport
        map(~fullMap) = NaN;
    end

end

function out = fitSource(img,targetSize)

    % Reproduce the two source pad/crop passes for baseline comparisons
    [nRow,nCol] = size(img);
    out = zeros(targetSize);
    if nRow < targetSize(1) || nCol < targetSize(2)
        top = max(0,floor((targetSize(1)-nRow)/2));
        left = max(0,floor((targetSize(2)-nCol)/2));
        out((1:nRow)+top,(1:nCol)+left) = img;
    elseif nRow > targetSize(1) || nCol > targetSize(2)
        top = max(0,floor((nRow-targetSize(1))/2));
        left = max(0,ceil((nCol-targetSize(2))/2));
        out = img(top+1:min(nRow,top+targetSize(1)),left+1:min(nCol,left+targetSize(2)));
    else
        out = img;
    end

end

