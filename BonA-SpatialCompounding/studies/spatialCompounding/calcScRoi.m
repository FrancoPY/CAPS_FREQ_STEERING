function geom = calcScRoi(x,z,cfgId)

    % Define nominal phantom truth and historical metrics ROIs on centered axes
    [X,Z] = meshgrid(x,z);
    mid = (x(1)+x(end))/2;
    geom.truth = 6*ones(size(X));
    switch cfgId
        case 'saminc_ba6116'
            geom.circles = [mid,22.5e-3,9e-3];
            geom.rects = [mid-5e-3,17.5e-3,10e-3,10e-3;mid-17.5e-3,17.5e-3,5e-3,10e-3;mid+12.5e-3,17.5e-3,5e-3,10e-3;mid-5e-3,32.5e-3,10e-3,10e-3];
            geom.targets = [11,6,6,6];
            geom.names = {'inclusion','leftBg','rightBg','lowerBg'};
        case 'sam2inc_difsize_ba6116116'
            centers = mid+[-6.4,6.4]*1e-3;
            geom.circles = [-6.54e-3,22.5e-3,6.4e-3;6.54e-3,22.5e-3,5.1e-3];
            d1 = 5e-3/sqrt(2);
            d2 = 2.5e-3;
            geom.rects = [centers(1)-d1,19.5e-3,2*d1,6.5e-3;centers(2)-d2,19.5e-3,2*d2,6.5e-3;mid-18e-3,19.5e-3,4e-3,6.5e-3;mid+12.5e-3,19.5e-3,5.5e-3,6.5e-3;mid,19.5e-3,1.3e-3,6.5e-3;centers(1)-d1,32e-3,2*d1,6.5e-3;centers(2)-d2,32e-3,2*d2,6.5e-3];
            geom.targets = [11,11,6,6,6,6,6];
            geom.names = {'leftInc','rightInc','leftBg','rightBg','centerBg','lowerLeftBg','lowerRightBg'};
        case 'sam4inc_samesize_ba6116'
            geom.circles = [-9.81e-3,17.5e-3,4.5e-3;9.81e-3,17.5e-3,4.5e-3;-9.81e-3,33.5e-3,4.5e-3;9.81e-3,33.5e-3,4.5e-3];
            geom.rects = [geom.circles(:,1:2)-3e-3,6e-3*ones(4,2);mid-2e-3,15e-3,4e-3,25e-3];
            geom.targets = [11,11,11,11,6];
            geom.names = {'upperLeftInc','upperRightInc','lowerLeftInc','lowerRightInc','centerBg'};
        case 'vlayers_ba696'
            geom.circles = zeros(0,3);
            width = (x(end)-x(1))/3;
            geom.rects = [x(1)+0.15*width,16.5e-3,0.7*width,27e-3;x(1)+1.15*width,16.5e-3,0.7*width,27e-3;x(1)+2.15*width,16.5e-3,0.7*width,27e-3];
            geom.targets = [6,9,6];
            geom.names = {'leftBg','layer','rightBg'};
            geom.truth(abs(X) <= 6.54e-3) = 9;
    end
    for iInc = 1:size(geom.circles,1)
        circ = geom.circles(iInc,:);
        geom.truth((X-circ(1)).^2+(Z-circ(2)).^2 <= circ(3)^2) = 11;
    end
    geom.masks = cell(1,size(geom.rects,1));
    for iRoi = 1:size(geom.rects,1)
        rect = geom.rects(iRoi,:);
        geom.masks{iRoi} = X > rect(1) & X < rect(1)+rect(3) & Z > rect(2) & Z < rect(2)+rect(4);
    end
    % Preserve combined historical metrics alongside each displayed ROI
    switch cfgId
        case 'saminc_ba6116'
            geom.groups = {geom.masks{1},geom.masks{2}|geom.masks{3},geom.masks{4}};
            geom.groupNames = {'inc','latBg','lowBg'};
            geom.groupTargets = [11,6,6];
        case 'sam2inc_difsize_ba6116116'
            geom.groups = {geom.masks{1}|geom.masks{2},geom.masks{3}|geom.masks{4}|geom.masks{5},geom.masks{6}|geom.masks{7}};
            geom.groupNames = {'inc','latBg','lowBg'};
            geom.groupTargets = [11,6,6];
        case 'sam4inc_samesize_ba6116'
            geom.groups = {geom.masks{1}|geom.masks{2}|geom.masks{3}|geom.masks{4},geom.masks{5}};
            geom.groupNames = {'inc','bg'};
            geom.groupTargets = [11,6];
        case 'vlayers_ba696'
            geom.groups = geom.masks;
            geom.groupNames = {'left','layer','right'};
            geom.groupTargets = [6,9,6];
    end
    geom.viewBox = [x(1),x(end),10e-3,50e-3];

end
