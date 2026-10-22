clear; close all; clc

%% Define settings
studyDir = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(studyDir)),studyDir)
paths = bonaPaths();
addpath(paths.core{:})
cfgId = 'saminc_ba6116';
fc = 5e6;
alphaCoeff = 0.1;
alphaPower = 2;
pressures = [80,400]*1e3;
angles = -15:0.5:15;
descAngles = [0:2.5:15,-15:2.5:-2.5];
demuscAngles = -15:5:15;
frameSets = {1,2,3,4,1:4};
frameNames = {'f1','f2','f3','f4','pooled'};
latWidth = 4;
latStep = 2;
avgCycles = 10;
filterOrder = 200;
filterBand = [3,7]*1e6;
descParam = struct('gridSize',0.6e-3,'zMin',10e-3,'mu',0.1,'pcgTol',1e-8,'maxIter',1000,'epsVal',0.1,'gradTol',1e-6,'maxOuter',1000);
demuscParam = struct('gridSize',0.6e-3,'zMin',10e-3,'zMax',50e-3,'mu',10^-0.5,'tau',0.1,'maxIter',7);
dataDir = fullfile(paths.data,'studies','spatialCompounding',cfgId);
refDir = fullfile(paths.data,'studies','spatialCompounding','ref_ba6');
outDir = fullfile(paths.data,'results','studies','spatialCompounding',cfgId);
mkdir(outDir)

%% Prepare RF envelopes and cumulative B/A
sam = matfile(fullfile(dataDir,'rfBf.mat'));
ref = matfile(fullfile(refDir,'rfBf.mat'));
shotSam = sam.shotKey;
shotRef = ref.shotKey;
samDone = sam.completed;
refDone = ref.completed;
rfSize = size(sam,'rf');
nFrame = rfSize(3);
nLine = rfSize(2);
samFrames = unique([frameSets{:}]);
fs = sam.fs;
c0 = sam.c0;
z = sam.z;
x = calcLatAxis(sam.x,latWidth,latStep);
h = fir1(filterOrder,filterBand/(fs/2));
L = round(avgCycles*(c0/fc)/(z(2)-z(1)));
alpha1 = alphaCoeff*(fc/1e6)^alphaPower*100/(20*log10(exp(1)));
alpha2 = alphaCoeff*(2*fc/1e6)^alphaPower*100/(20*log10(exp(1)));
med = struct('alpha1',alpha1,'alpha2',alpha2,'c0',c0,'ba',6);
f1 = nan(numel(z),numel(x),numel(angles),numel(frameSets));
bmode = zeros(size(f1));
for iAng = 1:numel(angles)
    loSam = find(shotSam(1,:) == pressures(1) & shotSam(2,:) == angles(iAng));
    hiSam = find(shotSam(1,:) == pressures(2) & shotSam(2,:) == angles(iAng));
    loRef = find(shotRef(1,:) == pressures(1) & shotRef(2,:) == angles(iAng));
    hiRef = find(shotRef(1,:) == pressures(2) & shotRef(2,:) == angles(iAng));
    assert(all(samDone(samFrames,[loSam,hiSam]),'all') && all(refDone(1,[loRef,hiRef]),'all'),'Run the study beamforming scripts for these angles.');
    envLo = zeros(numel(z),nLine,nFrame);
    envHi = zeros(size(envLo));
    for iFrame = samFrames
        envLo(:,:,iFrame) = calcRfEnvelope(sam.rf(:,:,iFrame,loSam),h,1);
        envHi(:,:,iFrame) = calcRfEnvelope(sam.rf(:,:,iFrame,hiSam),h,1);
    end
    lowRef = mean(calcRfEnvelope(ref.rf(:,:,1,loRef),h,1),2);
    highRef = mean(calcRfEnvelope(ref.rf(:,:,1,hiRef),h,1),2);
    for iSet = 1:numel(frameSets)
        frames = frameSets{iSet};
        for iCol = 1:numel(x)
            cols = (iCol-1)*latStep+(1:latWidth);
            lowSam = mean(envLo(:,cols,frames),[2,3]);
            highSam = mean(envHi(:,cols,frames),[2,3]);
            pSam = sqrt(abs(pressures(2)/pressures(1)*lowSam-highSam).*lowRef);
            pRef = sqrt(abs(pressures(2)/pressures(1)*lowRef-highRef).*lowSam);
            f1(:,iCol,iAng,iSet) = movmean(pSam,L)./movmean(pRef,L);
            bmode(:,iCol,iAng,iSet) = lowSam;
        end
    end
end
prepParam = struct('fc',fc,'alphaCoeff',alphaCoeff,'alphaPower',alphaPower,'filterBand',filterBand,'filterOrder',filterOrder,'avgCycles',avgCycles,'latWidth',latWidth,'latStep',latStep,'pressures',pressures);
save(fullfile(dataDir,'prepared.mat'),'f1','bmode','x','z','angles','frameSets','frameNames','prepParam','med','-v7.3')

%% Estimate separate frames and pooled envelopes
load(fullfile(dataDir,'prepared.mat'),'f1','bmode','x','z','angles','frameSets','frameNames','prepParam','med')
[~,descIdx] = ismember(descAngles,angles);
[~,demuscIdx] = ismember(demuscAngles,angles);
[~,broadIdx] = ismember(0,angles);
models = {'quadratic','general'};
modelNames = {'DM','EDM'};
for iSet = 1:numel(frameSets)
    baDm = calcDeplBa(f1(:,:,:,iSet),z,med,med,'quadratic');
    baEdm = calcDeplBa(f1(:,:,:,iSet),z,med,med,'general');
    modelDiff = max(abs(baDm-baEdm),[],'all','omitnan');
    for iModel = 1:numel(models)
        ba = calcDeplBa(f1(:,:,:,iSet),z,med,med,models{iModel});
        if iModel == 1 || ~isequaln(baDm,baEdm)
            desc = invertDeSc(ba(:,:,descIdx),x,z,descAngles,descParam);
            prob = prepDeMuSc(ba(:,:,demuscIdx),x,z,demuscAngles,demuscParam);
            demusc = invertDeMuSc(prob,demuscParam);
            broad = invertDeSc(ba(:,:,broadIdx),x,z,0,descParam);
        end
        xLoc = desc.x;
        zLoc = desc.z;
        geom = calcScRoi(xLoc,zLoc,cfgId);
        maps = {broad.image,desc.image,demusc.image};
        stats = cell(size(maps));
        for iMap = 1:numel(maps)
            for iRoi = 1:numel(geom.masks)
                stats{iMap}.(geom.names{iRoi}) = calcRoiStats(maps{iMap},geom.masks{iRoi},geom.targets(iRoi),false);
            end
            for iGroup = 1:numel(geom.groups)
                stats{iMap}.groups.(geom.groupNames{iGroup}) = calcRoiStats(maps{iMap},geom.groups{iGroup},geom.groupTargets(iGroup),false);
            end
            valid = isfinite(maps{iMap}) & zLoc(:) >= 10e-3 & zLoc(:) <= 50e-3;
            stats{iMap}.mae = mean(abs(maps{iMap}(valid)-geom.truth(valid)));
            stats{iMap}.maeCount = nnz(valid);
        end
        baBroad = ba(:,:,broadIdx);
        env = interp2(x,z,bmode(:,:,broadIdx,iSet),xLoc,zLoc(:));
        tag = [modelNames{iModel},'_',frameNames{iSet}];
        runParam = struct('cfgId',cfgId,'model',models{iModel},'frames',frameSets{iSet},'refFrames',1,'descAngles',descAngles,'demuscAngles',demuscAngles);
        save(fullfile(outDir,[tag,'.mat']),'desc','demusc','broad','stats','geom','baBroad','x','z','env','xLoc','zLoc','prepParam','descParam','demuscParam','modelDiff','runParam','-v7.3')
        plotScResult([{geom.truth},maps],{'Truth','Broadside','DeSC','DeMuSC'},xLoc,zLoc,env,geom,fullfile(outDir,[tag,'_maps.png']),false);
        plotScResult(maps,{'Broadside','DeSC','DeMuSC'},xLoc,zLoc,env,geom,fullfile(outDir,[tag,'_roi.png']),true);
        fprintf('%s %s: DM/EDM max difference %.3g\n',cfgId,tag,modelDiff);
    end
end
