clear; close all; clc

%% Define settings
studyDir = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(studyDir)))
paths = bonaPaths();
addpath(paths.core{:})
cfgId = 'saminc_ba6116';
angles = -15:0.5:15;
fNum = 3;
dataDir = fullfile(paths.data,'studies','spatialCompounding',cfgId);
rawFile = fullfile(dataDir,'rfSim.mat');
bfFile = fullfile(dataDir,'rfBf.mat');

%% Prepare the beamformed checkpoint
raw = matfile(rawFile);
fs = raw.fs;
c0 = raw.c0;
simParam = raw.simParam;
pitch = simParam.pitch;
shotKey = raw.shotKey;
rawDone = raw.completed;
[nSamp,nElem,nFrame,nShot] = size(raw,'rf');
bfParam = struct('method','linearZeroTimeV1','c0',c0,'fs',fs,'pitch',pitch,'fNum',fNum);
x = ((0:nElem-1)-(nElem-1)/2)*pitch;
z = (0:nSamp-1)'*c0/fs/2;
newFile = ~isfile(bfFile);
if ~newFile
    old = matfile(bfFile);
    oldSize = size(old,'rf');
    newFile = ~isequal(oldSize(1:2),[nSamp,nElem]);
    clear old
end
if newFile
    completed = false(nFrame,nShot);
    rfHash = cell(nFrame,nShot);
    save(bfFile,'completed','rfHash','bfParam','shotKey','fs','c0','x','z','-v7.3','-nocompression')
    initRfArray(bfFile,[nSamp,nElem,nFrame,nShot]);
end
out = matfile(bfFile,'Writable',true);
out.x = x;
out.z = z;
out.fs = fs;
out.c0 = c0;
completed = out.completed;
rfHash = out.rfHash;
oldKey = out.shotKey;
if size(completed,1) < nFrame || size(completed,2) < nShot
    completed(nFrame,nShot) = false;
    rfHash{nFrame,nShot} = '';
end
if ~isequal(out.bfParam,bfParam) || ~isequal(oldKey,shotKey(:,1:min(size(oldKey,2),nShot)))
    completed = false(nFrame,nShot);
end
out.completed = completed;
out.bfParam = bfParam;
out.shotKey = shotKey;
simState = raw.simState;
for iFrame = 1:nFrame
    densName = sprintf('BAPW_STD2_SAM2023_%d.mat',iFrame);
    assert(strcmp(hashSimFile(fullfile(simState.densDir,densName)),simState.densHash{iFrame}),'Run the simulation for the replaced density map.');
end

%% Reconstruct selected angles and retain matching frames
for iAng = 1:numel(angles)
    shots = find(shotKey(2,:) == angles(iAng));
    makeGeom = true;
    for iShot = shots
        for iFrame = 1:nFrame
            assert(rawDone(iFrame,iShot),'Run the matching simulation before beamforming this shot.');
            rf = raw.rf(:,:,iFrame,iShot);
            hash = hashSimArray(rf);
            if completed(iFrame,iShot) && strcmp(rfHash{iFrame,iShot},hash)
                continue
            end
            if makeGeom
                geom = prepPwGeom(nSamp,nElem,fs,c0,pitch,fNum,angles(iAng));
                makeGeom = false;
            end
            completed(iFrame,iShot) = false;
            out.completed = completed;
            out.rf(:,:,iFrame,iShot) = bfPlaneWaveLinear(rf,geom);
            rfHash{iFrame,iShot} = hash;
            completed(iFrame,iShot) = true;
            out.rfHash = rfHash;
            out.completed = completed;
        end
    end
    fprintf('%s: beamformed %.1f degrees\n',cfgId,angles(iAng));
end
out.bfParam = bfParam;
out.shotKey = shotKey;
fprintf('Beamformed checkpoint: %s\n',bfFile);
