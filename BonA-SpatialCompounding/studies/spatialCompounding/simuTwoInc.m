clear; close all; clc

%% Define the acquisition
studyDir = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(studyDir)))
paths = bonaPaths();
addpath(paths.core{:})
kWaveDir = 'C:/Users/eavil/AppData/Roaming/MathWorks/MATLAB Add-Ons/Collections/k-Wave/k-Wave';
addpath(kWaveDir)
cfgId = 'sam2inc_difsize_ba6116116';
angles = -15:0.5:15;
pressures = [80,400]*1e3;
nFrame = 4;
c0 = 1540;
rho0 = 1000;
fc = 5e6;
nCycle = 10;
alphaCoeff = 0.1;
alphaPower = 2;
gridSize = [944,655,140];
dx = 0.06e-3;
pml = [40,10,10];
densCv = 0.02;
pitch = 0.3e-3;
nElem = 128;
elemWidth = 5;
elemLen = 100;
elemSpace = 0;
txApo = 'Rectangular';
rxApo = 'Rectangular';
baBg = 6;
baInc = 11;
incCenters = [22.5e-3,1/3;22.5e-3,2/3];
incRadii = [6.4,5.1]*1e-3;
layerFrac = [];
dataDir = fullfile(paths.data,'studies','spatialCompounding',cfgId);
densDir = fullfile(paths.data,'studies','spatialCompounding','density');
rfFile = fullfile(dataDir,'rfSim.mat');
mkdir(dataDir)

%% Reuse densities and resume matched shots
simParam = struct('cfgId',cfgId,'c0',c0,'rho0',rho0,'fc',fc,'nCycle',nCycle,'alphaCoeff',alphaCoeff,'alphaPower',alphaPower,'gridSize',gridSize,'dx',dx,'pitch',pitch,'pml',pml);
simParam.phantom = struct('baBg',baBg,'baInc',baInc,'incCenters',incCenters,'incRadii',incRadii,'layerFrac',layerFrac);
simParam.probe = struct('nElem',nElem,'elemWidth',elemWidth,'elemLen',elemLen,'elemSpace',elemSpace,'txApo',txApo,'rxApo',rxApo);
assert(abs(pitch-(elemWidth+elemSpace)*dx) <= 4*eps(pitch),'Match the nominal pitch to the probe grid geometry.');
densNames = cell(1,nFrame);
densHash = cell(1,nFrame);
for iFrame = 1:nFrame
    densNames{iFrame} = sprintf('BAPW_STD2_SAM2023_%d.mat',iFrame);
    [~,densHash{iFrame}] = loadSimDensity(fullfile(densDir,densNames{iFrame}),gridSize,rho0,densCv,'density');
end
if ~isfile(rfFile)
    save(rfFile,'simParam','-v7.3','-nocompression')
end
[P,Theta] = ndgrid(pressures,angles);
txKey = [P(:)';Theta(:)'];
[out,state,txIdx] = syncSimRf(rfFile,simParam,densDir,densHash,txKey);
out.shotKey = state.txKey;
out.completed = state.done;
out.simParam = simParam;

%% Define the medium and transducer
kgrid = kWaveGrid(gridSize(1),dx,gridSize(2),dx,gridSize(3),dx);
kgrid.makeTime(c0,[],2.2*gridSize(1)*dx/c0);
fs = 1/kgrid.dt;
out.fs = fs;
out.c0 = c0;
if ~ismember('rf',who(out))
    initRfArray(rfFile,[numel(kgrid.t_array),nElem,nFrame,size(state.txKey,2)]);
end
medium.sound_speed = c0;
medium.sound_speed_ref = c0;
medium.alpha_coeff = alphaCoeff;
medium.alpha_power = alphaPower;
ba = baBg*ones(gridSize(1:2));
for iInc = 1:size(incCenters,1)
    inc = makeDisc(gridSize(1),gridSize(2),round(incCenters(iInc,1)/dx),round(incCenters(iInc,2)*gridSize(2)),round(incRadii(iInc)/dx));
    ba = ba+(baInc-baBg)*inc;
end
medium.BonA = repmat(single(ba),1,1,gridSize(3));
trans.number_elements = nElem;
trans.element_width = elemWidth;
trans.element_length = elemLen;
trans.element_spacing = elemSpace;
trans.radius = inf;
probeWidth = nElem*elemWidth+(nElem-1)*elemSpace;
trans.position = round([1,gridSize(2)/2-probeWidth/2,gridSize(3)/2-elemLen/2]);
trans.sound_speed = c0;
trans.focus_distance = inf;
trans.elevation_focus_distance = inf;
trans.transmit_apodization = txApo;
trans.receive_apodization = rxApo;
trans.active_elements = ones(nElem,1);
pulse = toneBurst(fs,fc,nCycle);

%% Simulate only missing frame and shot identities
for iFrame = 1:nFrame
    dens = load(fullfile(densDir,densNames{iFrame}),'density');
    medium.density = dens.density;
    for iTx = 1:size(txKey,2)
        iShot = txIdx(iTx);
        if state.done(iFrame,iShot)
            continue
        end
        trans.steering_angle = txKey(2,iTx);
        trans.input_signal = txKey(1,iTx)/(c0*rho0)*pulse;
        probe = kWaveTransducer(kgrid,trans);
        sensor = kspaceFirstOrder3DG(kgrid,medium,probe,probe,'PMLInside',false,'PMLSize',pml,'DataCast','gpuArray-single','DataRecast',true,'PlotSim',false);
        out.rf(:,:,iFrame,iShot) = double(sensor');
        state.done(iFrame,iShot) = true;
        out.simState = state;
        out.completed = state.done;
        fprintf('%s: frame %d, %.1f degrees, %.0f kPa\n',cfgId,iFrame,txKey(2,iTx),txKey(1,iTx)/1e3);
    end
end
