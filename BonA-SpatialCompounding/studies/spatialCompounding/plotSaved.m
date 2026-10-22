clear; close all; clc

%% Define the saved experiment
studyDir = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(studyDir)),studyDir)
paths = bonaPaths();
addpath(paths.core{:})
cfgId = 'sam2inc_difsize_ba6116116';
frameName = 'f2';
models = {'DM','EDM'};
outDir = fullfile(paths.data,'results','studies','spatialCompounding',cfgId);

%% Refresh figures without repeating estimation
for iModel = 1:numel(models)
    tag = [models{iModel},'_',frameName];
    data = load(fullfile(outDir,[tag,'.mat']));
    maps = {data.broad.image,data.desc.image,data.demusc.image};
    plotScResult([{data.geom.truth},maps],{'Truth','Broadside','DeSC','DeMuSC'},data.xLoc,data.zLoc,data.env,data.geom,fullfile(outDir,[tag,'_maps.png']),false);
    plotScResult(maps,{'Broadside','DeSC','DeMuSC'},data.xLoc,data.zLoc,data.env,data.geom,fullfile(outDir,[tag,'_roi.png']),true);
end
