function paths = bonaPaths()

    % Locate code, external datasets and disposable outputs
    paths.project = fileparts(mfilename('fullpath'));
    paths.data = 'D:\Research Datasets\BonA-Imaging';
    paths.tmp = fullfile(paths.project, 'tmp');

    % Locate the shared numerical functions
    names = {'beamforming', 'signal', 'geometry', 'metrics', 'model', 'simulation', 'plotting'};
    paths.core = cellfun(@(name) fullfile(paths.project, 'core', name), names, 'UniformOutput', false);

end
