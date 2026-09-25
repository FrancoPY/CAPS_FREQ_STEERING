function [lowFull, highFull, xAxis, zAxis, fs] = loadRfStack(pathIn, studyName, frameRange, noiseLevel)
    
    % loadRfStack  Load low- and high-pressure RF stacks
    %
    % Inputs
    %   pathIn      : folder with beamformed RF files
    %   studyName   : study identifier
    %   frameRange  : frames to load
    %   noiseLevel  : additive noise level
    %
    % Outputs
    %   lowFull     : low-pressure RF stack
    %   highFull    : high-pressure RF stack
    %   xAxis       : lateral axis [m]
    %   zAxis       : axial axis [m]
    %   fs          : sampling frequency [Hz]

    % Init
    nFrames = numel(frameRange);
    xAxis = [];
    zAxis = [];
    fs = [];

    % Load frames
    for iFrame = 1:nFrames

        frame = frameRange(iFrame);

        S = load(fullfile(pathIn, sprintf('%s_f%d_HP.mat', studyName, frame)), ...
            'rfHp', 'xAxis', 'zAxis', 'fs');
        highNow = S.rfHp;

        if isfield(S, 'xAxis')
            xAxis = S.xAxis;
        end
        if isfield(S, 'zAxis')
            zAxis = S.zAxis;
        end
        if isfield(S, 'fs')
            fs = S.fs;
        end

        S = load(fullfile(pathIn, sprintf('%s_f%d_LP.mat', studyName, frame)), 'rfLp');
        lowNow = S.rfLp;

        if iFrame == 1
            lowFull = zeros(size(lowNow, 1), size(lowNow, 2), nFrames, 'like', lowNow);
            highFull = zeros(size(highNow, 1), size(highNow, 2), nFrames, 'like', highNow);
        end

        if noiseLevel == 0
            lowFull(:, :, iFrame) = lowNow;
            highFull(:, :, iFrame) = highNow;
        else
            lowFull(:, :, iFrame) = lowNow + noiseLevel * randn(size(lowNow), 'like', lowNow);
            highFull(:, :, iFrame) = highNow + noiseLevel * randn(size(highNow), 'like', highNow);
        end

    end

end