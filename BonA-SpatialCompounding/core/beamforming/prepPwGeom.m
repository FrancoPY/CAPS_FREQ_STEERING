function geom = prepPwGeom(numSamp,numLine,fs,cBf,pitch,fNum,steerDeg)
    % Prepare delay indices and the receive aperture for plane-wave beamforming.

    % Verify the beamform parameters
    if numSamp < 2 || numSamp ~= round(numSamp)
        error('The RF sample count must be an integer greater than one.')
    end

    if numLine < 1 || numLine ~= round(numLine)
        error('The image-line count must be a positive integer.')
    end

    if fs <= 0 || cBf <= 0 || pitch <= 0 || fNum <= 0
        error('The sample rate, sound speed, pitch, and F-number must be positive.')
    end

    % Define the lateral and depth axes
    timeAxis = (0:numSamp - 1)' / double(fs);
    zAxis = timeAxis * double(cBf) / 2;
    xAxis = ((0:numLine - 1) - (numLine - 1) / 2) * double(pitch);
    steerRad = steerDeg * pi / 180;

    % Initialize the reusable delay arrays
    lowIdx = zeros(numSamp,numLine,numLine,'uint32');
    delayFrac = zeros(numSamp,numLine,numLine);
    validMask = false(numSamp,numLine,numLine);
    recvMask = false(numSamp,numLine,numLine);
    elemOff = uint32(reshape((0:numLine - 1) * numSamp,1,numLine));

    % Calculate the transmit and receive paths for each image line
    for iLine = 1:numLine
        xPos = xAxis(iLine);
        txDist = zAxis * cos(steerRad) + xPos * sin(steerRad);
        rxDist = sqrt(zAxis.^2 + (xPos - xAxis).^2);
        tof = (txDist + rxDist) / double(cBf);
        sampPos = tof * double(fs) + 1;
        lowNow = floor(sampPos);
        validNow = lowNow >= 1 & lowNow < numSamp;
        fracNow = sampPos - lowNow;
        lowNow(~validNow) = 1;
        fracNow(~validNow) = 0;
        lowIdx(:,:,iLine) = uint32(lowNow);
        delayFrac(:,:,iLine) = fracNow;
        validMask(:,:,iLine) = validNow;
        recvMask(:,:,iLine) = abs(xAxis - xPos) < zAxis / (2 * fNum);
    end

    % Organize the beamform geometry
    geom.xAxis = xAxis;
    geom.zAxis = zAxis;
    geom.lowIdx = lowIdx;
    geom.delayFrac = delayFrac;
    geom.validMask = validMask;
    geom.recvMask = recvMask;
    geom.elemOff = elemOff;

end
