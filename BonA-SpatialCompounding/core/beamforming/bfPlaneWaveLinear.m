function rfBf = bfPlaneWaveLinear(rfPrebf,geom)
    % Beamform one plane-wave channel with precomputed delay geometry.

    % Convert RF data to the delay-interpolation precision
    rf = double(rfPrebf);
    [numSamp,numLine] = size(rf);
    rfBf = zeros(numSamp,numLine);

    % Apply receive delays and the dynamic aperture to each image line
    for iLine = 1:numLine
        lowIdx = geom.lowIdx(:,:,iLine);
        valid = geom.validMask(:,:,iLine) & geom.recvMask(:,:,iLine);
        rfBf(:,iLine) = sumRfDelay(rf,lowIdx,geom.delayFrac(:,:,iLine),valid);
    end

end
