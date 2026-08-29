function bfData = bfPlaneWaveSimu(sensorDataKWave, fs, fNumber, steeringAngle, soundSpeed, offsetBaseSamples)

    % bfPlaneWaveSimu  Plane-wave DAS beamforming for simulated channel data
    %
    % Inputs
    %   sensorDataKWave    : channel data from k-Wave (SIN recortar; incluye
    %                        el retardo común offset_base al inicio, tal
    %                        como lo graba k-Wave)
    %   fs                 : sampling frequency [Hz]
    %   fNumber            : receive f-number
    %   steeringAngle      : transmit steering angle [deg]
    %   soundSpeed         : sound speed used for delays [m/s] (default
    %                        1500, debe coincidir con c0 de la simulación)
    %   offsetBaseSamples  : retardo común de transmisión (offset_base,
    %                        en muestras) usado para el steering en el
    %                        script de simulación. Todos los canales
    %                        arrancan a disparar recién en esta muestra,
    %                        así que el modelo de retardo (que asume t=0
    %                        en el instante real de disparo) se traduce a
    %                        índice de dato sumando este valor, NUNCA
    %                        recortando el array (eso borra ecos reales
    %                        del lado donde projDist es negativo).
    %                        Default 0 (para datos sin steering).
    %
    % Output
    %   bfData          : beamformed RF image [nSamplesOut x nLines]

    if nargin < 5 || isempty(soundSpeed)
        soundSpeed = 1500;
    end
    if nargin < 6 || isempty(offsetBaseSamples)
        offsetBaseSamples = 0;
    end

    % Data and geometry
    data = sensorDataKWave.';
    nSamplesFull = size(data, 1);
    nSamplesOut = nSamplesFull - offsetBaseSamples;
    nElem = size(data, 2);

    elementSpacing = 0.3e-3;
    elemPos = ((1:nElem) * elementSpacing) - mean((1:nElem) * elementSpacing);

    steerRad = steeringAngle * pi / 180;

    % Output
    depthAxis = (1:nSamplesOut) * (soundSpeed / fs) / 2;
    nLines = nElem;
    bfData = zeros(nSamplesOut, nLines);

    % DAS loop
    for ix = 1:nLines
        xPoint = elemPos(ix);

        for iz = 1:nSamplesOut
            zPoint = depthAxis(iz);

            % Delays
            projDist = zPoint * cos(steerRad) + xPoint * sin(steerRad);
            rxDist = sqrt(zPoint^2 + (xPoint - elemPos).^2);
            totalTau = (projDist + rxDist) / soundSpeed;

            idx = round(totalTau * fs) + offsetBaseSamples + 1;
            validIdx = idx >= 1 & idx <= nSamplesFull;
            idx(~validIdx) = 1;   % clamp only to allow indexing; excluded below via win

            % Samples
            linIdx = idx + (0:nElem-1) * nSamplesFull;
            sampVals = data(linIdx);

            apertureRadius = zPoint / (2 * fNumber);
            win = (abs(elemPos - xPoint) < apertureRadius) & validIdx;

            % Aperture and sum
            bfData(iz, ix) = sum(win .* sampVals);

        end
    end

end