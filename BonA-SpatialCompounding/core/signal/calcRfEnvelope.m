function env = calcRfEnvelope(rf, hFilter, L, useGpu, gpuClass)

    % Filter time-by-channel RF and smooth its analytic envelope with zero-padded edges
    if nargin < 4
        useGpu = false;
    end
    if nargin < 5
        gpuClass = 'single';
    end

    % Set the RF and filter precision on the GPU
    x = rf;
    h = hFilter(:);
    if useGpu
        x = gpuArray(cast(x, gpuClass));
        h = gpuArray(cast(h, gpuClass));
    end

    % Construct the FFT analytic signal along the time dimension
    xf = conv2(x, h, 'same');
    n = size(xf, 1);
    h = zeros(n, 1, 'like', xf);
    h(1) = 1;
    if mod(n, 2) == 0
        h(n / 2 + 1) = 1;
        h(2:n / 2) = 2;
    else
        h(2:(n + 1) / 2) = 2;
    end
    xa = ifft(fft(xf, [], 1) .* h, [], 1);
    avgKer = ones(L, 1, 'like', x) / L;
    env = conv2(abs(xa), avgKer, 'same');

    % Return host data in the RF input class
    if useGpu
        env = cast(gather(env), 'like', rf);
    end

end
