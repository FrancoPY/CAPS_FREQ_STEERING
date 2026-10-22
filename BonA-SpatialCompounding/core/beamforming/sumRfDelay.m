function bf = sumRfDelay(rf,idx,frac,valid)

    % Interpolate time-by-element-by-shot RF and sum valid receive channels
    % Supply safe lower sample indices and fractions for one image line
    nSamp = size(rf,1);
    nElem = size(rf,2);
    nShot = size(rf,3);
    linIdx = double(idx) + (0:nElem - 1) * nSamp;
    bf = zeros(size(idx,1),nShot);

    % Bound the temporary interpolation arrays for multi-shot acquisitions
    for first = 1:8:nShot
        shots = first:min(first + 7,nShot);
        shotIdx = linIdx + reshape((shots - 1) * nSamp * nElem,1,1,[]);
        vals = (1 - frac) .* rf(shotIdx) + frac .* rf(shotIdx + 1);
        vals = vals .* valid;
        bf(:,shots) = reshape(sum(vals,2),size(idx,1),numel(shots));
    end

end
