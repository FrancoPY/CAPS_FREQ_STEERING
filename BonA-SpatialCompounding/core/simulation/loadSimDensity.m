function [density,hash] = loadSimDensity(file,gridSize,rho0,densCv,varName)
    % Reuse saved scatterers or save a new realization and its RNG state.

    % Generate missing scatterers once for every matched acquisition
    if ~isfile(file)
        if ~isfolder(fileparts(file))
            mkdir(fileparts(file))
        end
        rng('shuffle')
        data.rngState = rng;
        data.(varName) = single(rho0 + rho0 * densCv * randn(gridSize));
        save(file,'-struct','data','-v7.3')
    end

    % Verify the grid before using the saved density realization
    data = load(file,varName);
    density = single(data.(varName));
    assert(isequal(size(density),gridSize),'The density map grid differs from this simulation.');
    hash = hashSimFile(file);

end
