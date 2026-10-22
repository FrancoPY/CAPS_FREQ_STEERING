function hash = hashSimArray(data)
    % Identify numerical RF or medium maps without MAT-file metadata.

    algo = System.Security.Cryptography.SHA256.Create();
    bytes = uint8(algo.ComputeHash(typecast(data(:),'uint8')));
    hash = lower(reshape(dec2hex(bytes,2).',1,[]));
    algo.Dispose();

end
