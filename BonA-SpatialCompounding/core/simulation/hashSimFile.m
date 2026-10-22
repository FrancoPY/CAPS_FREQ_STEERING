function hash = hashSimFile(file)
    % Identify saved RF or density content independently of its folder.

    % Hash the file in a stream without loading the acquisition into memory
    stream = System.IO.File.OpenRead(file);
    cleanup = onCleanup(@()stream.Close());
    algo = System.Security.Cryptography.SHA256.Create();
    bytes = uint8(algo.ComputeHash(stream));
    hash = lower(reshape(dec2hex(bytes,2).',1,[]));
    algo.Dispose();

end
