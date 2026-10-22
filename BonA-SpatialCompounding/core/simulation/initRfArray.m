function initRfArray(file,dims)

    % Create an extensible double RF array in an initialized v7.3 MAT file
    % Chunk complete time/element frames so individual shots can resume cheaply
    h5create(file,'/rf',[dims(1),dims(2),Inf,Inf],'Datatype','double','ChunkSize',[dims(1),dims(2),1,1]);
    h5writeatt(file,'/rf','MATLAB_class','double');
    h5write(file,'/rf',zeros(dims(1),dims(2)),[1,1,dims(3),dims(4)],[dims(1),dims(2),1,1]);

end
