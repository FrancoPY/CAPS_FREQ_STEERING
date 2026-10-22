function [out,state,txIdx,changed] = syncSimRf(file,key,densDir,densHash,txKey)
    % Resume RF(time,element,frame,shot), allowing new frames and shots.

    out = matfile(file,'Writable',true);
    nFrame = numel(densHash);
    oldState = struct();
    if ismember('simState',who('-file',file))
        saved = load(file,'simState');
        state = saved.simState;
        oldState = state;
        assert(isequaln(state.key,key),'The acoustic settings changed. Delete the incorrect RF before rerunning this simulation.');
    else
        state.key = key;
        state.densDir = densDir;
        state.densHash = {};
        state.txKey = zeros(size(txKey,1),0);
        state.done = false(0,0);
    end

    % Append new shot identities without changing existing RF indices
    [known,txIdx] = ismember(txKey.',state.txKey.','rows');
    txIdx = reshape(txIdx,1,[]);
    for shot = find(~known).'
        state.txKey(:,end + 1) = txKey(:,shot);
        txIdx(shot) = size(state.txKey,2);
    end
    nFrame = max(nFrame,size(state.done,1));
    done = false(nFrame,size(state.txKey,2));
    done(1:size(state.done,1),1:size(state.done,2)) = state.done;
    state.done = done;

    % Invalidate every shot that used a replaced density realization
    for frame = 1:numel(densHash)
        if frame > numel(state.densHash) || ~strcmp(state.densHash{frame},densHash{frame})
            state.done(frame,:) = false;
        end
    end
    state.densDir = densDir;
    state.densHash(1:numel(densHash)) = densHash;
    changed = ~isequaln(oldState,state);
    if changed
        out.simState = state;
    end

end
