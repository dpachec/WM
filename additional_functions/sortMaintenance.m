function [ids neuralRDMs] = sortMaintenance(ids, neuralRDMs)

    ids0 = cellfun(@(x) strsplit(x), ids, 'UniformOutput', false);
    cues = double(string(cellfun(@(x) x(2), ids0, 'UniformOutput', false)));
    
    
    clear idsF ids23
    for cuei = 1:length(cues)
        idsF(cuei,:) = double(string(ids0{cuei}(12+cues(cuei))));
        ids0{cuei}(3) = ids0{cuei}(12+cues(cuei));
        ids23{cuei,:} = char(join(ids0{cuei}, '  '));
    end

    [id1 id2] = sort(idsF);
    ids = ids(id2);
    neuralRDMs = neuralRDMs(id2, id2,:,:); 












end