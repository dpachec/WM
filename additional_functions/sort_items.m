function [cfg_contrasts] = sort_items(cfg_contrasts)

    oneListIDs = cfg_contrasts.oneListIds; 

    ids = cellfun(@(x) strsplit(string(x)), oneListIDs, 'UniformOutput', false);
    ids0 = double(string(cellfun(@(x) x(3), ids, 'UniformOutput', false)));

    [si1 si2] = sort(ids0);

    cfg_contrasts.oneListIds =  cfg_contrasts.oneListIds(si2);
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(si2, :, :, :);

end