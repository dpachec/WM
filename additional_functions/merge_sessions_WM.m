function [cfg_contrasts] = merge_sessions_WM(d2m, sb)


nContr = length(d2m);

for ctr = 1:nContr
    loadedC2M(ctr) = load(d2m{ctr});  % stores all blocks and sesions from 1 subject    
end

if isfield(loadedC2M(1).cfg_contrasts, 'oneListPow')

    clear ids_sess oneListIds oLP
    for i = 1:nContr
        ids_sess{i,:} = cellfun(@(x) [x ['    ' num2str(sb(i, 1))] ['    ' num2str(sb(i, 2))]   ['    ' num2str(i)] ], ...
                    loadedC2M(i).cfg_contrasts.oneListIds_c, 'UniformOutput', false);
        oLP{i} =  loadedC2M(i).cfg_contrasts.oneListPow; 
    end

    oneListIds = cat(1, ids_sess{:});
    oneListPow = cat(1, oLP{:});

    cfg_contrasts= []; 
    cfg_contrasts.oneListIds = oneListIds; 
    cfg_contrasts.oneListPow = oneListPow; 
    cfg_contrasts.chanNames =  loadedC2M(i).cfg_contrasts.chanNames; 

    
    
else % for the traces only
   
    clear ids_sess oneListIds oLP
    for i = 1:nContr
        ids_sess{i,:} = cellfun(@(x) [x ['    ' num2str(sb(i, 1))] ['    ' num2str(sb(i, 2))]   ['    ' num2str(i)] ], ...
                    loadedC2M(i).cfg_contrasts.oneListIds, 'UniformOutput', false);
        oLT{i} =  loadedC2M(i).cfg_contrasts.oneListTraces; 
    end

    oneListIds = cat(1, ids_sess{:});
    oneListTraces = cat(3, oLT{:});

    cfg_contrasts= []; 
    cfg_contrasts.oneListIds = oneListIds; 
    cfg_contrasts.oneListTraces = oneListTraces; 
    cfg_contrasts.chanNames =  loadedC2M(i).cfg_contrasts.chanNames; 
    
    
end
