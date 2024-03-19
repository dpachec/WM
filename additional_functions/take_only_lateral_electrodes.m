function [cfg_contrasts] = take_only_lateral_electrodes(cfg_contrasts)
    

    electrodes = cfg_contrasts.chanNames; 
    ids2keep = zeros(size(electrodes, 1), 1); 

    for chani = 1:size(electrodes, 1)
        
        elecName = electrodes{chani, 5}; 
        elecCoordX = electrodes{chani, 2};
        elecCoordZ = electrodes{chani, 4};

        if (elecCoordX > 35 | elecCoordX <-35) & elecCoordZ > -15
            
            ids2keep(chani, :) = 1; 
        
        end



    end

    ids2keep = logical(ids2keep); 
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, ids2keep, :,:); 
    cfg_contrasts.chanNames = cfg_contrasts.chanNames(ids2keep, :); 




end