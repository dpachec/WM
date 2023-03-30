function[cfg_contrasts] = getIdsWM(period, cfg_contrasts)

    if strcmp(period, 'M123')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'MALL')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'MALL1') | strcmp(period, 'MALL2') | strcmp(period, 'MALL3') 
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'M1')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'M21')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '1'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E1')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '1'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E2')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '3'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E3')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E123')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '1') | strcmp(x(1), '3') | strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E41') % HACK ONLY FOR BL NEXT NOW 
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '1') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E42') % HACK ONLY FOR BL NEXT NOW 
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '3') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E43') % HACK ONLY FOR BL NEXT NOW 
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '5') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(period, 'E44') % HACK ONLY FOR BL NEXT NOW 
        ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '3') | strcmp(x(1), '5' )) & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end



    cfg_contrasts.oneListTraces   = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 

end