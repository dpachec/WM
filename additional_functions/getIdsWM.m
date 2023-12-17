function[cfg_contrasts] = getIdsWM(period, cfg_contrasts)

    oneListIds = cellfun(@(x) strsplit(x, ' '), cfg_contrasts.oneListIds_c, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    if strcmp(period, 'M123')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4; 
    elseif strcmp(period, 'MALL')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) == 4; 
    elseif strcmp(period, 'M123C')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) == 4 & oneListIds(:, 8) == 1;
    elseif strcmp(period, 'M123I')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) == 4 & oneListIds(:, 8) == 0;
    elseif strcmp(period, 'E123')
        ids = oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5; 
    elseif strcmp(period, 'E123C')
        ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & oneListIds(:, 8) == 1; 
    elseif strcmp(period, 'E123I')
        ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & oneListIds(:, 8) == 0; 
    elseif strcmp(period, 'M11') | strcmp(period, 'M12') | strcmp(period, 'M13')
        ids = oneListIds(:, 1) == 5; 
    end

    cfg_contrasts.oneListTraces   = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds      = cfg_contrasts.oneListIds_c(ids); 

    if isfield(cfg_contrasts, 'oneListIds_c')
        cfg_contrasts = rmfield(cfg_contrasts, 'oneListIds_c'); 
    end

     
%     elseif strcmp(period, 'MALL1') | strcmp(period, 'MALL2') | strcmp(period, 'MALL3') 
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));

%     elseif strcmp(period, 'M21')
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '1'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E1')
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '1'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E2')
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '3'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E3')
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
%     
%     elseif strcmp(period, 'E41') % HACK ONLY FOR BL NEXT NOW 
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '1') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E42') % HACK ONLY FOR BL NEXT NOW 
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '3') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E43') % HACK ONLY FOR BL NEXT NOW 
%         ids = cell2mat(cellfun(@(x) strcmp(x(1), '5') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
%     elseif strcmp(period, 'E44') % HACK ONLY FOR BL NEXT NOW 
%         ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '3') | strcmp(x(1), '5' )) & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
%     end

    




end