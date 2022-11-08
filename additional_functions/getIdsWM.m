function[ids] = getIdsWM(period, cfg_contrasts)

    if strcmp(period, 'M')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(period, 'E1')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '1'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(period, 'E2')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '3'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(period, 'E3')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(period, 'E123')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '1') | strcmp(x(1), '3') | strcmp(x(1), '5'), cfg_contrasts.oneListIds_c, 'un', 0));
    end

end