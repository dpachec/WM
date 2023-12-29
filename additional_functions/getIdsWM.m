function[cfg_contrasts] = getIdsWM(period, cfg_contrasts)

    %oneListIds = cellfun(@(x) strsplit(x, ' '), cfg_contrasts.oneListIds_c, 'un', 0);
    %oneListIds = double(string(cat(1, oneListIds{:})));
    oneListIds = str2num(cell2mat(cfg_contrasts.oneListIds_c));
    oneListIdsa = cellfun(@(x) strsplit(x, ' '), cfg_contrasts.oneListIds_c, 'un', 0);
    

    isCued = ((oneListIds(:, 1) == 1 & oneListIds(:, 2) == 1) | (oneListIds(:, 1) == 3 & oneListIds(:, 2) == 2) | (oneListIds(:, 1) == 5 & oneListIds(:, 2) == 3) ); 

    if strcmp(period, 'M123')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4; 
    elseif strcmp(period, 'MALL')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) == 4; 
    elseif strcmp(period, 'M123CI')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4 & oneListIds(:, 8) == 1;
    elseif strcmp(period, 'M123II')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4 & oneListIds(:, 8) == 0;
    elseif strcmp(period, 'M123CC')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4 & oneListIds(:, 9) == 1;
    elseif strcmp(period, 'M123IC')
        ids = oneListIds(:, 1) == 7 & oneListIds(:, 2) ~= 4 & oneListIds(:, 9) == 0;        
    elseif strcmp(period, 'E123')
        ids = oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5; 
    elseif strcmp(period, 'E123CI')
        ids = isCued & oneListIds(:, 8) == 1; 
    elseif strcmp(period, 'E123II')
        ids = isCued & oneListIds(:, 8) == 0; 
    elseif strcmp(period, 'E123CC')
        ids = isCued & oneListIds(:, 9) == 1; 
    elseif strcmp(period, 'E123IC')
        ids = isCued & oneListIds(:, 9) == 0;         
    elseif strcmp(period, 'M11') 
        ids = oneListIds(:, 1) == 5;  
        for i = 1:length(oneListIds)
            oneListIdsa{i}(3)=oneListIdsa{i}(13); 
        end
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:}));
    elseif strcmp(period, 'M12') 
         ids = oneListIds(:, 1) == 5;  
        for i = 1:length(oneListIds)
            oneListIdsa{i}(3)=oneListIdsa{i}(14); 
        end
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:}));
    elseif strcmp(period, 'M13') 
        ids = oneListIds(:, 1) == 5;  
        for i = 1:length(oneListIds)
            oneListIdsa{i}(3)=oneListIdsa{i}(15); 
        end
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:}));
    elseif strcmp(period, 'E123CAT1') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 1;
    elseif strcmp(period, 'E123CAT2') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 2;
    elseif strcmp(period, 'E123CAT3') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 3;
    elseif strcmp(period, 'E123CAT4') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 4;
    elseif strcmp(period, 'E123CAT5') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 5;
    elseif strcmp(period, 'E123CAT6') 
         ids = (oneListIds(:, 1) == 1 | oneListIds(:, 1) == 3 | oneListIds(:, 1) == 5) & ...
                floor(oneListIds(:, 3) / 100) == 6;
         
    elseif strcmp(period, 'M123CAT1') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 1;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 
        
    elseif strcmp(period, 'M123CAT2') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 2;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 
        
    elseif strcmp(period, 'M123CAT3') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 3;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 
        
    elseif strcmp(period, 'M123CAT4') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 4;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 

    elseif strcmp(period, 'M123CAT5') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 5;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 
        
    elseif strcmp(period, 'M123CAT6') 
        for i = 1:length(oneListIds)
            cue = double(string(oneListIdsa{i}(2)));
            oneListIdsa{i}(3)=oneListIdsa{i}(12+cue); 
        end
        oneListIdsb = double(string(cat(1, oneListIdsa{:})));
        ids = ( oneListIdsb(:, 1) == 7 & oneListIdsb(:, 2) ~= 4 ) & floor(oneListIdsb(:, 3) / 100) == 6;
        cfg_contrasts.oneListIds_c = join(cat(1, oneListIdsa{:})); 
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