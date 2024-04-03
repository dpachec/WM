function [cfg_contrasts] = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
% 
 if isfield(cfg_contrasts, 'oneListIds_c')
     cfg_contrasts.oneListIds = cfg_contrasts.oneListIds_c; 
 end

oneListPow = cfg_contrasts.oneListPow;
%oneListIds = cfg_contrasts.oneListIds;
%oneListIds = str2num(cell2mat(oneListIds));
oneListIds = cellfun(@(x) strsplit(x, ' '), cfg_contrasts.oneListIds, 'un', 0);
oneListIds = double(string(cat(1, oneListIds{:})));

if strcmp (zScType, 'sess')
    z2u = 21;
    %z2u = 30;
elseif strcmp(zScType, 'blo')
    z2u = 23;
    %z2u = 32;
end

if acrossTrials
   if ~strcmp(zScType, 'allTrials')
        numSB = max(oneListIds(:, z2u));
   else 
       numSB = 1;
       idSess = ones(length(oneListIds), 1);
   end
   for sessi = 1:numSB
       trltyp = [1 3 5 7];
       if ~strcmp(zScType, 'allTrials')
          idSess = oneListIds(:,z2u) == sessi; 
       end
       for iti = 1:4
           ids = oneListIds(:,1) == trltyp(iti); 
           ids = (ids+idSess ==2); 
           if iti ==4
            ids2 = oneListIds(:,2) ~= 4; 
            ids = ids + ids2 == 2; 
           end

           data = oneListPow(ids, :, :, :);
           mT = mean(data, 1);
           stdT = std(data, [], 1);
           data = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);
           oneListPow(ids, :, :, :) = data; 
       end      
   end 
else
    mT = mean(oneListPow(:,:,:,bline(1):bline(2)),4, 'omitnan');
    stdT = std(oneListPow(:, :, :, bline(1):bline(2)), [], 4, 'omitnan');
    oneListPow = bsxfun(@rdivide, oneListPow - mT, mT); 
end

%disp ('done')

cfg_contrasts.oneListPow = oneListPow;





%%



























                