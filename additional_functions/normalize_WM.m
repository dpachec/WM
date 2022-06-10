function [cfg_contrasts] = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)

oneListPow = cfg_contrasts.oneListPow;
oneListIds = cfg_contrasts.oneListIds_c;
oneListIds = str2num(cell2mat(oneListIds));

if strcmp (zScType, 'sess')
    z2u = 21;
elseif strcmp(zScType, 'blo')
    z2u = 23;
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



























                