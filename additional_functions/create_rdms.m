function [out_rdms] = create_rdms (cfg_contrasts, f, it, win_width, mf, avTW)
 
if strcmp(it, 'all')
    if isfield(cfg_contrasts, 'oneListIds_enc')
        oneListIds              =       cfg_contrasts.oneListIds_enc;
        oneListPow              =       cfg_contrasts.oneListPow_enc;
    else
        oneListIds              =       cfg_contrasts.oneListIds_c;
        oneListPow              =       cfg_contrasts.oneListPow;
    end
else
    if isfield(cfg_contrasts, 'oneListIds_enc')
        oneListIds              =       cfg_contrasts.oneListIds_maint;
        oneListPow              =       cfg_contrasts.oneListPow_maint;
    else
        oneListIds              =       cfg_contrasts.oneListIds_c;
        oneListPow              =       cfg_contrasts.oneListPow;
    end
end
    
chanNames               =       cfg_contrasts.chanNames;


countALL_EE = 1;
ALL_EE = [];

clear new_all_ee
for i = 1:length(oneListIds)
    evei = strsplit(oneListIds{i});
    if strcmp(it, 'all') | strcmp(it, 'eCue') 
       if strcmp(evei(1), '1') | strcmp(evei(1), '3') | strcmp(evei(1), '5') | strcmp(evei(1), '*')
            %disp (['ALL_EE> ' oneListIds{i}]);   
            new_all_ee{countALL_EE} = [i];
            countALL_EE = countALL_EE+1;
       end

    else
        if strcmp(evei(1), it) | strcmp(evei(1), '7') | strcmp(evei(1), '7*') | strcmp(evei(1), '-7')  | strcmp(evei(1), '9') | strcmp(evei(1), '8') 
           %disp (['ALL_EE> ' oneListIds{i}]);   
           new_all_ee{countALL_EE} = [i];
           countALL_EE = countALL_EE+1;
        end
    end
end

new_all_ee = cell2mat(new_all_ee)';

if ndims(oneListPow) > 3
    ALL_EE = oneListPow(new_all_ee, :, f, :);
else
    ALL_EE = oneListPow(new_all_ee, f, :);
end

ALL_EE_IDS = oneListIds(new_all_ee);

% % % % take id of cued item
if strcmp(it, '7')  | strcmp(it, '9')  % 7 = maintenance 123

    ids = cellfun(@(x) strsplit(x), ALL_EE_IDS, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(2), ids, 'UniformOutput', false);
    ids1 = cellfun(@(x) x(1), ids, 'UniformOutput', false);
    x2check  = char(string(ids0));
    x2check2  = char(string(ids1));
    idx2 = x2check ~= '4' & ( x2check2 == '7' | x2check2 == '7*' | x2check2 == '-7' | x2check2 == '9'  );
    ALL_EE_IDS = ALL_EE_IDS(logical(sum(idx2, 2)));

    clear ALL_EE_IDS2;
    ALL_EE_IDS2 = ALL_EE_IDS;
    for i = 1:length(ALL_EE_IDS)
        idh = strsplit(ALL_EE_IDS{i});
        toSum = double(string(idh(2)));
        idh{3} = idh{12+toSum};
        ALL_EE_IDS2{i}=  join(idh, '    ');
    end
    ALL_EE_IDS = ALL_EE_IDS2;
    ALL_EE = ALL_EE(logical(sum(idx2, 2)),:,:,:);
end

if strcmp(it, '8') % 8 = maintenance-all trials

    ids = cellfun(@(x) strsplit(x), ALL_EE_IDS, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(2), ids, 'UniformOutput', false);
    ids1 = cellfun(@(x) x(1), ids, 'UniformOutput', false);
    x2check  = char(string(ids0));
    x2check2  = char(string(ids1));
    idx2 = x2check == '4' & x2check2 == '7';
    ALL_EE_IDS = ALL_EE_IDS(idx2);

    clear ALL_EE_IDS2;
    ALL_EE_IDS2 = ALL_EE_IDS;
    for i = 1:length(ALL_EE_IDS)
        idh = strsplit(ALL_EE_IDS{i});
        idh{3} = [idh{13} '_' idh{14} '_' idh{15}];
        ALL_EE_IDS2{i}=  join(idh, '    ');
    end
    ALL_EE_IDS = ALL_EE_IDS2;
    ALL_EE = ALL_EE(idx2,:,:,:);
end


if strcmp(it, 'eCue') 
    ids = cellfun(@(x) strsplit(x), ALL_EE_IDS, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(2), ids, 'UniformOutput', false);
    ids1 = cellfun(@(x) x(1), ids, 'UniformOutput', false);
    x2check  = char(string(ids0));
    idx2 = x2check ~= '4';
    ALL_EE_IDS = ALL_EE_IDS(idx2);

    for i = 1:length(ALL_EE_IDS)
        idh = strsplit(ALL_EE_IDS{i});
        toSum = double(string(idh(2)));
        if idh{3} == idh{12+toSum}
            idx3(i) = 1; 
        else
            idx3(i) = 0; 
        end
    end
    idx3 = logical(idx3);
    ALL_EE_IDS = ALL_EE_IDS(idx3);
    ALL_EE = ALL_EE(idx3,:,:,:);
end


if ndims(ALL_EE) > 3
    nTrials = size(ALL_EE, 1); nChans = size(ALL_EE, 2); nFreq = size(ALL_EE, 3); nTimes = size(ALL_EE, 4); 
    ALLER = reshape (ALL_EE, [nTrials, nChans*nFreq, nTimes]);
else
    nTrials = size(ALL_EE, 1); nFreq = size(ALL_EE, 2); nTimes = size(ALL_EE, 3); 
    ALLER = ALL_EE; 
end

%covF = mean(ALLER, 3, 'omitnan');
%sigma =  covdiag(covF); % estimate covariance matrix 

clear rdm
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
parfor timei = 1:bins 
    %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    if avTW
        x = mean(ALLER(:, :, timeBins), 3, 'omitnan');
    else
        x = ALLER(:, :, timeBins);
        x = reshape(x, nTrials, []); 
    end
    rdm(:, :, timei) = corr(x', 'type', 's');
end


out_rdms.rdm = rdm;
out_rdms.ids = ALL_EE_IDS; 















