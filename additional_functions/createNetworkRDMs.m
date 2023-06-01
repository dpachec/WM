function [networkRDMs ids2rem]  = createNetworkRDMs(cfg, oneListIDs, sessi, paths)
         


if strcmp(cfg.brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(cfg.brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(cfg.brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

if strcmp (cfg.net2load , 'BLNETi')
    [ACT] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'BLNETe')
    [ACT] = load_BLNETe(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'Alex')
    [ACT] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORrt')
    [ACT] = load_COR_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORs')
    [ACT] = load_CORs_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'Res18-2') | strcmp (cfg.net2load , 'Res18-4') | strcmp (cfg.net2load , 'Res18-6') | strcmp (cfg.net2load , 'Res18-8') 
    [ACT] = load_Res18_activ(cfg, sessi, paths, oneListIDs);
elseif strcmp (cfg.net2load , 'CAT')
    [ACT] = load_CATMODEL_activ(cfg, sessi, paths, oneListIDs);
else
    %disp ('loading BLnext'); 
    [ACT] = load_blnext(cfg, sessi, paths, oneListIDs);
end


if strcmp(cfg.period, 'MALL1')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        ids2(i,:) = idh{13};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
end
if strcmp(cfg.period, 'MALL2')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        ids2(i,:) = idh{14};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
end
if strcmp(cfg.period, 'MALL3')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        ids2(i,:) = idh{15};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
end

if strcmp(cfg.period(1), 'M')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        toSum = double(string(idh(2)));
        ids2(i,:) = idh{12+toSum};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
elseif strcmp(cfg.period(1), 'E') 
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(3));
        ids3(i,:) = double(string(idT([1 3])));
        idx = find(~mod(ids3, 10)); 
        ids4 = ids3-10; ids4(idx) = ids3(idx);
    end
end



if strcmp (cfg.net2load , 'BLNETi') | strcmp (cfg.net2load , 'BLNETe') | strcmp (cfg.net2load , 'Alex') | strcmp (cfg.net2load , 'CORrt') ...
        | strcmp (cfg.net2load , 'CORs') 
    networkRDMs = ACT(:, ids4, ids4);
else
    networkRDMs = ACT; 
end
