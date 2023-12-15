function [networkRDMs ids2rem]  = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths)

if isfield(cfg_contrasts, 'oneListIds_c')
    oneListIDs = cfg_contrasts.oneListIds_c; 
else
    oneListIDs = cfg_contrasts.oneListIds; 
end


if strcmp(cfg.brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(cfg.brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(cfg.brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

if strcmp (cfg.net2load , 'BLNETi')
    [ACT] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'BNETi')
    [ACT] = load_BNETi(cfg, sessi, subj_ch_fr, paths);    
elseif strcmp (cfg.net2load , 'BFNETi')
    [ACT] = load_BFNETi(cfg, sessi, subj_ch_fr, paths);    
elseif strcmp (cfg.net2load , 'BKNETi')
    [ACT] = load_BKNETi(cfg, sessi, subj_ch_fr, paths);    
elseif strcmp (cfg.net2load , 'BDNETi')
    [ACT] = load_BDNETi(cfg, sessi, subj_ch_fr, paths);        
elseif strcmp (cfg.net2load , 'BLNETe')
    [ACT] = load_BLNETe(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'BLNETeBatchNorm')
    [ACT] = load_BLNETeBatchNorm(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'Alex')
    [ACT] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'AlexEco')
    [ACT] = load_alexECO_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORz')
    [ACT] = load_CORz(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORrt')
    [ACT] = load_CORrt_TV(cfg, sessi, subj_ch_fr, paths);
    %[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORr')
    %[ACT] = load_CORr_TL(cfg, sessi, subj_ch_fr, paths);
    [ACT] = load_CORr_TV(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORrtRELU')
    [ACT] = load_CORrtRELU_activ(cfg, sessi, subj_ch_fr, paths);
elseif strcmp (cfg.net2load , 'CORs')
    [ACT] = load_CORs_TV(cfg, sessi, subj_ch_fr, paths);
elseif length(cfg.net2load) > 6 & ( strcmp (cfg.net2load(1:7), 'Res18-2') | strcmp (cfg.net2load(1:7), 'Res18-4') | strcmp (cfg.net2load(1:7), 'Res18-6') | ...
        strcmp (cfg.net2load(1:7), 'Res18-8') )
    [ACT] = load_Res18_activ(cfg, sessi, paths, oneListIDs);
elseif length(cfg.net2load) > 6 & (strcmp (cfg.net2load(1:7), 'Res34-2') | strcmp (cfg.net2load(1:7), 'Res34-4') | strcmp (cfg.net2load(1:7), 'Res34-6') | ...
        strcmp (cfg.net2load(1:7), 'Res34-8') )
    [ACT] = load_Res34_activ(cfg, sessi, paths, oneListIDs);
elseif strcmp (cfg.net2load , 'CAT')
    [ACT] = load_CATMODEL_activ(cfg, sessi, paths, oneListIDs);
elseif strcmp (cfg.net2load, 'AE-t00')  | strcmp (cfg.net2load, 'AE-t06') | strcmp (cfg.net2load, 'AE-t10')
    net2load = strsplit(cfg.net2load, '-'); 
    net2load = net2load{2};
    [ACT] = load_AE_activ(cfg, sessi, subj_ch_fr, paths, net2load);
else
    %disp ('loading BLnext'); 
    [ACT] = load_blnext(cfg, sessi, paths, oneListIDs);
end


if strcmp(cfg.period, 'M123')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        toSum = double(string(idh(2)));
        ids2(i,:) = idh{12+toSum};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
elseif strcmp(cfg.period, 'M11') 
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(13));
        ids3(i,:) = double(string(idT([1 3])));
        idx = find(~mod(ids3, 10)); 
        ids4 = ids3-10; ids4(idx) = ids3(idx);
    end   
elseif strcmp(cfg.period, 'M12') 
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(14));
        ids3(i,:) = double(string(idT([1 3])));
        idx = find(~mod(ids3, 10)); 
        ids4 = ids3-10; ids4(idx) = ids3(idx);
    end     
elseif strcmp(cfg.period, 'M13') 
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(15));
        ids3(i,:) = double(string(idT([1 3])));
        idx = find(~mod(ids3, 10)); 
        ids4 = ids3-10; ids4(idx) = ids3(idx);
    end 
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
        | strcmp (cfg.net2load , 'CORs')  | strcmp (cfg.net2load , 'CORrtRELU')  | strcmp (cfg.net2load , 'BLNETeBatchNorm') ...
        | strcmp (cfg.net2load , 'AlexEco') | strcmp (cfg.net2load , 'CORr') | strcmp (cfg.net2load , 'AE-t00')  | strcmp (cfg.net2load , 'AE-t10') ...
        | strcmp (cfg.net2load , 'AE-t06') | strcmp (cfg.net2load , 'BKNETi') | strcmp (cfg.net2load , 'BNETi') | strcmp (cfg.net2load , 'BDNETi') ...
        | strcmp (cfg.net2load , 'BFNETi') | strcmp (cfg.net2load , 'CORz') 


    networkRDMs = ACT(:, ids4, ids4);
else
    networkRDMs = ACT; 
end
