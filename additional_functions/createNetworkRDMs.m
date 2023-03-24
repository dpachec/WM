function [networkRDMs ids2rem]  = createNetworkRDMs(cfg, oneListIDs, sessi, paths)
         


if strcmp(cfg.brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(cfg.brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(cfg.brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

if strcmp (cfg.net2load , 'RNN')
    %[ACT] = load_rnn(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
    [ACT] = load_rnn_eco(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
elseif strcmp (cfg.net2load , 'Alex')
    [ACT] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
else
    disp ('loading BLnext'); 
    [ACT] = load_blnext(cfg, sessi, paths, oneListIDs);
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



if strcmp (cfg.net2load , 'RNN') | strcmp (cfg.net2load , 'Alex')
    networkRDMs = ACT(:, ids4, ids4);
else
    networkRDMs = ACT; 
end
