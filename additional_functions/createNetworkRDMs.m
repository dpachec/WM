function [networkRDMs ids2rem] = createNetworkRDMs(cfg, oneListIDs, sessi, paths)

net2load = cfg.net2load; 
lays2load = cfg.lays2load; 
brainROI = cfg.brainROI; 
period = cfg.period; 
                                    


if strcmp(brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

if strcmp (net2load , 'RNN')
    [ACT] = load_rnn(lays2load, sessi, subj_ch_fr, paths.activations);%load network if not loaded yet
elseif strcmp (net2load , 'Alex')
    [ACT] = load_alex_activ(lays2load, subji, subj_ch_fr, paths.stim);%load network if not loaded yet
else
    disp ('loading BLnext'); 
    [ACT ids2rem] = load_blnext(net2load, lays2load, sessi, subj_ch_fr, paths.multi_item_activations, brainROI, oneListIDs);
end


if strcmp(period(1), 'M')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        toSum = double(string(idh(2)));
        ids2(i,:) = idh{12+toSum};
    end
    ids3 = double(string(ids2(:,[1 3])));
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
elseif strcmp(period(1), 'E') 
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(3));
        ids3(i,:) = double(string(idT([1 3])));
        idx = find(~mod(ids3, 10)); 
        ids4 = ids3-10; ids4(idx) = ids3(idx);
    end
end



if strcmp (net2load , 'RNN') | strcmp (net2load , 'Alex')
    networkRDMs = ACT(:, ids4, ids4);
else
    networkRDMs = ACT; 
end
