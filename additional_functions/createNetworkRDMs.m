function [networkRDMS] = createNetworkRDMs(oneListIDs, net2load, lays2load, brainROI, subji, paths, period)

if strcmp(brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

if strcmp (net2load , 'RNN')
    [ACT] = load_rnn(lays2load, subji, subj_ch_fr, paths.activations);%load network if not loaded yet
else
    [ACT] = load_alex_activ(lays2load, subji, subj_ch_fr, paths.stim);%load network if not loaded yet
end

if strcmp(period, 'M')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        toSum = double(string(idh(2)));
        ids2(i,:) = idh{12+toSum};
    end
    ids3 = double(string(ids2(:,[1 3])));
elseif strcmp(period, 'E')
    for i = 1:length(oneListIDs)
        idh = strsplit(oneListIDs{i});
        idT = char(idh(3));
        ids3(i,:) = double(string(idT([1 3])));
        
    end
end

idx = find(~mod(ids3, 10)); 
ids4 = ids3-10; ids4(idx) = ids3(idx);

networkRDMS = ACT(:, ids4, ids4); 
