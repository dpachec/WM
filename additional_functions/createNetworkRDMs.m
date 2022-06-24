function [networkRDMS] = createNetworkRDMs(oneListIDs, lays2load, brainROI, subji)

if strcmp(brainROI, 'vvs')  subj_ch_fr = 17; end %%VVS
if strcmp(brainROI, 'pfc')  subj_ch_fr = 7; end %%PFC
if strcmp(brainROI, 'hipp') subj_ch_fr = 8; end %%HIPP

[ACT] = load_rnn(lays2load, subji, subj_ch_fr);%load network if not loaded yet


for i = 1:length(oneListIDs)
    idh = strsplit(oneListIDs{i});
    toSum = double(string(idh(2)));
    ids2(i,:) = idh{12+toSum};
end

ids3 = double(string(ids2(:,[1 3])));
idx = find(~mod(ids3, 10)); 
ids4 = ids3-10; ids4(idx) = ids3(idx);

networkRDMS = ACT(:, ids4, ids4); 
