%% PLV inter area (PFV-VVS)
%% particular time period and band

clearvars 


f2u = [3 8];
tP = 2001:2800;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)


    t_vvs = squeeze(c_vvs.oneListTraces(:,tP,:));
    EEG = []; EEG.data = t_vvs; EEG.trials = size(t_vvs, 3); EEG.srate=1000; EEG.nbchan=size(t_vvs, 1); EEG.pnts=size(t_vvs,2);EEG.event=[];
    EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
    data_vvs        = squeeze(EEG.data); 
    dataHA_vvs      = angle(hilbert(data_vvs));
    t_pfc = squeeze(c_pfc.oneListTraces(:,tP,:));
    EEG = []; EEG.data = t_pfc; EEG.trials = size(t_pfc, 3); EEG.srate=1000; EEG.nbchan=size(t_pfc, 1); EEG.pnts=size(t_pfc,2);EEG.event=[];
    EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
    data_pfc        = squeeze(EEG.data); 
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        parfor chanj = 1:size(c_pfc.chanNames,1)
            diffPha = angle(hilbert(squeeze(data_vvs(chani, :,:)))) - angle(hilbert(squeeze(data_pfc(chanj, :,:))));
            PLV2U(:, chani, chanj) = abs(mean(exp(1i*(diffPha))));
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


save([paths.results.PLV 'time_PLV_ALL_' num2str(f2u(1)) '-' num2str(f2u(2)) '_' num2str(tP(1)) '-' num2str(tP(end))], 'PLV_ALL');




toc

%% process PLV_ALL (across time)
clear SI_TR MI_TR

for subji = 1:10

    allPLV = PLV_ALL{subji, 1};
    allIDs = PLV_ALL{subji, 2};
    ids0 = cellfun(@(x) strsplit(x), allIDs, 'un', 0);
    ids1 = cell2mat(cellfun(@(x) double(string(x(1))), ids0, 'un', 0));
    ids2 = cell2mat(cellfun(@(x) double(string(x(2))), ids0, 'un', 0));

    
    SI_TR{subji,1} = allPLV(ids1 == 7 & ids2 ~= 4,:,:,:); 
    MI_TR{subji,1} = allPLV(ids1 == 7 & ids2 == 4,:,:,:); 

end

%% plot 

d2pSI = cellfun(@(x) squeeze(mean(mean(mean(x)))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(mean(x)))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 
%c1 = logit(c1)
%c2 = logit(c2)
md2pSI = mean(c1)
md2pMI = mean(c2)

[h p ci ts] = ttest(c1, c2)

%% 2Bar 
data.data = [c1 c2]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.32 .45] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


%% Across trials

clearvars 

f2u = [3 8];

paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)
    
    clear PLV2U_SI PLV2U_MI
    for chani = 1:size(c_vvs.chanNames,1)
        parfor chanj = 1:size(c_pfc.chanNames,1)
            id2u = cellfun(@(x) strsplit(x, ' '), c_vvs.oneListIds_c, 'un', 0);
            id2u = cat(1, id2u{:});
            id2u_SI = strcmp(id2u(:, 1), '7') & ~strcmp(id2u(:, 2), '4'); 
            id2u_MI = strcmp(id2u(:, 1), '7') & strcmp(id2u(:, 2), '4'); 

            % % 
            t_vvsSI = squeeze(c_vvs.oneListTraces(chani,:,id2u_SI));
            t_pfcSI = squeeze(c_pfc.oneListTraces(chanj,:,id2u_SI));
            t_vvsMI = squeeze(c_vvs.oneListTraces(chani,:,id2u_MI));
            t_pfcMI = squeeze(c_pfc.oneListTraces(chanj,:,id2u_MI));            
            
            %match trial numbers
            % get condition with more trials 
            nTrC1 = size(t_vvsSI,2); nTrC2 = size(t_vvsMI, 2); 
            id2uSample = randperm(nTrC2);
            if nTrC1 > nTrC2
                t_vvsSI = t_vvsSI(:, id2uSample);
                t_pfcSI = t_pfcSI(:, id2uSample);
            elseif nTrC2 > nTrC1
                t_vvsMI = t_vvsMI(:, id2uSample);
                t_pfcMI = t_pfcMI(:, id2uSample);
            end

            % % % SINGLE ITEM TRIALS
            EEG = []; 
            EEG.data    = t_vvsSI;
            EEG.trials  = size(t_vvsSI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_vvsSI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_vvs        = squeeze(EEG.data); 
            dataHA_vvs      = angle(hilbert(data_vvs));
            EEG = []; 
            EEG.data    = t_pfcSI;
            EEG.trials  = size(t_pfcSI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_pfcSI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_pfc        = squeeze(EEG.data); 
            dataHA_pfc      = angle(hilbert(data_pfc));
            diffPha = dataHA_vvs - dataHA_pfc;
            PLV2U_SI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2));

            % % % % MULTI ITEM TRIALS
            EEG = []; 
            EEG.data    = t_vvsMI;
            EEG.trials  = size(t_vvsMI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_vvsMI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_vvs        = squeeze(EEG.data); 
            dataHA_vvs      = angle(hilbert(data_vvs));
            EEG = []; 
            EEG.data    = t_pfcMI;
            EEG.trials  = size(t_pfcMI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_pfcMI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_pfc        = squeeze(EEG.data); 
            dataHA_pfc      = angle(hilbert(data_pfc));
            diffPha = dataHA_vvs - dataHA_pfc;
            PLV2U_MI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2));

        end
    end
    PLV_ALL{subji,1} = PLV2U_SI; %pfc or vvs 
    PLV_ALL{subji,2} = PLV2U_MI; %pfc or vvs 
end


save([paths.results.PLV 'trials_PLV_ALL'], 'PLV_ALL');



toc




%% Process across trials 

SI_TR = PLV_ALL(:, 1);
MI_TR = PLV_ALL(:, 2);

d2pSI = cellfun(@(x) squeeze(mean(mean(x))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(x))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 
%c1 = logit(c1)
%c2 = logit(c2)
md2pSI = mean(c1)
md2pMI = mean(c2)

[h p ci ts] = ttest(c1, c2)

%% Plot
d2SI = cell2mat(d2pSI')';
d2MI = cell2mat(d2pMI')';

mSI = squeeze(mean(d2SI));
mMI = squeeze(mean(d2MI));

hb = h; hb(hb==0) = nan; hb(hb==1) = 0.1; 
figure;
plot(mSI, 'r', 'LineWidth', 2); hold on
plot(mMI, 'b', 'LineWidth', 2)
plot(hb, 'LineWidth', 4)

%% check for the time period of network fit

c1R = mean(c1(:, 2400:3200), 2)
c2R = mean(c2(:, 2400:3200), 2)
[h p ci ts] = ttest(c1R, c2R)


%% 2Bar 
data.data = [c1R c2R]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.07 .24] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

%% Across time frequency resolved 
clearvars 


f2u = [3 8];
tP = 2001:2800;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)

    clear PLV2U
    for freqi = [3:29 30:5:150]
        t_vvs = squeeze(c_vvs.oneListTraces(:,tP,:));
        EEG = []; EEG.data = t_vvs; EEG.trials = size(t_vvs, 3); EEG.srate=1000; EEG.nbchan=size(t_vvs, 1); EEG.pnts=size(t_vvs,2);EEG.event=[];
        EEG         = pop_eegfiltnew (EEG, freqi,freqi);
        data_vvs        = squeeze(EEG.data); 
        dataHA_vvs      = angle(hilbert(data_vvs));
        t_pfc = squeeze(c_pfc.oneListTraces(:,tP,:));
        EEG = []; EEG.data = t_pfc; EEG.trials = size(t_pfc, 3); EEG.srate=1000; EEG.nbchan=size(t_pfc, 1); EEG.pnts=size(t_pfc,2);EEG.event=[];
        EEG         = pop_eegfiltnew (EEG, freqi,freqi);
        data_pfc        = squeeze(EEG.data); 
        for chani = 1:size(c_vvs.chanNames,1)
            for chanj = 1:size(c_pfc.chanNames,1)
                diffPha = angle(hilbert(squeeze(data_vvs(chani, :,:)))) - angle(hilbert(squeeze(data_pfc(chanj, :,:))));
                PLV2U(:, freqi, chani, chanj) = abs(mean(exp(1i*(diffPha))));
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


save([paths.results.PLV 'time_PLV_ALL_FR_' num2str(tP(1)) '-' num2str(tP(end))], 'PLV_ALL');




toc



%% TIME resolved in a particular band

clearvars 

f2u = [3 8]

currentF = pwd;
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

nTimes = 5000; 
win_width = 500; 
mf = 100; 
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    %restrict time 
    c_vvs.oneListTraces = c_vvs.oneListTraces(:, 1001:6000,:);
    c_pfc.oneListTraces = c_pfc.oneListTraces(:, 1001:6000,:);
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            parfor triali = 1:size(c_vvs.oneListTraces,3)
                id2u = strsplit(c_vvs.oneListIds_c{triali});
                t_vvs = squeeze(c_vvs.oneListTraces(chani,:,triali));
                t_pfc = squeeze(c_pfc.oneListTraces(chanj,:,triali));
                EEG_vvs.data    = t_vvs;
                EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                EEG_vvs         = pop_eegfiltnew (EEG_vvs, f2u(1),f2u(2));
                data_vvs        = squeeze(EEG_vvs.data); 
                dataHA_vvs      = angle(hilbert(data_vvs));
                EEG_pfc.data    = t_pfc;
                EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                EEG_pfc         = pop_eegfiltnew (EEG_pfc, f2u(1),f2u(2));
                data_pfc        = squeeze(EEG_pfc.data); 
                dataHA_pfc      = angle(hilbert(data_pfc));
                diffPha = dataHA_vvs - dataHA_pfc;
                for timei = 1:bins 
                    timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                    diffPhaBIN = diffPha(timeBins);                                    
                    PLV2U(triali, chani, chanj, timei, :) = abs(mean(exp(1i*(diffPhaBIN))));
                end
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


cd (currentF)
save('PLV_ALL', 'PLV_ALL');



toc

%% Frequency resolved for cluster period

clearvars 

currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc';


tP = 2401:3200; 

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    %restrict time 
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            for triali = 1:size(c_vvs.oneListTraces,3)
                triali
                id2u = strsplit(c_vvs.oneListIds_c{triali});
                t_vvs = squeeze(c_vvs.oneListTraces(chani,tP,triali));
                t_pfc = squeeze(c_pfc.oneListTraces(chanj,tP,triali));
                for freqi = 1:54
                    EEG_vvs.data    = t_vvs;
                    EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                    EEG_vvs         = pop_eegfiltnew (EEG_vvs, freqi,freqi);
                    data_vvs        = squeeze(EEG_vvs.data); 
                    dataHA_vvs      = angle(hilbert(data_vvs));
                    EEG_pfc.data    = t_pfc;
                    EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                    EEG_pfc         = pop_eegfiltnew (EEG_pfc, freqi,freqi);
                    data_pfc        = squeeze(EEG_pfc.data); 
                    dataHA_pfc      = angle(hilbert(data_pfc));
                    diffPha = dataHA_vvs - dataHA_pfc;
                    PLV2U(triali, chani, chanj, freqi, :) = abs(mean(exp(1i*(diffPha))));
                end
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


cd (currentF)
save('PLV_ALL', 'PLV_ALL');


toc


%% convert to logit
clear PLV_ALL_LG
for subji = 1:10
    PLV2U = PLV_ALL{subji, 1};
    for triali = 1:size(PLV2U, 1)
        for chani = 1:size(PLV2U, 2)
            for chanj = 1:size(PLV2U, 3)

                x = PLV2U(triali, chani, chanj, :);
                xLG = logit(x);
                PLV2U_LG(triali, chani, chanj, :) = xLG;

            end
        end

    end
    PLV_ALL_LG{subji,1} = PLV2U_LG;
    PLV_ALL_LG{subji,2} = PLV_ALL{subji, 2};
end



%% process PLV_ALL
clear SI_TR MI_TR

for subji = 1:10

    allPLV = PLV_ALL{subji, 1};
    allIDs = PLV_ALL{subji, 2};
    ids0 = cellfun(@(x) strsplit(x), allIDs, 'un', 0)
    ids1 = cell2mat(cellfun(@(x) double(string(x(1))), ids0, 'un', 0));
    ids2 = cell2mat(cellfun(@(x) double(string(x(2))), ids0, 'un', 0));

    
    SI_TR{subji,1} = allPLV(ids1 == 7 & ids2 == 4,:,:,:); 
    MI_TR{subji,1} = allPLV(ids1 == 7 & ids2 ~= 4,:,:,:); 

end

%% plot 

d2pSI = cellfun(@(x) squeeze(mean(mean(mean(x)))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(mean(x)))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 
md2pSI = mean(c1)
md2pMI = mean(c2)

plot(md2pSI, 'r'); hold on; 
plot(md2pMI, 'b')

[h p ci ts] = ttest(c1, c2)





%% TIME FREQUENCY reSOLVED

clearvars 

currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc';


pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

nTimes = 5000; 
win_width = 500; 
mf = 100; 
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    %restrict time 
    c_vvs.oneListTraces = c_vvs.oneListTraces(:, 1001:6000,:);
    c_pfc.oneListTraces = c_pfc.oneListTraces(:, 1001:6000,:);

    count = 0;
    for triali = 1:size(c_vvs.oneListTraces,3)
        id2u = strsplit(c_vvs.oneListIds_c{triali});
        if strcmp(id2u{1},'7') & strcmp(id2u{2},'4') 
            count = count+1
            for chani = 1:size(c_vvs.chanNames,1)
                for chanj = 1:size(c_pfc.chanNames,1)
                    for freqi = 1:54 
                        for timei = 1:bins %parfor is optimal here 
                            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                            t_vvs = squeeze(c_vvs.oneListTraces(chani,timeBins,triali));
                            t_pfc = squeeze(c_pfc.oneListTraces(chanj,timeBins,triali));
        
                            EEG_vvs.data    = t_vvs;
                            EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                            EEG_vvs         = pop_eegfiltnew (EEG_vvs, freqi,freqi);
                            data_vvs        = squeeze(EEG_vvs.data); 
                            dataHA_vvs      = angle(hilbert(data_vvs));
        
                            EEG_pfc.data    = t_pfc;
                            EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                            EEG_pfc         = pop_eegfiltnew (EEG_pfc, freqi,freqi);
                            data_pfc        = squeeze(EEG_pfc.data); 
                            dataHA_pfc      = angle(hilbert(data_pfc));
        
                            diffPha = dataHA_vvs - dataHA_pfc;
                                    
                            myPLV(count, chani, chanj, freqi, timei ) = abs(mean(exp(1i*(diffPha))));
                        end
                    end
                end
            end
        end
    end

    PLVch{subji,:} = myPLV; 
end


cd (currentF)
%save('PLVch', 'PLVch');



toc



















%% 


clear
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_1-8Hz_2500-5500

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 1); 
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];
x = cellfun(@size, pfc, 'un', 0)

clear allR
for subji = 1:length(pfc)
   
    plvT = squeeze(mean(mean(PLVch{subji},3),2));
    plvT(plvT==0) = [];
    tmp = squeeze(pfc{subji}(7,:,30:52,13:25)); %last layer in beta and time period of interest
    tmp = squeeze(mean(mean(tmp,3),2));
    ch2(subji,:) = mean(tmp);
    %tmp = reshape (tmp, size(tmp, 1), []);
    
    
    allR(subji,:) = corr(plvT, tmp, 'type','s');
         
    
end

[h p ci ts] = ttest(allR);
disp(['t >>  ' num2str(ts.tstat)]);



%% 





















%% build 2 compare each frequency
clearvars 
currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\vvs\allTrials';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\pfc\allTrials';

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

freqs2use = [3:29 30:5:150];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    count = 1;
    for triali = 1:size(c_vvs.oneListTraces,3)
        id2u = strsplit(c_vvs.oneListIds_c{triali});
        if strcmp(id2u{1},'7') & ~strcmp(id2u{2},'4') 
            for freqi = 1:length(freqs2use)
                f = freqs2use(freqi);
                for chani = 1:size(c_vvs.chanNames,1)
                    for chanj = 1:size(c_pfc.chanNames,1)
                       t_vvs = squeeze(c_vvs.oneListTraces(chani,2500:5500,triali));
                       t_pfc = squeeze(c_pfc.oneListTraces(chanj,2500:5500,triali));

                       EEG_vvs.data    = t_vvs;
                       EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                       EEG_vvs         = pop_eegfiltnew (EEG_vvs, f,f);
                       data_vvs        = squeeze(EEG_vvs.data); 
                       phase_data(1,:)      = angle(hilbert(data_vvs));

                       EEG_pfc.data    = t_pfc;
                       EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                       EEG_pfc         = pop_eegfiltnew (EEG_pfc, f,f);
                       data_pfc        = squeeze(EEG_pfc.data); 
                       phase_data(2,:) = angle(hilbert(data_pfc));

                       PLVch{subji}(freqi, count, chani, chanj ) = abs(mean(exp(1i*(diff(phase_data,1)))));
                    end
                end
            end
            count = count+1;
        end
    end

    
end


cd (currentF)
save('PLVch_allF_2500-5500', 'PLVch');



toc



%% plot for each freq
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_allF_2500-5500

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 1); 
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];


clear allR
for subji = 1:length(pfc)
   
    for freqi = 1:52
        conT = squeeze(mean(mean(PLVch{subji}(freqi,:,:,:),4),3))';
        tmp = squeeze(pfc{(subji)}(1,:,11:27,7:12)); %last layer in beta and time period of interest
        tmp = squeeze(mean(mean(tmp,3),2));
        %ch2(subji,:) = mean(tmp);
        %tmp = reshape (tmp, size(tmp, 1), []);
        allR(subji,freqi) = corr(conT, tmp, 'type','s');
    end 
    
end

[h p ci ts] = ttest(allR)
 
plot(h)


%% plot for each freq 2
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_allF_2500-5500
%load PLVch_allF

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 2); 
s2e_pfc(4) = 1;
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];



clear allR
for subji = 1:length(pfc)
   for freqi = 1:52
    for freqi2 = 1:52
        conT = squeeze(mean(mean(PLVch{subji}(freqi,:,:,:),4),3))';
        %tmp = squeeze(pfc{(subji)}(1,:,freqi2:freqi2,7:12)); %last layer in beta and time period of interest
        
        %nornmalized to baseline
        tmp = squeeze(pfc{(subji)}(1,:,freqi2:freqi2,7:12)); %last layer in beta and time period of interest
        %tmpb = mean(squeeze(pfc{(subji)}(1,:,freqi2:freqi2,1:5)),2); %last layer in beta and time period of interest
        %tmp = tmp-tmpb;
        
        tmp = squeeze(mean(mean(tmp,3),2));
        %ch2(subji,:) = mean(tmp);
        %tmp = reshape (tmp, size(tmp, 1), []);
        allR(subji,freqi, freqi2) = corr(conT, tmp, 'type','s');
    end 
   end
end

[h p ci ts] = ttest(allR)
 

figure()
contourf(myresizem(squeeze(ts.tstat), 10), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(squeeze(h), 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);



%% 
d2p = squeeze(PLVch{1}(22,:,:));
figure();
imagesc(d2p); colorbar



%% load trialIDs

d2p = squeeze(all_r_Times(16,56,:,:));
imagesc(d2p)



%% all_r-times rho values show expected increases


d2p = squeeze(mean(all_r_Times(pfc_ids,56,:,:),'omitnan'));
figure()
imagesc(d2p); colorbar


d2p2 = all_r_Times_Trials(pfc_ids)









































%%