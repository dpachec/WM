%% 
% % TRIAL LEVEL ANAYLSIS FOR ALL FREQUENCY BANDS

clear , clc
paths = load_paths_WM('pfc');
main_path_pfc = paths.currSessEE; 
paths = load_paths_WM('vvs');
main_path_vvs = paths.currSessEE;  


fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

for bandi = 1:length(fnames)

    flist = dir([main_path_vvs fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_vvs fnames{bandi} '/' f2load])
    ALLE_VVS_B{bandi,:}  = ALL_EE([7 9 13 18 19 20 21 23 27 28]); 
    ids_VVS_B{bandi,:}  = all_IDs([7 9 13 18 19 20 21 23 27 28]); 
    
    flist = dir([main_path_pfc fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_pfc fnames{bandi} '/' f2load])
    ALLE_PFC_B{bandi,:}  = ALL_EE([2 3  5  9 10 11 12 14 15 16]);
    ids_PFC_B{bandi,:}  = all_IDs([2 3  5  9 10 11 12 14 15 16]);
    


end

%% analysis for all frequency bands (precomputed diagonal)

clearvars -except ALLE_PFC_B ALLE_VVS_B ids_VVS_B ids_PFC_B

timeL = 8:15; 

tic

for bandi = 1:length(ALLE_VVS_B)

    bandi 

    ids_VVS = ids_VVS_B{bandi};
    ids_PFC = ids_PFC_B{bandi}; 
    ALLE_VVS = ALLE_VVS_B{bandi}; 
    ALLE_PFC = ALLE_PFC_B{bandi}; 

    for subji = 1:10
    
        ids = ids_VVS{subji}; 
        its2cmp = cellfun(@(x) x(9:11), ids(:,2), 'un', 0);
        cats2cmp = cellfun(@(x) x(9), ids(:,2), 'un', 0);
        
        clear corrTR_VVS corrTR_PFC
        for triali = 1:length(ids)
            idx = ids(triali,1); 
            it2cmp = idx{1}(9:11);
            cat2cmp = idx{1}(9);

            %ids2testSameCat = ids(strcmp(ids(:, 1), idx) & ~strcmp(it2cmp, its2cmp) & strcmp(cat2cmp, cats2cmp),:)
            %ids2testDiffCat = ids(strcmp(ids(:, 1), idx) & ~strcmp(cat2cmp, cats2cmp),:)
            cc1VVS = ALLE_VVS{subji}(strcmp(ids(:, 1), idx) & ~strcmp(it2cmp, its2cmp) & strcmp(cat2cmp, cats2cmp), timeL); 
            cc2VVS = ALLE_VVS{subji}(strcmp(ids(:, 1), idx) & ~strcmp(cat2cmp, cats2cmp), timeL); 
            cc1PFC = ALLE_PFC{subji}(strcmp(ids(:, 1), idx) & ~strcmp(it2cmp, its2cmp) & strcmp(cat2cmp, cats2cmp), timeL); 
            cc2PFC = ALLE_PFC{subji}(strcmp(ids(:, 1), idx) & ~strcmp(cat2cmp, cats2cmp), timeL); 
            
    
            corrTR_VVS(triali,:) = mean(cc1VVS, 'all', 'omitnan') - mean(cc2VVS, 'all', 'omitnan') ;
            corrTR_PFC(triali,:) = mean(cc1PFC, 'all', 'omitnan') - mean(cc2PFC, 'all', 'omitnan') ;
            
        end
    
    
        
        rhoALL(subji, :) = corr(corrTR_PFC, corrTR_VVS, 'type', 's');
    
    end


    disp('done')
    
    [h(bandi) p(bandi) ci ts] = ttest(rhoALL);
    t(bandi) = ts.tstat ;
    rhoALL_B(:, bandi) = rhoALL; 



end


toc


save ('rhoALL_B', 'rhoALL_B', 't', 'p', 'h')


%% LOAD all conditions
%% ENCODNG
clear, clc

region = 'vvs'; 
paths = load_paths_WM(region);

%contrasts = { 'DISC_EE' 'DIDC_EE'};
contrasts = { 'DISC_EE' };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEE c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    idData{i,:} = all_IDs;
end

noAv = 1;
[out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
for i = 1:length(out_c) 
    eval([c{i} ' = out_c{i};']);
end

DISC_EE_VVS = DISC_EE; 
%DIDC_EE_VVS = DIDC_EE; 

clearvars -except DISC_EE_VVS DIDC_EE_VVS

region = 'pfc'; 
paths = load_paths_WM(region);

%contrasts = { 'DISC_EE' 'DIDC_EE'};
contrasts = { 'DISC_EE' };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEE c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    idData{i,:} = all_IDs;
end

noAv = 1;
[out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
for i = 1:length(out_c) 
    eval([c{i} ' = out_c{i};']);
end

DISC_EE_PFC = DISC_EE; 
%DIDC_EE_PFC = DIDC_EE; 



clearvars -except DISC_EE_VVS DIDC_EE_VVS DISC_EE_PFC DIDC_EE_PFC


%%

DISC_EE_PFC = DISC_EE_PFC([2 3  5  9 10 11 12 14 15 16]);
%DIDC_EE_PFC = DIDC_EE_PFC([2 3  5  9 10 11 12 14 15 16]);
DISC_EE_VVS = DISC_EE_VVS([7 9 13 18 19 20 21 23 27 28]); 
%DIDC_EE_VVS = DIDC_EE_VVS([7 9 13 18 19 20 21 23 27 28]); 



%%
clc
clear dVVS dPFC allR

tiledlayout(5,2); set(gcf, 'Position', [10 10 1000 1000])
for subji = 1:10
    ax= nexttile
    
    vvsC1 = DISC_EE_VVS{subji};
    %vvsC2 = DIDC_EE_VVS{subji};
    pfcC1 = DISC_EE_PFC{subji};
    %pfcC2 = DIDC_EE_PFC{subji};

    for triali = 1:size(vvsC1, 1)
        
        %take only diag
         diagVVS = diag(squeeze(vvsC1(triali, :, :)));
         dVVS(triali,:) = mean(diagVVS(6:15)); 
         diagPFC = diag(squeeze(pfcC1(triali, :, :)));
         dPFC(triali,:) = mean(diagPFC(6:15)); 

        % % % take full period
        %dVVS(triali,:) = mean(vvsC1(triali, 6:15, 6:15), 'all', 'omitnan'); 
        %dPFC(triali,:) = mean(pfcC1(triali, 6:15, 6:15), 'all', 'omitnan'); 

    end

    scatter(dVVS, dPFC,'.')
    h2 = lsline(ax);h2.LineWidth = 2;h2.Color = [.5 .5 .5 ]
    B = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
    allSlopes(subji, :) = B(2);
    allIntercepts(subji, :) = B(1);
    
    allR(subji, :) = corr(dVVS, dPFC, 'type', 's'); 
end


[h p ci ts] = ttest(allR);
t = ts.tstat; 
disp(['p = ' num2str(p) '   t = ' num2str(t)])

[h p ci ts] = ttest(allSlopes);
t = ts.tstat; 
disp(['p = ' num2str(p) '   t = ' num2str(t)])

%% plot one bar
clear data
ylim = [-.1 .4];


data.data = allR;
figure(); set (gcf, 'Position', [300 300 520 650]);
mean_S = mean(data.data, 1);std_S = std(data.data, [], 1); h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on; 
set(hb, 'Marker', '.', 'MarkerSize',60);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);set(hb,'linestyle','none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{''}, 'FontSize', 30, 'linew',2, 'ylim', ylim, 'xlim', [0 2]);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
ylabel('Rho')
exportgraphics(gcf, 'myP.png', 'Resolution', 100)





%% ENCODING MAINTENANCE 
clear, clc

region = 'vvs'; 
paths = load_paths_WM(region);

contrasts = { 'DISC_EM2' 'DIDC_EM2'};

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEM2 c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    idData{i,:} = all_IDs;
end

noAv = 1;
[out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
for i = 1:length(out_c) 
    eval([c{i} ' = out_c{i};']);
end

DISC_EM2_VVS = DISC_EM2; 
DIDC_EM2_VVS = DIDC_EM2; 

clearvars -except DISC_EM2_VVS DIDC_EM2_VVS

region = 'pfc'; 
paths = load_paths_WM(region);

contrasts = { 'DISC_EM2' 'DIDC_EM2'};

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEM2 c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    idData{i,:} = all_IDs;
end

noAv = 1;
[out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
for i = 1:length(out_c) 
    eval([c{i} ' = out_c{i};']);
end

DISC_EM2_PFC = DISC_EM2; 
DIDC_EM2_PFC = DIDC_EM2; 



clearvars -except DISC_EM2_VVS DIDC_EM2_VVS DISC_EM2_PFC DIDC_EM2_PFC


%%

DISC_EM2_PFC = DISC_EM2_PFC([2 3  5  9 10 11 12 14 15 16]);
DIDC_EM2_PFC = DIDC_EM2_PFC([2 3  5  9 10 11 12 14 15 16]);
DISC_EM2_VVS = DISC_EM2_VVS([7 9 13 18 19 20 21 23 27 28]); 
DIDC_EM2_VVS = DIDC_EM2_VVS([7 9 13 18 19 20 21 23 27 28]); 



%%
clc
clear dVVS dPFC allR

for subji = 1:10
    
    vvsC1 = DISC_EM2_VVS{subji};
    vvsC2 = DIDC_EM2_VVS{subji};
    pfcC1 = DISC_EM2_PFC{subji};
    pfcC2 = DIDC_EM2_PFC{subji};

    for triali = 1:size(vvsC1, 1)
        
        %take only diag
        dVVS(triali,:) = mean(vvsC1(triali, 6:15, 6:40), 'all'); 
        dPFC(triali,:) = mean(pfcC1(triali, 6:15, 6:40), 'all'); 

    end
    
    allR(subji, :) = corr(dVVS, dPFC, 'type', 's'); 
end

[h p ci ts] = ttest(allR)
t = ts.tstat; 
disp(['p = ' num2str(p) '   t = ' num2str(t)])




%% plot one bar
clear data
ylim = [-.1 .4];


data.data = allR;
figure(); set (gcf, 'Position', [300 300 520 650]);
mean_S = mean(data.data, 1);std_S = std(data.data, [], 1); h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on; 
set(hb, 'Marker', '.', 'MarkerSize',60);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);set(hb,'linestyle','none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{''}, 'FontSize', 30, 'linew',2, 'ylim', ylim, 'xlim', [0 2]);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
ylabel('Rho')
%exportgraphics(gcf, 'myP.png', 'Resolution', 300)




%% correlate representational geometries in both regions for all frequencies
clear, clc
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
paths_PFC = load_paths_WM('pfc');
paths_VVS = load_paths_WM('vvs');
filelistSessPFC = getFiles(paths_PFC.traces);
filelistSessVVS = getFiles(paths_VVS.traces);

filelistSessPFC = filelistSessPFC([2 3  5  9 10 11 12 14 15 16]);
filelistSessVVS = filelistSessVVS([7 9 13 18 19 20 21 23 27 28]); 

cfg.period      = 'M123';
cfg.timeRes     = 0.1; 
cfg.freqs       = [3:54]; 
cfg.win_width   = 5; 
cfg.mf          = 1; 
cfg.fR          = 1; 
cfg.avTFV       = 0; 

for sessi= 1:length(filelistSessVVS) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths_PFC.traces filelistSessPFC{sessi}]);   
    [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    %cfg_contrasts               = average_repetitions(cfg_contrasts);
    neuralRDMsPFC               = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
    
    load([paths_VVS.traces filelistSessVVS{sessi}]);   
    [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    %cfg_contrasts               = average_repetitions(cfg_contrasts);
    neuralRDMsVVS                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);

    nFreqs = length(cfg.freqs); 
    nTimes = size(neuralRDMsVVS, 4); 
    all_r_Times = zeros(nFreqs, nTimes);
    parfor freqi = 1:nFreqs
        for timei = 1:nTimes
            rdmVVS = squeeze(neuralRDMsVVS(:, :, freqi, timei));
            rdmVVS = vectorizeRDM(rdmVVS);
            rdmPFC = squeeze(neuralRDMsPFC(:, :, freqi, timei));
            rdmPFC = vectorizeRDM(rdmPFC);
            allTEst = corr(rdmVVS', rdmPFC', 'type', 's');
            all_r_Times(freqi,timei) = allTEst;  
        end
    end
    
    allRSESS{sessi} = all_r_Times; 
end

save('allRSeSS', 'allRSESS')


%% 


allRS = cat(3, allRSESS{:}); allRS = permute(allRS, [3, 1, 2]);

[h p ci ts] = ttest(allRS); 
h = squeeze(h); t = squeeze(ts.tstat); 

figure; 
%times = 1:210; 
times = 1:460; 
freqs = 1:520
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);





%% correlate representational geometries in both regions for all frequencies TRIAL LEVEL
clear, clc
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
paths_PFC = load_paths_WM('pfc');
paths_VVS = load_paths_WM('vvs');
filelistSessPFC = getFiles(paths_PFC.traces);
filelistSessVVS = getFiles(paths_VVS.traces);

filelistSessPFC = filelistSessPFC([2 3  5  9 10 11 12 14 15 16]);
filelistSessVVS = filelistSessVVS([7 9 13 18 19 20 21 23 27 28]); 

cfg.period      = 'M123';
cfg.timeRes     = 0.1; 
cfg.freqs       = [3:54]; 
cfg.win_width   = 5; 
cfg.mf          = 1; 
cfg.fR          = 1; 
cfg.avTFV       = 0; 

for sessi= 1:length(filelistSessVVS) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths_PFC.traces filelistSessPFC{sessi}]);   
    [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    %cfg_contrasts               = average_repetitions(cfg_contrasts);
    neuralRDMsPFC               = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
    
    load([paths_VVS.traces filelistSessVVS{sessi}]);   
    [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    %cfg_contrasts               = average_repetitions(cfg_contrasts);
    neuralRDMsVVS                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);

    nFreqs = length(cfg.freqs); 
    nTimes = size(neuralRDMsVVS, 4); 
    all_r_Times = zeros(nFreqs, nTimes);
    nTrials = size(neuralRDMsVVS, 1);
    all_r_Times = zeros(nTrials, nFreqs, nTimes);


    for triali = 1:size(neuralRDMsVVS, 1)
        for freqi = 1:nFreqs
            parfor timei = 1:nTimes
                rdmVVS = squeeze(neuralRDMsVVS(:, triali, freqi, timei));
                rdmVVS(triali) = []; 
                rdmPFC = squeeze(neuralRDMsPFC(:, triali, freqi, timei));
                rdmPFC(triali) = []; 
                allTEst = corr(rdmVVS, rdmPFC, 'type', 's');
                all_r_Times(triali, freqi,timei) = allTEst;  
            end
        end
    end
        
    allRSESS{sessi} = all_r_Times; 
end

save('allRSeSS', 'allRSESS')




%% correlate results with PFC fit
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%cd D:\_WM\analysis\pattern_similarity\across_regions
%load allRSeSS_M123_noAvRepet_TFV_TRIALS

clearvars -except allRSESS
f2sav = 'RNN_pfc_M123_[56]_3-54_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);
nnFit_PFC = nnFit([2 3  5  9 10 11 12 14 15 16]);; 

for subji = 1:length(nnFit_PFC)
    %nexttile
    trlsPFC = nnFit_PFC{subji};
    %mPFC = squeeze(mean(mean(trlsPFC(1, :, 15:29, 8:13), 3), 4))'; 


    for freqi = 1:52
        for timei = 1:46
            allRSH = allRSESS{subji};
            tBinTr = allRSH(:, freqi, timei);
            psBinTr = squeeze(trlsPFC(1, :, freqi, timei))'; 
            
            %rhoAll(subji, freqi, timei) = corr(tBinTr, mPFC, 'type', 's');
            rhoAll(subji, freqi, timei) = corr(tBinTr, psBinTr, 'type', 's');
        end
    end
    
    
    
end

%[h p ci ts] = ttest(rhoAll); 

%disp(['p = ' num2str(p) '   t = ' num2str(ts.tstat)])
%exportgraphics(gcf, [paths.results.DNNs 'myCorr.png'], 'Resolution', 300); 


%%

[h p ci ts] = ttest(rhoAll); 
h = squeeze(h); t = squeeze(ts.tstat);

figure; imagesc(h)



%% plot mean just 2 check 

d2p = squeeze(cell2mat(cellfun(@mean, nnFit_PFC, 'un', 0)))
md2p = squeeze(mean(d2p)); 

imagesc(md2p)












%%