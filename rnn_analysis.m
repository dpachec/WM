%% Calculate from epoched raw traces
%% first load traces
clear
%ROI__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
%example f2sav = 'RNN_pfc_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
%f2sav = 'RNN_vvs_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
f2sav = 'RNN_vvs_M_[56]_30-38_1_0_0_0_.1_5_1.mat'; 
cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        %this version is incorrect (not considering all trials but only
        %single -item... )
        ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
    
    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% IN LOOP 
clear
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf


listF2sav = {
                'Alex_vvs_E_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'Alex_pfc_E_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'Alex_hipp_E_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'Alex_vvs_M_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'Alex_pfc_M_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'Alex_hipp_M_[1-8]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_vvs_E_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_pfc_E_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_hipp_E_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_vvs_M_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_pfc_M_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
                'RNN_hipp_M_[1-56]_3-54_1_1_1_0_.1_5_1.mat'; 
             };
    

for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}
    cfg = getParams(f2sav);
    f2t = strsplit(f2sav, '_');
    region = f2t{2};
    paths = load_paths_WM(region);
    filelistSess = getFiles(paths.traces);
    
    t1 = datetime; 
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   
       
        if strcmp(cfg.period, 'M')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(cfg.period, 'E')
            ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
        end
        oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
        cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
        cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if (cfg.avRep)
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
        networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
        
    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav], 'nnFit');
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%% load file to plot frequency resolved
%ROI__layers__freqs__avRepet__avTFV__fRes(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
clear 
f2sav = 'Alex_vvs_E_[1-8]_3-54_1_0_1_0_.1_5_1.mat'
sub2exc = [];
lay2plot = 7; 

f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);



clear nnH
for subji = 1:length(nnFit)
    
   nnH(subji, : ,:) = nnFit{subji, 1}(lay2plot,:,:);
   %nnH(subji, : ,:) = nnFit{subji, 1}(7,:);
   %nnH(subji, : ,:,:) = squeeze(nnFit{subji, 1});
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);



%%plot frequency resolved
d2p = squeeze(mean(nnH, 'omitnan'));
freqs = 1:520; 
times = -1.75:.01:6.849; 
clustinfo = bwconncomp(h);
%clen = max(cellfun(@length, clustinfo.PixelIdxList))


% for pixi = 1:length(clustinfo.PixelIdxList)
%    h(clustinfo.PixelIdxList{pixi}) = 0;   
% end
% 
% % %h(clustinfo.PixelIdxList{25}) = 1; % PFC
% % %tObs = sum(t(clustinfo.PixelIdxList{25})) % PFC
% % 
% % h(clustinfo.PixelIdxList{18}) = 1; % VVS
% % h(clustinfo.PixelIdxList{10}) = 1; % VVS
% % %h(clustinfo.PixelIdxList{34}) = 1; % VVS
% % tObs = sum(t(clustinfo.PixelIdxList{18})) % VVS
% 
% h(clustinfo.PixelIdxList{29}) = 1; % PFC
% tObs = sum(t(clustinfo.PixelIdxList{29})) % PFC


%  h(clustinfo.PixelIdxList{7}) = 1; % VVS ENV
%  h(clustinfo.PixelIdxList{11}) = 1; % VVS ENV
%  tObs = sum(t(clustinfo.PixelIdxList{7})) % VVS ENC
 
%  
% h(clustinfo.PixelIdxList{15}) = 1; % HIPP ENV
% tObs = sum(t(clustinfo.PixelIdxList{15})) % HIPP ENC


myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
set(gca, 'xlim', [-1 4], 'FontSize', 40, 'clim', [-5 5])
%set(gca, 'xlim', [-.5 1.5], 'FontSize', 40, 'clim', [-5 5])
plot([0 0],get(gca,'ylim'), 'k:','lineWidth', 6);
%set(gca, 'clim', [-.025 .025])

%clustinfo = bwconncomp(h);


%exportgraphics(gcf, 'myP.png', 'Resolution', 300); 

%%  plot all layers 
clear 
f2sav = 'Alex_vvs_E_[1-8]_3-54_1_0_1_0_.1_5_1.mat'
sub2exc = [];
lay2plot = 7; 

f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);

tiledlayout(3,3)
set(gcf, 'Position', [100 100 1200 1000])
for layi = 1:size(nnFit{1}, 1)
    nexttile
    clear nnH
    for subji = 1:length(nnFit)
       nnH(subji, : ,:) = nnFit{subji, 1}(layi,:,:);       
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    times = -1.75:.1:6.849; 
    clustinfo = bwconncomp(h);
    
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 4);
    if strcmp(f2t{3}, 'E')
        set(gca, 'xlim', [-.5 1], 'clim', [-5 5], 'xtick', [], 'xticklabels', [])
    else
        set(gca, 'xlim', [-.5 4], 'clim', [-5 5], 'xtick', [], 'xticklabels', [])
    end
    plot([0 0],get(gca,'ylim'), 'k:','lineWidth', 6);
    
    %exportgraphics(gcf, 'myP.png', 'Resolution', 300); 

end
%% load file to plot BANDS
%ROI__layers__freqs__avRepet__avTFV__fRes(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
clear 
f2sav = 'RNN_vvs_E_[56]_30-38_1_0_0_0_.1_5_1.mat'; 
f2t = strsplit(f2sav, '_');
region = f2t{2};

paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);

if strcmp(region , 'pfc')
    sub2exc = [1];
elseif strcmp(region , 'vvs')
    sub2exc = [18 22];
elseif strcmp(region , 'hipp')
    sub2exc = [2];
end

for subji = 1:length(nnFit)
    
   %nnH(subji, : ,:) = nnFit{subji, 1}(1,:,:);
   %nnH(subji, : ,:) = nnFit{subji, 1}(7,:);
   nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
tObs = sum(t(clustinfo.PixelIdxList{2}))
d2p = squeeze(mean(nnH, 'omitnan'));

times = -1.75:.1:6.849; 
figure
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 

mART = squeeze(mean(nnH)); 
stdART = squeeze(std(nnH)); 
seART = stdART/ sqrt(size(nnH, 1));

shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
plot (times, hb, 'Linewidth', 4)
%set(gca, 'FontSize', 40, 'xlim', [-.5 1])
set(gca, 'FontSize', 40, 'xlim', [-.5 4])
%exportgraphics(gcf, 'myP.png', 'Resolution', 300); 




%% 

figure()
for layi = 1:56
    d2p = squeeze(mean(nnH(:, layi,:), 'omitnan'));
    plot(d2p); hold on; 
end
%h(h==0) = nan; h(h==1) = .02;
%plot(h, 'lineWidth', 2)



%% PERMUTATIONS

clear
nPerm = 100;
%ROI__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
%example f2sav = 'pfc_8-16-24-32-40-48-56_13-29_0_1_500_1_1'; 
%f2sav = ['RNN_pfc_M_56_3-54_1_0_1_0_.1_5_1_p' num2str(nPerm) '.mat']; 

f2sav = ['RNN_hipp_E_[56]_3-54_1_0_1_0_.1_5_1_p' num2str(nPerm) '.mat']; 

cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
    
    if ndims(neuralRDMs) == 4
        disp ('time resolved')
        neuralRDMs = neuralRDMs(:,:,:,21:30);
    else
        disp ('no time resolved')
        neuralRDMs = neuralRDMs(:,:,21:30);
    end
    
    for permi = 1:nPerm
        sC = size(networkRDMs, 2);
        ids = randperm(sC);
        networkRDMs = networkRDMs(:, ids, ids); 
        nnFitPerm(permi, sessi,:,:, :)              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    end
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFitPerm');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% PERMUTATIONS IN LOOP

clear
nPerm = 100;

listF2sav = {
                'Alex_vvs_E_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'Alex_pfc_E_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'Alex_hipp_E_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'Alex_vvs_M_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'Alex_pfc_M_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'Alex_hipp_M_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_vvs_E_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_pfc_E_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_hipp_E_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_vvs_M_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_pfc_M_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                'RNN_hipp_M_[1-56]_3-54_1_0_1_0_.1_5_1.mat'; 
                
             };
    

for listi = 2:length(listF2sav)
    
    clearvars -except listF2sav listi nPerm 
        
    f2sav       = listF2sav{listi}
        
    cfg = getParams(f2sav);
    f2t = strsplit(f2sav, '_');
    region = f2t{2};
    paths = load_paths_WM(region);
    filelistSess = getFiles(paths.traces);
    
    t1 = datetime; 
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   
       
        if strcmp(cfg.period, 'M')
            ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
        elseif strcmp(cfg.period, 'E')
            ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
        end
        oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
        cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
        cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if (cfg.avRep)
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
        networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
        
        if ndims(neuralRDMs) == 4
            disp ('time resolved')
            neuralRDMs = neuralRDMs(:,:,:,21:30);
        else
            disp ('no time resolved')
            neuralRDMs = neuralRDMs(:,:,21:30);
        end
        
        for permi = 1:nPerm
            sC = size(networkRDMs, 2);
            ids = randperm(sC);
            networkRDMs = networkRDMs(:, ids, ids); 
            nnFitPerm(permi, sessi,:,:, :)              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        end
        
    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav(1:end-4) '_' num2str(nPerm) 'p.mat'], 'nnFitPerm');
    

end
   
t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% compute clusters in each permutation frequency resolved

nPerm = 200; 

sub2exc = [2]; 

for permi = 1:nPerm
    
    dataP = squeeze(nnFitPerm(permi, :,1, :,:));
    dataP(sub2exc, :, :) = []; 
    [h p ci ts] = ttest(dataP);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 

end

%% compute clusters in each permutation BANDS
sub2exc = [2]; 
nPerm = 200; 
for permi = 1:nPerm
    
    dataP = squeeze(nnFitPerm(permi, :,:)); %maintenance
    %dataP = squeeze(nnFitPerm(permi, :,1:10)); %encoding
    
    dataP(sub2exc, :) = []; 
    [h p ci ts] = ttest(dataP);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 

end

%%
obsT = tObs; 

figure
histogram(max_clust_sum_perm); hold on; 
scatter(obsT, 0, 175,'filled','r');

mcsR = obsT; 
mcsP = max_clust_sum_perm; 
allAb = mcsP(abs(mcsP) > abs(mcsR));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm




%% MULTI - ITEM TRIALS
clear
%...__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
f2sav = 'BLnext4_vvs_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '5') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    
    

    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    [networkRDMs ids2rem]       = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, ...
                                    sessi, paths, cfg.period);
    neuralRDMs(ids2rem, :, : ,:) = []; 
    neuralRDMs(:,ids2rem, : ,:) = []; 
    
    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))








%% load file to plot multi-item
%ROI__layers__freqs__avRepet__avTFV__fRes(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
clear 
f2sav = 'BLnext12_hipp_M_[6:7]_3-54_0_0_1_0_.1_5_1.mat'; 
sub2exc = [];

f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);

clear nnH
for subji = 1:length(nnFit)
    
   nnH(subji, : ,:) = nnFit{subji, 1}(2,:,:);
   %nnH(subji, : ,:) = nnFit{subji, 1}(7,:);
   %nnH(subji, : ,:,:) = squeeze(nnFit{subji, 1});
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);



%%plot frequency resolved
d2p = squeeze(mean(nnH, 'omitnan'));
freqs = 1:520; 
times = -1.75:.01:6.849; 
clustinfo = bwconncomp(h);
%clen = max(cellfun(@length, clustinfo.PixelIdxList))


% for pixi = 1:length(clustinfo.PixelIdxList)
%    h(clustinfo.PixelIdxList{pixi}) = 0;   
% end
%  h(clustinfo.PixelIdxList{11}) = 1; % VVS ENV
%  tObs = sum(t(clustinfo.PixelIdxList{7})) % VVS ENC
 
 

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
set(gca, 'xlim', [-1 4], 'FontSize', 40, 'clim', [-5 5])
%set(gca, 'xlim', [-.5 1.5], 'FontSize', 40, 'clim', [-5 5])
plot([0 0],get(gca,'ylim'), 'k:','lineWidth', 6);
%set(gca, 'clim', [-.025 .025])


%exportgraphics(gcf, 'myP.png', 'Resolution', 300); 
%% MULTI - ITEM in loop 

clear


listF2sav = {
                'BLnext2_hipp_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_hipp_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_hipp_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_hipp_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext2_vvs_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_vvs_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_vvs_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_vvs_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext2_pfc_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_pfc_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_pfc_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_pfc_E_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext2_hipp_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_hipp_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_hipp_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_hipp_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext2_vvs_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_vvs_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_vvs_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_vvs_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext2_pfc_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext4_pfc_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext8_pfc_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 
                'BLnext12_pfc_M_[6:7]_3-54_0_0_1_0_.1_5_1_AV.mat'; 




             };
    

for listi = 1:length(listF2sav)
    
clearvars -except listF2sav listi 

f2sav       = listF2sav{listi}
cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '5') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    
    

    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    [networkRDMs ids2rem]       = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, ...
                                    sessi, paths, cfg.period);
    neuralRDMs(ids2rem, :, : ,:) = []; 
    neuralRDMs(:,ids2rem, : ,:) = []; 
    
    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))




end





%% load file to plot multi-item
%ROI__layers__freqs__avRepet__avTFV__fRes(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
clear 
f2sav = 'BLnext8_vvs_M_[6:7]_3-54_0_0_1_0_.1_5_1.mat'
sub2exc = [];

f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);

clear nnH
for subji = 1:length(nnFit)
    
   nnH(subji, : ,:) = nnFit{subji, 1}(2,:,:);
   %nnH(subji, : ,:) = nnFit{subji, 1}(7,:);
   %nnH(subji, : ,:,:) = squeeze(nnFit{subji, 1});
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);



%%plot frequency resolved
d2p = squeeze(mean(nnH, 'omitnan'));
freqs = 1:520; 
times = -1.75:.01:6.849; 
clustinfo = bwconncomp(h);
%clen = max(cellfun(@length, clustinfo.PixelIdxList))


% for pixi = 1:length(clustinfo.PixelIdxList)
%    h(clustinfo.PixelIdxList{pixi}) = 0;   
% end
%  h(clustinfo.PixelIdxList{11}) = 1; % VVS ENV
%  tObs = sum(t(clustinfo.PixelIdxList{7})) % VVS ENC
 
 
figure()
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
set(gca, 'xlim', [-1 4], 'FontSize', 40, 'clim', [-5 5])
%set(gca, 'xlim', [-.5 1.5], 'FontSize', 40, 'clim', [-5 5])
plot([0 0],get(gca,'ylim'), 'k:','lineWidth', 6);
%set(gca, 'clim', [-.025 .025])


%exportgraphics(gcf, 'myP.png', 'Resolution', 300); 

%% PLOT OBS DATA 

sub2exc = [];

for subji = 1:length(nnFit)
   
   nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);

sub2exc =[];

all_r_Times = nnH; 

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52;

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   + 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          %hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end



figure()
layT = tiledlayout(3, 3);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    elseif size(all_r_Times, 4) == 546
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:5.009;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:52;
        contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        %plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        %plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        %set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = -1.75:.01:6.849; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'clim', [-3 3]); % 'xlim', [0 3.5],
        set(gca, 'xlim', [-1 6])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

%f2p = [num2str(layi, '%02.f') '.png'];
%exportgraphics(gcf, f2p, 'Resolution', 300)







%% STIMULUS REPRESENTATION IN DNNs
% % % 
clear
f2sav       = 'RNN';
loadNet_WM;


%% STIMULUS REPRESENTATION IN DNNs
% % % 
clear
f2sav       = 'RNN';
[act_CH act_FR] = load_DNNs;
















%% RNN all RDMS

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   title(num2str(layi))
end

%exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% RNN all RDMS

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   %title(num2str(layi)) % check order
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+coli/9 0+rowi/9 n n ])
        count = count+-1;
        
    end
end

%exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+rowi/9 0+coli/9 n n ])
        count = count-1;
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer / timepoint


clearvars -except act_CH act_FR
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(act_CH, 1)
    d2p = squeeze(act_CH(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


figure()
plot (CCI, 'Linewidth', 2)
set(gca, 'FontSize', 20)

figure
lw = 4;
plot (CCI(1:8), 'Linewidth', lw, 'Color', cols(1, :)); hold on; 
plot (CCI(9:16), 'Linewidth', lw, 'Color', cols(2, :)); hold on; 
plot (CCI(17:24), 'Linewidth', lw, 'Color', cols(3, :)); hold on; 
plot (CCI(25:32), 'Linewidth', lw, 'Color', cols(4, :)); hold on; 
plot (CCI(33:40), 'Linewidth', lw, 'Color', cols(5, :)); hold on; 
plot (CCI(41:48), 'Linewidth', lw, 'Color', cols(6, :)); hold on; 
plot (CCI(49:56), 'Linewidth', lw, 'Color', cols(7, :)); hold on; 
set(gca, 'FontSize', 30, 'xlim', [0 9], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);




%%  permutations
clearvars -except act_CH act_FR

nPerm = 1000; 

for permi = 1:nPerm
%create cateogy model
idSh = randperm(60*60);
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 
M = M(idSh); M = reshape(M, [60 60]);

    for layi = 1:size(act_CH, 1)
        d2p = squeeze(act_FR(layi, :,:)); 
        d2p = 1- d2p;
        rdmMDS = d2p; 
        rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
        mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
        mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
        %CCIP(permi, layi) = (mAcross - mWithin) / (mAcross + mWithin);
        CCIP(permi, layi) = (mAcross - mWithin);

    end
    
end

%% 
figure()
histogram(CCIP(:, 49))



%% plot only last time point

figure(); set(gcf, 'Position', [100 100 500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (3, 3, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25
 set(ha(1),'position',[0 0 n n ])
 set(ha(2),'position',[.55 .27 n n ])
 set(ha(3),'position',[.275 .27 n n])
 set(ha(4),'position',[.0 .27 n n ])
 set(ha(5),'position',[.55 .54 n n ])
 set(ha(6),'position',[.275 .54 n n])
 set(ha(7),'position',[.0 .54 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot only last time point one line

figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[ 6/9 0 n n ])
 set(ha(2),'position',[ 5/9 0 n n ])
 set(ha(3),'position',[ 4/9 0 n n])
 set(ha(4),'position',[ 3/9 0 n n ])
 set(ha(5),'position',[ 2/9 0 n n ])
 set(ha(6),'position',[ 1/9 0 n n])
 set(ha(7),'position',[ 0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% plot only last time point one line vertical

figure(); set(gcf, 'Position', [100 100 700 700]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .1; 
 set(ha(1),'position',[0 6/9 n n ])
 set(ha(2),'position',[0 5/9 n n ])
 set(ha(3),'position',[0 4/9 n n])
 set(ha(4),'position',[0 3/9 n n ])
 set(ha(5),'position',[0 2/9 n n ])
 set(ha(6),'position',[0 1/9 n n])
 set(ha(7),'position',[0 .0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% plot only last time point one line ALEXNET

figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% plot MDS only last time point one line
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 56];
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *6  0.1 n n ])
 set(ha(2),'position',[.09 * 5 0.1 n n ])
 set(ha(3),'position',[.09 * 4 0.1 n n])
 set(ha(4),'position',[.09 * 3 0.1 n n ])
 set(ha(5),'position',[.09 * 2 0.1 n n ])
 set(ha(6),'position',[.09 0.1 n n])
 set(ha(7),'position',[.0 0.1 n n ])
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot MDS only last time point one line ALEXNET
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% MDS all layers (black and yellow circles plot)

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
c1 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c2 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c3 = repmat(zeros(1), 56, 1);

cols = [c1 c2 c3];



d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);
figure()
plot(rdmMDS(1:8,1),rdmMDS(1:8,2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS(9:16,1),rdmMDS(9:16,2),'k', 'linewidth', 2)
plot(rdmMDS(17:24,1),rdmMDS(17:24,2),'k', 'linewidth', 2)
plot(rdmMDS(25:32,1),rdmMDS(25:32,2),'k', 'linewidth', 2)
plot(rdmMDS(33:40,1),rdmMDS(33:40,2),'k', 'linewidth', 2)
plot(rdmMDS(41:48,1),rdmMDS(41:48,2),'k', 'linewidth', 2)
plot(rdmMDS(49:56,1),rdmMDS(49:56,2),'k', 'linewidth', 2)
scatter(rdmMDS(:,1),rdmMDS(:,2),3500,cols, '.'); 
set(gca,'FontSize', 26);



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

%exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% MDS all layers (pink and blue triangle-circles plot) only for first and last time points

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');%allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
%c2 = repmat(zeros(1), 7, 1)';

c2 = repmat(zeros(1), 7, 1)';

c3 = repmat(ones(1), 7, 1)';

cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);


allM = allM([1 8 9 16 17 24 25 32 33 40 41 48 49 56], [1 8 9 16 17 24 25 32 33 40 41 48 49 56]); 

d2p = 1- allM;

[rdmMDS] = cmdscale(d2p);
figure()
plot(rdmMDS([1 2],1),rdmMDS([1 2],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([3 4],1),rdmMDS([3 4],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([5 6],1),rdmMDS([5 6],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([7 8],1),rdmMDS([7 8],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([9 10],1),rdmMDS([9 10],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([11 12],1),rdmMDS([11 12],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([13 14],1),rdmMDS([13 14],2),'k', 'linewidth', 2);hold on; 



fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 7);

for i = 1:14
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gca,'FontSize', 30, 'xlim', [-.4 .4], 'ylim', [-.2 .4]); 
end



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% Representational consistency all layers / time points

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 15)


%% Representational consistency only first and last time point

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 


act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 20)
set(gca, 'xtick', [1:7], 'xticklabel', {[1:7]}, 'ytick', [1:7], 'yticklabel', {[1:7]},'clim', [.5 1])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% MDS Representational consistency only first and last time point

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 


act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);

[rdmMDS] = cmdscale(allMS);

clear c1 c2 c3 cols
c1 = (1:7)/7';
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';

figure()
for i = 1:7
   scatter(rdmMDS(i,1),rdmMDS(i,2),5350, cols(i, :), '.'); hold on; axis square
   
end
set(gca,'FontSize', 40, 'xlim', [-.5 .5], 'ylim', [-.5 .5]); 

exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% Representational consistency ALEXNET

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 20)
set(gca, 'xtick', [1:8], 'xticklabel', {[1:8]}, 'ytick', [1:8], 'yticklabel', {[1:8]},'clim', [0.3 1])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Changes in representational consistency 

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  

c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


clear allM2
for tlyi = 1:55

    allM2(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');allMS = allM;
    
end

figure()
lw = 3;

plot(abs(allM2(1:7)), 'Linewidth', lw, 'Color', [cols(1,:)]); hold on; 
plot(abs(allM2(9:15)), 'Linewidth', lw, 'Color', [cols(2,:)])
plot(abs(allM2(17:23)), 'Linewidth', lw, 'Color', [cols(3,:)])
plot(abs(allM2(25:31)), 'Linewidth', lw, 'Color', [cols(4,:)])
plot(abs(allM2(33:39)), 'Linewidth', lw, 'Color', [cols(5,:)])
plot(abs(allM2(41:47)), 'Linewidth', lw, 'Color', [cols(6,:)])
plot(abs(allM2(49:55)), 'Linewidth', lw, 'Color', [cols(7,:)])
%legend({'Layer 1' 'Layer 2' 'Layer 3' 'Layer 4' 'Layer 5' 'Layer 6' 'Layer 7' })
set(gca, 'ylim', [0.8 1], 'xlim', [0 8], 'xtick', [1:7], 'xticklabels', {'1-2' '2-3' '3-4' '4-5' '5-6' '6-7' '7-8'}, 'FontSize', 30)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% example pytorch DQN

layers = importONNXLayers('test_model.onnx','OutputLayerType', 'classification', 'ImportWeights', true)
layers = removeLayers(layers,'Transpose_0')
layers = removeLayers(layers,'Input_input')
inputlayer = imageInputLayer([84 84 3],'Name','input', 'Normalization', 'none')
layers = addLayers(layers,inputlayer)
layers = connectLayers(layers,'input', 'Conv_1')
analyzeNetwork(layers)
%%
net = assembleNetwork(layers);

%% 
analyzeNetwork(net)

%%
im = rand(84, 84, 3);

act1 = activations(net,im,'Conv_1');


%% START HERE WITH THE ALL TRIALS ANALYSIS
% 1) bands
clearvars -except act_CH act_FR 
f2sav       = 'Alex_pfc_CueAll_noAv_54_Real_nT_1-49'
f           = 3:54; %in case B
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    f2t = strsplit(f2sav, '_'); 
    if strcmp(f2t{5}, '54') | strcmp(f2t{5}, '150') 
        for freqi = 1:length(freqs2test)
            f  = freqs2test(freqi);
            rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
        end
    elseif strcmp(f2t{5}, 'B') 
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
        if subji == 1 % change name only once
            f2sav = [f2sav num2str(f(1)) '-' num2str(f(end)) 'Hz'];
        end
    end
    
    [allS ids act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR); 
   
    
    for layi = 1:nLays
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(layi,:,:)); 
        else
            M =  squeeze(act_CH2(layi,:,:)); 
        end
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;

        allTEst = corr(allS, M', 'type', 's');

        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
 
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');


%% PLOT OBS DATA (all layers) - case with mix negative and positive clusters

sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52; % from 3 to 54

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   %+ 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end



figure()
layT = tiledlayout(8, 7);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 56:56 %1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = 0:.01:3.49;
        times = 0:.01:4.49;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-3 3])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

f2p = 'myFig.png';
exportgraphics(gcf, f2p, 'Resolution', 300)





%% plot NICELY

h = zeros(52, 45);

figure(); set(gcf, 'Position', [100 100 1500 1000]);
% Scheme|'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'|
myCmap = colormap(brewermap([],'*Spectral'));
times = 0:.01:4.49;
%freqs = [.1:.1:29.9 30:.5:150]; 
freqs = 1:520;
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 10); %colorbar
plot([.45 .45], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot([.95 .95], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
%plot([1.25 1.25], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 12);
set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-4 4]); colorbar
set(gca,'XTick',[], 'YTick', [])
set(gca,'xticklabel',[], 'FontSize', 12)
    
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Frequency resolved RNN analysis for all trials and 8th time point only 

clearvars -except act_CH act_FR 

%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 56]; 
f2sav       = 'RNN_pfc_Maint_Av_54_Real_nT_loi_8th'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, length(lois), nFreqs, nTimes);

tic

for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
    end
    
    [allS ids act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR); 
   
    
    parfor layi = 1:length(lois)
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(lois(layi),:,:)); 
        else
            M =  squeeze(act_CH2(lois(layi),:,:)); 
        end
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;

        allTEst = corr(allS, M', 'type', 's');

        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
 
 
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');

%% Plot in 3 different time periods vertical separately for first and last time points in each layer

%allTC= squeeze(mean(all_r_Times(:, [1 9 17 25 33 41 49], :, :), 4, 'omitnan')); % first layer
%allTC= squeeze(mean(all_r_Times(:, [8 16 24 32 40 48 56], :, :), 4, 'omitnan')); % lastlayer
allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [100 100 100 500])
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])




exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 








%% Plot in 3 different time periods vertical

nLays = 56; 

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';


figure()
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal
set(gca, 'xlim', [0 560])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal







%% export sequences

clear
%ROI__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
%example f2sav = 'RNN_pfc_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
%f2sav = 'RNN_vvs_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
f2sav = 'RNN_pfc_M_[56]_30-38_1_0_0_0_.1_5_1.mat'; 
cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3

    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
    ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    ids2save = cfg_contrasts.oneListIds_c(ids); 
    id3 = cellfun(@(x) strsplit(x), ids2save, 'un', 0);
    id4 = cellfun(@(x) double(string(x(:, 13:15))), id3, 'un', 0);
    id5 = cat(1, id4{:});
    all_seqs{sessi,:} = id5; 
    filename = ['sub_' num2str(sessi,  '%02.f') '_image_seqs_ALL.csv']
    csvwrite(filename, id5);

end

%%

sublist = dir('*.csv'); sublist = {sublist.name}'; 
clear all_seqs1
for subji = 1:length(sublist)
    all_seqs1{subji,:} = csvread(sublist{subji}); 

end




 
