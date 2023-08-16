%% Calculate from epoched raw traces
%% 
clear
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf

%f2sav = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORrtRELU_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BLNETeBatchNorm_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_E123_[1-8]_3-54_1_0_0_0_.1_5_1';
%
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-8_1_0_0_0_.1_5_1';
%f2sav = 'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'Res18-2_pfc_MALL_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'CORrt_vvs_M123_[8]_3-54_1_0_1_0_.1_5_1';
f2sav = 'BLNETe_pfc_M123NC_[8-8-56]_3-54_1_0_1_0_.1_5_1';



cfg = getParams(f2sav);
cfg.DNN_analysis = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
filelistSess = getFilesWM(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   

    [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 

    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);

    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    
    neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
    networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);

    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav '.mat'], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% IN LOOP 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {
                    'CAT_vvs_E123_[1]_3-54_1_0_1_0_.1_5_1'
                    'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1'
                    'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1'
             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.traces filelistSess{sessi}]);   
       
        [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if (cfg.avRep)
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
        
    end
    
    save([paths.results.DNNs f2sav '.mat'], 'nnFit');

end




%% PERMUTATIONS IN LOOP
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

clear
nPerm = 1000;

listF2sav = {

                    'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
                    'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
                    'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
                    
             };
    

for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi nPerm 
        
    f2sav       = listF2sav{listi}
        
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    t1 = datetime; 
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   

        nChans = size(cfg_contrasts.chanNames, 1); 
    
        if nChans > 1

            cfg_contrasts = getIdsWM(cfg.period, cfg_contrasts);
            cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
            
            cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
            if (cfg.avRep)
                cfg_contrasts               = average_repetitions(cfg_contrasts);
            end
        
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
    
    
            % % % restrict time for permutation data
            if strcmp(cfg.period(1), 'M')
                if ndims(neuralRDMs) == 4
                    neuralRDMs = neuralRDMs(:,:,:,6:40); %frequency-resolved
                else
                    neuralRDMs = neuralRDMs(:,:,6:40); %band analysis
                end
            else
                if ndims(neuralRDMs) == 4
                    neuralRDMs = neuralRDMs(:,:,:,6:15); %frequency-resolved
                else
                    neuralRDMs = neuralRDMs(:,:,6:15); %band analysis
                end
            end
            
            for permi = 1:nPerm
                sC = size(networkRDMs, 2);
                ids = randperm(sC);
                networkRDMs = networkRDMs(:, ids, ids); 
                nnFitPerm(permi, sessi,:,:, :)              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
            end

        end

    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav '_' num2str(nPerm) 'p.mat'], 'nnFitPerm', '-v7.3');
    

end
   
t2 = datetime; 
etime(datevec(t2), datevec(t1))


%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
%f2sav = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'BLNETe_vvs_E123_[56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
f2sav = 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1';

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,8, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:39));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 39], 'clim', [-5 5], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 

%% plot final figure only last layer / time point ENCODING
times = 1:150;
freqs = 1:520; 
h = zeros(52, 15); 

h(clustinfo.PixelIdxList{1}) = 1; % Alexnet VVS Layer 1

%h(clustinfo.PixelIdxList{5}) = 1; % BLNETi PFC
%h(clustinfo.PixelIdxList{1}) = 1; % BLNETi VVS

%h(clustinfo.PixelIdxList{1}) = 1; % BLNETe cluster 1
%h(clustinfo.PixelIdxList{2}) = 1; % BLNETe cluster 2

%h(clustinfo.PixelIdxList{2}) = 1; %Cornet
%h(clustinfo.PixelIdxList{4}) = 1; %Cornet cluster 2


figure; set(gcf, 'Position', [100 100 200 400])
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);

set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 150], 'clim', [-5 5], 'FontSize', 10);
set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 

%% plot final figure only last layer MAINTENANCE
times = 1:390;
freqs = 1:520; 
%h = zeros(52, 39); 

%h(clustinfo.PixelIdxList{2}) = 1; %BLNETi PFC

%h(clustinfo.PixelIdxList{3}) = 1; %CORNET PFC
%h(clustinfo.PixelIdxList{15}) = 1; %CORNET VVS

%h(clustinfo.PixelIdxList{5}) = 1; %category model

%h(clustinfo.PixelIdxList{2}) = 1; %pfc - Cornet
%h(clustinfo.PixelIdxList{23}) = 1; %pfc - Cornet

%h(clustinfo.PixelIdxList{8}) = 1; %pfc - ecoset
%h(clustinfo.PixelIdxList{10}) = 1; %vvs1 - ecoset
%h(clustinfo.PixelIdxList{23}) = 1; %vvs2 - ecoset


figure; set(gcf, 'Position', [1000 918 560 420])
myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)
%contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
%contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);

set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 390], 'clim', [-4 4], 'FontSize', 10);
%set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% COMPUTE CLUSTERS in each permutation FREQUENCY RESOLVED
clc
clearvars -except allTObs f2sav
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

% use the same name as the previously plotted file
f2sav = [f2sav '_1000p.mat'];

cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

cd ([paths.results.DNNs])
load(f2sav);

if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end

f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{13}(1:end-5))));
nLays = size(nnFitPerm, 3);

for permi = 1:nPerm
    
    for layi = 1:nLays
        
        dataP = atanh(squeeze(nnFitPerm(permi, :,layi, :,:)));
        %dataP = squeeze(nnFitPerm(permi, :,layi, :,:));
        dataP(sub2exc, :, :) = []; 
        [h p ci ts] = ttest(dataP);
        h = squeeze(h); t = squeeze(ts.tstat);
        h(:, 1:5) = 0; % only sum t-values in cluster after the baseline
        
        
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        %sort all clusters 
        if exist('allSTs')
            [max2u id] = max(allSTs);
        else
            allSTs = 0; 
            id = 1;
        end
        
        max_clust_sum_perm(permi,layi,:) = allSTs(id); 
    end

end

cd (paths.github)


%% compute p value bands for all layers FREQ RES
clc
clear p

for layi = 1:size(allTObs, 1)
    clear mcsR mcsP
    for pixi = 1:size(allTObs, 2)
        mcsR = allTObs(layi, pixi); 

        if mcsR ~= 0
            mcsP = squeeze(max_clust_sum_perm(:,layi));
        
            %allAb = mcsP(abs(mcsP) > abs(mcsR));
            allAb = mcsP(mcsP > mcsR);
            p(layi, pixi,:) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
        else
            p(layi, pixi,:) = nan; 
        end
    end
end

p (p==1.0010 | p == 1) = nan; 
p_ranked = p; p_ranked(isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));
p

%% p last layer only
p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))



%%
x = p(4,:); 
x = sort(x)


%% load file to plot BANDS (ALL LAYERS RNN and Alex) THIS IS NOT THE FINAL FIGURE PLOTS (CODE BELOW)
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear, clc 

f2sav = 'BLNETi_pfc_E123_[8-8-56]_9-12_1_0_0_0_.1_5_1'; 


cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,8, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    
    clear nnH
    for subji = 1:length(nnFit)
        if strcmp (cfg.period(1), 'M')
            nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:39));
        else
            nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:15));
        end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    h = squeeze(h); t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end
    %tObs = sum(t(clustinfo.PixelIdxList{2}))
    d2p = squeeze(mean(nnH, 'omitnan'));
    
    %times = -1.75:.1:6.849; 
    
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    
    mART = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(mART, 2);
        shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
        scatter(times, hb, 'Linewidth', 4)
        set(gca, 'xlim', [1 37]); 
        set(gca, 'FontSize', 12, 'ylim', [-.03 .035]);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth',1);
        plot(get(gca,'xlim'), [0 0],'k:','lineWidth',1);

    elseif strcmp(cfg.period(1), 'E')
        times = 1:size(mART, 2);
        shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
        scatter(times, hb, 200,'.','Linewidth', 4)
        set(gca, 'FontSize', 12, 'ylim', [-.015 .1]);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth',1);
        plot(get(gca,'xlim'), [0 0],'k:','lineWidth',1);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 



%% load file to plot BANDS (ALL LAYERS RNN and Alex) -- IN one PLOT ONLY 
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear %, clc
%f2sav = 'BLNETi_vvs_E123_[8-8-56]_39-54_1_0_0_0_.1_5_1'; 
f2sav = 'Alex_vvs_M123_[1-8]_9-12_1_0_0_0_.1_5_1';
%f2sav = 'CORrt_pfc_M123_[2-2-8]_39-54_1_0_0_0_.1_5_1';
%f2sav = 'CAT_pfc_E123_[1]_39-54_1_0_0_0_.1_5_1';


cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);

clear hbL
for layi = 1:size(nnFit{1}, 1)
    clear nnH
    for subji = 1:length(nnFit)
        if strcmp (cfg.period(1), 'M')
            nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:39));
        else
            nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:15));
        end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    %h(1:5) = 0; 
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(layi, :, :) = 0;
    end
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    hbL(layi, :) = hb; 
    
    mART(layi, :) = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    

end


if strcmp(cfg.period(1), 'E')

     figure()
    set(gcf, 'Position', [100 100 150 300])
    times = 1:15;
    plot(times, mART, 'Linewidth', 3); hold on; 
    
    % ALEXNET
    if strcmp(cfg.net2load, 'Alex')%ALEXNET
        myCmap = colormap(brewermap(20,'RdPu'));
        myCmap = myCmap([6 8 10 12 14 16 18 20],:);
        
        if strcmp(cfg.brainROI, 'pfc') 
            x = (-.010:-.0022:-.026)';
        elseif strcmp(cfg.brainROI, 'vvs') 
            x = (-.010:-.0035:-.035)';
        end
        
        %hbL([1:8], :) = nan; % ALL
       

    elseif strcmp(cfg.net2load, 'BLNETi')%BLNET% BLNET
         myCmap = colormap(brewermap(18,'RdPu'));
         myCmap = myCmap([6 8 10 12 14 16 18],:);
         if strcmp(cfg.brainROI, 'pfc') 
            x = (-.012:-.002:-.025)'; % if PFC 
        elseif strcmp(cfg.brainROI, 'vvs') 
            x = (-.012:-.0034:-.035)'; % if VVS
        end
     
      if strcmp(cfg.brainROI, 'vvs') & cfg.freqs(1) == 39 %VVS HIGH GAMMA
            
            hbL([4], 1:5) = nan; % VVS HIG GAMMA
            %hbL([5:6], :) = nan; % VVS HIG GAMMA
            %hbL(7, 12:15) = nan; % VVS HIG GAMMA
      end

     %hbL([1:7],:) = nan; % ALL 
    
    %hbL([1:2 5:6 ], :) = nan; % PFC THETA
    %hbL([7], 13:15) = nan; % PFC THETA
    %hbL([1:7], :) = nan; % PFC GAMNMA = ALL 
    

      elseif strcmp(cfg.net2load, 'CORrt')%BLNET%CorNET
            x = (-.015:-.0035:-.028)';
            myCmap = colormap(brewermap(18,'RdPu'));
            myCmap = myCmap([6 10 14 18],:);
            if strcmp(cfg.brainROI, 'pfc') & cfg.freqs(1) == 9 %VVS THETA
                hbL(1:4, :) = nan; % ALL
            end
            if strcmp(cfg.brainROI, 'pfc') & cfg.freqs(1) == 9 %VVS ALPHA
                hbL(2, :) = nan; % ALL
            end
            %hbL(1, 10:15) = nan; % VVS ALPHA
            %hbL([1:2 4], 10:15) = nan; % VVS HIGH GAMMA
    
        elseif strcmp(cfg.net2load, 'CAT') % CATEGORY MODEL
            x = -.025;
            myCmap = colormap(brewermap(18,'RdPu'));
            myCmap = myCmap([6 10 14 18],:);
            %hbL(1, :) = nan; % ALL

    end
    
    
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 15]); 
    if strcmp(cfg.brainROI, 'vvs')
        set(gca, 'FontSize', 12, 'ylim', [-.0375 .1]);
    elseif strcmp(cfg.brainROI, 'pfc')
        set(gca, 'FontSize', 12, 'ylim', [-.0275 .0375]);
    end
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);
    
    %legend

    
% % % % % % MAINTENANCE
elseif strcmp(cfg.period(1), 'M')
    
    set(gcf, 'Position', [100 100 405 163])
    times = 1:39;
    plot(times, mART, 'Linewidth', 3); hold on; 

    
    if strcmp(cfg.net2load, 'Alex')%ALEXNET

        myCmap = colormap(brewermap(20,'RdPu'));
        myCmap = myCmap([6 8 10 12 14 16 18 20],:);

        x = (-.012:-.0035:-.038)';

        
        if strcmp(cfg.brainROI, 'vvs') & cfg.freqs(1) == 3 %VVS THETA
            hbL(8, 37:39) = nan; 
        elseif strcmp(cfg.brainROI, 'vvs') & cfg.freqs(1) == 9 %VVS ALPHA
            %hbL([1 4], :) = nan; 
            hbL(8, 1:20) = nan; 
        elseif strcmp(cfg.brainROI, 'pfc') & cfg.freqs(1) == 13 %VVS ALPHA
            %hbL([1 4], :) = nan; 
            hbL([3], 25:32) = nan;  
            hbL([1:2 5:7], :) = nan;  
            
        else
            %hbL(1:8, :) = nan; % ALL
        end
        
    elseif strcmp(cfg.net2load, 'BLNETi')%BLNET

        myCmap = colormap(brewermap(18,'RdPu'));
        myCmap = myCmap([6 8 10 12 14 16 18],:);
        x = (-.015:-.0035:-.038)';
        %
        %hbL([1:4 7], :) = nan; % VVS THETA
        %hbL(6, 1:20) = nan; % VVS THETA
        %hbL([1:4], 1:20) = nan; % VVS ALPHA
        %hbL([7], 21:end) = nan; % VVS ALPHA
        if strcmp(cfg.brainROI, 'pfc') & cfg.freqs(1) == 13 % PFC BETA
            hbL(1:6, :) = nan; 
        else %% all 
            %hbL(1:7, :) = nan; % ALL
        end

    elseif strcmp(cfg.net2load, 'CORrt')%BLNET

        x = (-.015:-.0035:-.028)';
        myCmap = colormap(brewermap(18,'RdPu'));
        myCmap = myCmap([6 10 14 18],:);
        if strcmp(cfg.brainROI, 'vvs') & cfg.freqs(1) == 9 % VVS THETA
            hbL(1:2, 1:25) = nan; 
            hbL(3, 1:10) = nan; 
        elseif strcmp(cfg.brainROI, 'pfc') & cfg.freqs(1) == 13 % VVS THETA %% all 
            hbL(1:3, 1:39) = nan; % ALL
        else 
            hbL(1:4, 1:39) = nan; % ALL
        end
        %hbL(1, 1:37) = nan; % VVS ALPHA
        %hbL(1:4, 1:25) = nan; % VVS ALPHA
        %hbL(1:4, 1:25) = nan; % PFC BETA

    end

    
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 39]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);
   

end


set(gca, 'ColorOrder', myCmap)



exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 



%% compute clusters in each permutation BANDS for all layers

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%f2sav =  'RNN_pfc_M123_[8-8-56]_13-29_1_1_0_0_.1_5_1_1000p'; 
%f2sav = 'CORrt_vvs_E123_[1-8]_39-54_1_0_0_0_.1_5_1_1000p';
%f2sav = 'Alex_vvs_E123_[1-8]_3-8_1_0_0_0_.1_5_1_1000p';
%f2sav = 'BLNETi_pfc_M123_[8-8-56]_13-29_1_0_0_0_.1_5_1_1000p'; 
%f2sav = 'CAT_vvs_M123_[1]_3-8_1_0_0_0_.1_5_1_1000p';
%f2sav = 'CORrt_vvs_M123_[2-2-8]_39-54_1_0_0_0_.1_5_1_1000p';

f2sav = [f2sav '_1000p'];

cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
cd ([paths.results.DNNs])
load([f2sav '.mat']);

if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end

f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{13}(1:end-1))));
if ndims (nnFitPerm) == 3
    nnFitPerm1(:, :, 1,:) = nnFitPerm;
    nnFitPerm = nnFitPerm1;

end


for layi = 1:size(nnFitPerm, 3)
    for permi = 1:nPerm
        
        dataP = atanh(squeeze(nnFitPerm(permi, :,layi, :)));
        dataP(sub2exc, :,:) = []; 
        [h p ci ts] = ttest(dataP);
        h = squeeze(h); t = squeeze(ts.tstat);
    
        
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        %sort all clusters 
        if exist('allSTs')
            [max2u id] = max(allSTs);
            max_clust_sum_perm(permi,layi, :) = allSTs(id); 
        else
            max_clust_sum_perm(permi,layi, :) = 0; 
        end
    
    
    end
end


cd (paths.github)


%% compute p value bands for all layers
clc
clear p

for layi = 1:size(allTObs, 1)
    clear mcsR mcsP
    for pixi = 1:size(allTObs, 2)
        mcsR = allTObs(layi, pixi); 

        if mcsR ~= 0
            mcsP = squeeze(max_clust_sum_perm(:,layi));
        
            %allAb = mcsP(abs(mcsP) > abs(mcsR));
            allAb = mcsP(mcsP > mcsR)
            p(layi, pixi,:) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
        else
            p(layi, pixi,:) = nan; 
        end
    end
end

p (p==1.001 | p == 1) = 0; 
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));

p


%%
%figure
%histogram(mcsP); hold on; 
%catter(mcsR,0, 'filled','r');
    

%% compute p for one layer only
clear p mcsR mcsP

mcsR =    65.4844891668492
mcsP = squeeze(max_clust_sum_perm(:, 7));

for clusti = 1:length(mcsR)
    
    %allAb = mcsP(abs(mcsP) > abs(mcsR));
    allAb = mcsP(mcsP > mcsR);
    p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

end


%% load files to plot all CATEGORY RESULTS  -- IN one PLOT ONLY 

clear, clc 

listF2sav = {       
%     
                    'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1';
                    'CAT_pfc_E123_[1]_9-12_1_0_0_0_.1_5_1';
                    'CAT_pfc_E123_[1]_13-29_1_0_0_0_.1_5_1';
                    'CAT_pfc_E123_[1]_30-38_1_0_0_0_.1_5_1';
                    'CAT_pfc_E123_[1]_39-54_1_0_0_0_.1_5_1';


%                     'CAT_pfc_M123_[1]_3-8_1_0_0_0_.1_5_1';
%                     'CAT_pfc_M123_[1]_9-12_1_0_0_0_.1_5_1';
%                     'CAT_pfc_M123_[1]_13-29_1_0_0_0_.1_5_1';
%                     'CAT_pfc_M123_[1]_30-38_1_0_0_0_.1_5_1';
%                     'CAT_pfc_M123_[1]_39-54_1_0_0_0_.1_5_1';


% %                     'CAT_vvs_E123_[1]_3-8_1_0_0_0_.1_5_1';
% %                     'CAT_vvs_E123_[1]_9-12_1_0_0_0_.1_5_1';
% %                     'CAT_vvs_E123_[1]_13-29_1_0_0_0_.1_5_1';
% %                     'CAT_vvs_E123_[1]_30-38_1_0_0_0_.1_5_1';
% %                     'CAT_vvs_E123_[1]_39-54_1_0_0_0_.1_5_1';

%                     'CAT_vvs_M123_[1]_3-8_1_0_0_0_.1_5_1';
%                     'CAT_vvs_M123_[1]_9-12_1_0_0_0_.1_5_1';
%                     'CAT_vvs_M123_[1]_13-29_1_0_0_0_.1_5_1';
%                     'CAT_vvs_M123_[1]_30-38_1_0_0_0_.1_5_1';
%                     'CAT_vvs_M123_[1]_39-54_1_0_0_0_.1_5_1';


             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi allFnnH
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.results.DNNs);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.results.DNNs listF2sav{listi} '.mat']);   
        
        nnH = cat(1, nnFit{:,1});
        
    end

    allFnnH(listi, :, :) = atanh(nnH);


end

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end
allFnnH(:, sub2exc, :) = []; 
%allFnnH = squeeze(mean(allFnnH, 2));


%% 
for freqi = 1:size(allFnnH, 1)
    if strcmp(cfg.period(1), 'M')
        nnH = allFnnH(freqi, :, 1:39);
    else
        nnH = allFnnH(freqi, :, 1:15);
    end
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    %h(1:5) = 0; 
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(freqi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(freqi, :, :) = 0;
    end
    
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    clear hbl
    hbL(freqi, :) = hb; 
    
    mART(freqi, :) = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    

end


if strcmp(cfg.period(1), 'M')
    %set(gcf, 'Position', [100 100 405 163])
    set(gcf, 'Position', [100 100 400 260])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    times = 1:39;
    x = (-.015:-.003:-.028)';
    hbL([1:5], :) = nan;     
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 39]); 
    set(gca, 'FontSize', 12, 'ylim', [-.04 .04]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    %colormap(jet(5));
    times = 1:15;
    plot(times, mART, 'Linewidth', 5); hold on; 

    if strcmp(cfg.brainROI, 'pfc')
        x = (-.018:-.0025:-.030)';%PFC
        hbL([2], :) = nan; %PFC
    else
        x = (-.012:-.0050:-.032)';%VVS
        hbL(:, 1:4) = nan; %VVS
    end

    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 5); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 15]); 
    if strcmp(cfg.brainROI, 'pfc')
        set(gca, 'FontSize', 12, 'ylim', [-.03 .04]);
    elseif strcmp(cfg.brainROI, 'vvs')
        set(gca, 'FontSize', 12, 'ylim', [-.04 .15]);
    end
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);
    %legend

end


set(gca, 'ColorOrder', myCmap)



exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 








%% MULTI - ITEM TRIALS

clear
%...__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials;1:Trials)__timeRes__win-width__mf_FST
f2sav = 'BLnext2_vvs_MALL_[6]_3-54_0_0_1_0_.1_5_1'; 
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
filelistSess = getFiles(paths.traces);

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    nChans = size(cfg_contrasts.chanNames, 1); 
    
    if nChans > 1
        [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if cfg.avRep
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts.oneListIds, sessi, paths);
        
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    end
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav '.mat'], 'nnFit');







%%  plot all layers MULTI-ITEM
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials;1:Trials)__timeRes__win__mf_FST
clear 
f2sav = 'BLnext2_pfc_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 

cfg = getParams(f2sav);


if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end
 

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(4,12)
if strcmp(cfg.period(1), 'E')
    set(gcf, 'Position', [100 100 1400 800])
elseif strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 2000 500])
end

for layi = 1:size(nnFit{2}, 1)-4
    nexttile
    clear nnH
    for subji = 1:length(nnFit)
       if ~isempty(nnFit{subji, 1})
            nnH(subji, : ,:) = nnFit{subji, 1}(layi,:,:);       
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    
    endT = size(nnFit{2}, 3);
    times = 1:endT; 

    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'E')
        set(gca, 'xlim', [1 15], 'clim', [-5 5], 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}, 'FontSize', 8);
        set(gca, 'xtick', [1 6.5 15], 'xticklabels', {'-.5' '0' '1.5'}, 'FontSize', 8);
        plot([6.5 6.5],get(gca,'ylim'), 'k:','lineWidth', 2);
    elseif strcmp(cfg.period(1), 'M')
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}, 'FontSize', 8);
        set(gca, 'xtick', [1 6.5 40], 'xticklabels', {'-.5' '0' '3.5'}, 'FontSize', 8);
        plot([6.5 6.5],get(gca,'ylim'), 'k:','lineWidth', 2);
        if strcmp(cfg.period(1:2), 'M1')
            set(gca, 'xlim', [12 40], 'clim', [-5 5], 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}, 'FontSize', 8);
            set(gca, 'xtick', [12 40], 'xticklabels', {'0' '1.5'}, 'FontSize', 8);
            plot([6.5 6.5],get(gca,'ylim'), 'k:','lineWidth', 2);
        end
    end
    
    
    
end

 exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 



%% MULTI - ITEM in loop 

clear

listF2sav = {   
                'BLnext12_pfc_E41_[7]_3-54_0_0_1_0_.1_5_1_AV.mat' 
                'BLnext12_pfc_E42_[7]_3-54_0_0_1_0_.1_5_1_AV.mat' 
                'BLnext12_pfc_E43_[7]_3-54_0_0_1_0_.1_5_1_AV.mat' 
             };
    

for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi 


    f2sav = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        
        load([paths.traces filelistSess{sessi}]);   
        [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
        
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts.oneListTraces, cfg.timeRes); % 
        cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if (cfg.avRep)
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
        [networkRDMs ids2rem]       = createNetworkRDMs(cfg, cfg_contrasts.oneListIds, sessi, paths);
        
        neuralRDMs(ids2rem, :, : ,:) = []; 
        neuralRDMs(:,ids2rem, : ,:) = []; 
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
        
    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav], 'nnFit');
    

end


