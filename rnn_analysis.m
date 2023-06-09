%% Calculate from epoched raw traces
%% first load traces
clear
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf

%f2sav = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'Res18-2_pfc_MALL_[3]_3-54_0_0_1_0_.1_5_1.mat'; 
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'CORrtRELU_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';
%f2sav = 'BLNETeBatchNorm_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';
f2sav = 'AlexEco_pfc_E123_[1-8]_3-54_1_0_0_0_.1_5_1.mat';

cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);
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

    
    neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
    networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts.oneListIds_c, sessi, paths);

    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% IN LOOP 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {
                    %'Res18-2_vvs_MALL_[4]_3-54_0_0_1_0_.1_5_1.mat';                 
                    %'Res18-2_pfc_MALL_[4]_3-54_0_0_1_0_.1_5_1.mat';                 
                    %'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1.mat'; 
                    %'CAT_vvs_M123_[1]_3-54_1_0_1_0_.1_5_1.mat'; 
                    %'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';
                    %'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';
                    %'CORrt_pfc_E123_[1-8]_30-38_1_0_0_0_.1_5_1.mat'; 
                    %'CORrt_pfc_E123_[1-8]_39-54_1_0_0_0_.1_5_1.mat'; 
                    
                    'CAT_vvs_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    'CAT_vvs_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';
                    

                    


             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.traces filelistSess{sessi}]);   
       
        [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
        tic
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        toc
        
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
        if (cfg.avRep)
            cfg_contrasts               = average_repetitions(cfg_contrasts);
        end
    
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts.oneListIds_c, sessi, paths);
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
        
    end
    
    save([paths.results.DNNs f2sav], 'nnFit');

end




%%  plot all layers ALEX frequency resolved
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
clear , clc

%f2sav = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat';
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1.mat';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end


paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);


tiledlayout(3,3)
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1000])
end
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
    clustinfo = bwconncomp(h);
    if strcmp(cfg.period(1), 'M')
        if size(t, 2) > 40
            times = 1:46; 
        else
            times = 1:40; 
        end
    else
        times = 1:15; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}, 'FontSize', 16);
        set(gca, 'xtick', [ 3 40], 'xticklabels', {'0' '3.5'}, 'FontSize', 24);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'clim', [-5 5], 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}); 
        set(gca, 'xtick', [1 5.5 15], 'xticklabels', {'-.5' '0' '1'}, 'FontSize', 22);
        plot([5.5 5.5],get(gca,'ylim'), 'k:','lineWidth', 2);
    end
    
    
    
end

 exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


 %% plot final figure only last layer MAINTENANCE ALEXNET
times = 1:400;
freqs = 1:520; 
h = zeros(52, 40); 


%h(clustinfo.PixelIdxList{2}) = 1; %pfc - Cornet
%h(clustinfo.PixelIdxList{23}) = 1; %pfc - Cornet

%h(clustinfo.PixelIdxList{8}) = 1; %pfc - ecoset
%h(clustinfo.PixelIdxList{10}) = 1; %vvs1 - ecoset
%h(clustinfo.PixelIdxList{23}) = 1; %vvs2 - ecoset

myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)
%contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
%contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);


if layi == 49
    set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
    set(gca, 'xtick', [3  37], 'xticklabels', {'0' '3.5'});
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
end
set(gca, 'xlim', [1 370], 'clim', [-5 5], 'FontSize', 10);
%set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 

%%  plot all layers RNN frequency resolved
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
%f2sav = 'RNN_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'BLnext2_pfc_MALL_[6]_3-54_0_0_1_0_.1_5_1.mat'; 
%f2sav = 'CORrt_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat';
%f2sav = 'Res18-2_pfc_MALL_[3]_3-54_0_0_1_0_.1_5_1.mat'; 
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav =  'BLNETi_vvs_M123_[1-56]_3-54_1_0_1_0_.1_5_1.mat';
%f2sav = 'CORrtRELU_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1.mat'; 
%f2sav = 'BLNETeBatchNorm_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';
%f2sav = 'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1.mat';
f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1.mat';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);


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
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:37));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:11));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    
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
        times = 1:11; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        if layi == 49
            set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
            set(gca, 'xtick', [3  37], 'xticklabels', {'0' '3.5'});
        else
            set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        end
        set(gca, 'xlim', [1 37], 'clim', [-5 5], 'FontSize', 10);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        if layi == 49
            set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}); 
            set(gca, 'xtick', [1 5.5 15], 'xticklabels', {'-.5' '0' '1'})
        else
            set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        end
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 

%% plot final figure only last layer / time point RNN ENCODING
times = 1:110;
freqs = 1:520; 
h = zeros(52, 11); 
%h(clustinfo.PixelIdxList{3}) = 1; % Ecoset 
%h(clustinfo.PixelIdxList{4}) = 1; % Ecoset cluster 2

%h(clustinfo.PixelIdxList{4}) = 1; %Cornet
%h(clustinfo.PixelIdxList{5}) = 1; %Cornet cluster 2


figure; set(gcf, 'Position', [100 100 200 400])
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);


if layi == 49
    set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
    set(gca, 'xtick', [3  37], 'xticklabels', {'0' '3.5'});
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
end
set(gca, 'xlim', [1 110], 'clim', [-5 5], 'FontSize', 10);
set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 

%% check frequencies 
x1(:,1) = 30:5:150
x1(:,2) = 30:1:54

%% plot final figure only last layer MAINTENANCE
times = 1:370;
freqs = 1:520; 
h = zeros(52, 37); 

h(clustinfo.PixelIdxList{9}) = 1; %BL-NET PFC

%h(clustinfo.PixelIdxList{7}) = 1; %BL-NET PFC

%h(clustinfo.PixelIdxList{5}) = 1; %category model

%h(clustinfo.PixelIdxList{2}) = 1; %pfc - Cornet
%h(clustinfo.PixelIdxList{23}) = 1; %pfc - Cornet

%h(clustinfo.PixelIdxList{8}) = 1; %pfc - ecoset
%h(clustinfo.PixelIdxList{10}) = 1; %vvs1 - ecoset
%h(clustinfo.PixelIdxList{23}) = 1; %vvs2 - ecoset

myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)
%contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
%contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);


if layi == 49
    set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
    set(gca, 'xtick', [3  37], 'xticklabels', {'0' '3.5'});
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
end
set(gca, 'xlim', [1 370], 'clim', [-4 4], 'FontSize', 10);
%set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% COMPUTE CLUSTERS in each permutation FREQUENCY RESOLVED

clearvars -except allTObs
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%f2sav = 'RNN_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_1000p.mat';
f2sav = 'CORrt_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1_1000p.mat';
%f2sav =  'RNN_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_1000p.mat'; 
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1_1000p.mat'; 
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_1000p.mat';
%f2sav =  'BLNETe_vvs_E123_[56]_3-54_1_0_1_0_.1_5_1_1000p.mat';
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1_1000p.mat';

cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);

%cd ([paths.results.DNNs 'permutations\full_period\'])
cd ([paths.results.DNNs ])
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
        dataP = squeeze(nnFitPerm(permi, :,layi, :,:));
        dataP(sub2exc, :, :) = []; 
        [h p ci ts] = ttest(dataP);
        h = squeeze(h); t = squeeze(ts.tstat);

        
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        %sort all clusters 
        %[max2u id] = max(abs(allSTs));
        [max2u id] = max(allSTs);
        
        max_clust_sum_perm(permi,layi,:) = allSTs(id); 
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
            allAb = mcsP(mcsP > mcsR);
            p(layi, pixi,:) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
        else
            p(layi, pixi,:) = nan; 
        end
    end
end

p (p==1.0010 | p == 1) = 0; 
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))

%% p last layer only
p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))




%% load file to plot BANDS (one Layer - time Point)
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
f2sav = 'CAT_pfc_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';


cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);

if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
    %sub2exc = [1 7 13 16]; %at least 3 electrodres
    %sub2exc = [1 3 4 7 10 13 16]; %at least 4 electrodres
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end

for subji = 1:length(nnFit)
    % 8    16    24    32    40    48    56
   nnH(subji, : ,:) = nnFit{subji, 1}(1,:);
   %nnH(subji, : ,:) = nnFit{subji, 1}(7,:);
   %nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
%tObs = sum(t(clustinfo.PixelIdxList{1}))
d2p = squeeze(mean(nnH, 'omitnan'));

%times = -1.75:.1:6.849; 

figure
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 

mART = squeeze(mean(nnH)); 
stdART = squeeze(std(nnH)); 
seART = stdART/ sqrt(size(nnH, 1));

if strcmp(cfg.period(1), 'M')
    times = 1:46;
    %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
    shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
    plot (times, hb, 'Linewidth', 4)
    %set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'xtick', [3 8 37], 'xticklabels', {'0' '0.5' '3.5'}, 'xlim', [1 37])
    set(gca, 'FontSize', 18, 'ylim', [-.03 .035]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',1);
    plot([8 8],get(gca,'ylim'), 'k:','lineWidth',1);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',1);
elseif strcmp(cfg.period(1), 'E')
    times = 1:21; 
    %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}); 
    shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
    plot (times, hb, 'Linewidth', 4)
    set(gca, 'xtick', [1 5.5 15], 'xticklabels', {'-.5' '0' '1'}, 'xlim', [1 15])
    set(gca, 'FontSize', 12);
    plot([5.5 5.5],get(gca,'ylim'), 'k:','lineWidth',2);
    
end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 





%% load file to plot BANDS (ALL LAYERS RNN and Alex)
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear, clc 
%f2sav = 'RNN_vvs_M123_[8-8-56]_9-12_1_0_0_0_.1_5_1.mat'; 
%f2sav = 'RNNe_pfc_E123_[8-8-56]_3-8_1_0_0_0_.1_5_1.mat'; 
f2sav = 'CORrt_pfc_E123_[1-8]_30-38_1_0_0_0_.1_5_1.mat';




cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end

paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);


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
       %nnH(subji, : ,:) = nnFit{subji, 1}(layi,1:37);
       nnH(subji, : ,:) = nnFit{subji, 1}(layi,:);       
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
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
     
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    
    mART = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    
    if strcmp(cfg.period(1), 'M')
        times = 1:37;
        %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});
        shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
        %plot (times, hb, 'Linewidth', 4)
        scatter(times, hb, 'Linewidth', 4)
        set(gca, 'xlim', [1 37]); 
        %set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
        set(gca, 'FontSize', 12, 'ylim', [-.03 .035]);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth',1);
        plot(get(gca,'xlim'), [0 0],'k:','lineWidth',1);

    elseif strcmp(cfg.period(1), 'E')
        times = 1:21;
        %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}); 
        shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
        %plot (times, hb, 'Linewidth', 4)
        scatter(times, hb, 200,'.','Linewidth', 4)
        %set(gca, 'xtick', [1 5.5 15], 'xticklabels', {'-.5' '0' '1'}, 'xlim', [1 15])
        %set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 12]); 
        set(gca, 'FontSize', 12, 'ylim', [-.015 .1]);
        plot([3 3],get(gca,'ylim'), 'k:','lineWidth',1);
        plot(get(gca,'xlim'), [0 0],'k:','lineWidth',1);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 



%% load file to plot BANDS (ALL LAYERS RNN and Alex) -- IN one PLOT ONLY 
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear, clc 
%f2sav = 'RNNe_pfc_E123_[8-8-56]_3-8_1_0_0_0_.1_5_1.mat'; 
f2sav = 'CORrt_pfc_M123_[1-8]_39-54_1_0_0_0_.1_5_1.mat';

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end

paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);

% % % % % only for CORNET NETWORK
for subji = 1:length(nnFit)
    nnH = nnFit{subji, 1}([2 4 6 8],:);
    nnFit{subji, 1} = nnH; 
end

clear hbL
for layi = 1:size(nnFit{1}, 1)
    clear nnH
    for subji = 1:length(nnFit)
       nnH(subji, : ,:) = nnFit{subji, 1}(layi,1:37);       
       %nnH(subji, : ,:) = nnFit{subji, 1}(layi,:);
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
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
     
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    hbL(layi, :) = hb; 
    
    mART(layi, :) = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    

end


if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 500/1.5 300/1.7])
    myCmap = colormap(brewermap(4,'RdBu'));
    times = 1:37;
    %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'});

    %x = (-.015:-.003:-.038)';
    %hbL([1:5], :) = nan; % VVS THETA
    %hbL(8, 35:end) = nan; % VVS THETA
    %hbL(1, :) = nan; % VVS ALPHA
    %hbL(:, 1:10) = nan; % VVS ALPHA
    %hbL(8, 1:20) = nan; % VVS ALPHA


    %x = (-.015:-.003:-.035)';
    %hbL([1:5 7], :) = nan; % VVS THETA
    %hbL(6, 1:20) = nan; % VVS THETA

    %hbL([1:4], 1:20) = nan; % VVS ALPHA
    %hbL([7], 21:45) = nan; % VVS ALPHA
    
    %hbL(1:7, :) = nan; % ALL
    
    %hbL(1:6, :) = nan; % PFC BETA


    % % % % cornet
    x = (-.015:-.003:-.025)';
    hbL([1:4], :) = nan; % VVS THETA, BETA, LOW-GAMMA, HIGH-GAMMA
    %hbL(:, 10:21) = nan; % VVS ALPHA
    %hbL([1:3], :) = nan; % PFC BETA
    

    
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(4,'RdBu'));
    %myCmap = colormap(jet(8));
    times = 1:21;
    %set(gca, 'ytick', [1 29 52], 'yticklabels', {'3' '30' '150'}); 
    %shadedErrorBar(times, mART, seART, 'r', 1); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    %x = (-.015:-.003:-.038)';
    %hbL([2:3], :) = nan; % VVS HIGH GAMMA ALEX
    %hbL([3:8], :) = nan; % PFC GAMMA ALEX
    %hbL([1:8], :) = nan; % ALL

    %x = (-.015:-.003:-.035)';
    %hbL([1 5:7], :) = nan; % VVS HIG GAMMA
    %hbL([1 1:6], :) = nan; % PFC THETA
    %hbL([1:7], :) = nan; % PFC GAMNMA = ALL 
    

    %CorNET
    x = (-.015:-.003:-.025)';
    hbL(1, 10:20) = nan; % VVS ALPHA

    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    %scatter(times, hb, 200,'.','Linewidth', 4)
    %set(gca, 'xtick', [1 5.5 15], 'xticklabels', {'-.5' '0' '1'}, 'xlim', [1 15])
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 12]); 
    %set(gca, 'FontSize', 12, 'ylim', [-.0375 .1]);
    set(gca, 'FontSize', 12, 'ylim', [-.0385 .1]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);
    %legend

end


set(gca, 'ColorOrder', myCmap)



exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% load files to plot all CATEGORY RESULTS  -- IN one PLOT ONLY 

clear, clc 

listF2sav = {       
    
                    'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';


%                     'CAT_pfc_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
%                     'CAT_pfc_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
%                     'CAT_pfc_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
%                     'CAT_pfc_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
%                     'CAT_pfc_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';


%                     'CAT_vvs_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';

%                     'CAT_vvs_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';


             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi allFnnH
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI);
    filelistSess = getFiles(paths.results.DNNs);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.results.DNNs listF2sav{listi}]);   
        
        nnH = cat(1, nnFit{:,1});
        
    end

    allFnnH(listi, :, :) = nnH; 


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
        nnH = allFnnH(freqi, :, 1:37);
    else
        nnH = allFnnH(freqi, :, 1:11);
    end
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
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
    set(gcf, 'Position', [100 100 500/1.5 300/1.7])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    times = 1:37;
    x = (-.015:-.003:-.028)';
    hbL([1:5], :) = nan;     
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    %colormap(jet(5));
    times = 1:11;
    plot(times, mART, 'Linewidth', 5); hold on; 
    x = (-.015:-.0053:-.039)';
    hbL(4, 1:3) = nan; 

    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 5); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 11]); 
    set(gca, 'FontSize', 12, 'ylim', [-.04 .15]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);
    %legend

end


set(gca, 'ColorOrder', myCmap)



exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% 

figure()
for layi = 1:56
    d2p = squeeze(mean(nnH(:, layi,:), 'omitnan'));
    plot(d2p); hold on; 
end
%h(h==0) = nan; h(h==1) = .02;
%plot(h, 'lineWidth', 2)




%% PERMUTATIONS IN LOOP
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

clear
nPerm = 1000;

listF2sav = {

%                 'AlexEco_pfc_M123_[1-8]_9-12_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_pfc_M123_[1-8]_13-29_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_pfc_M123_[1-8]_30-38_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_pfc_M123_[1-8]_39-54_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_vvs_M123_[1-8]_3-8_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_vvs_M123_[1-8]_9-12_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_vvs_M123_[1-8]_13-29_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_vvs_M123_[1-8]_30-38_1_0_0_0_.1_5_1.mat';
%                 'AlexEco_vvs_M123_[1-8]_39-54_1_0_0_0_.1_5_1.mat';

                  'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1'; 
                  'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1'; 

             };
    

for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi nPerm 
        
    f2sav       = listF2sav{listi}
        
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI);
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
        
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts.oneListPow);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts.oneListIds_c, sessi, paths);
    
    
            % % % restrict time for permutation data
            if strcmp(cfg.period(1), 'M')
                if ndims(neuralRDMs) == 4
                    neuralRDMs = neuralRDMs(:,:,:,3:37); %frequency-resolved
                else
                    neuralRDMs = neuralRDMs(:,:,3:37); %band analysis
                end
            else
                if ndims(neuralRDMs) == 4
                    neuralRDMs = neuralRDMs(:,:,:,3:12); %frequency-resolved
                else
                    neuralRDMs = neuralRDMs(:,:,3:12); %band analysis
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
    save([paths.results.DNNs f2sav '_' num2str(nPerm) 'p.mat'], 'nnFitPerm');
    

end
   
t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% COMPUTE TOBS FOR BANDS (ALL LAYERS RNN and Alex) 
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf

clear, clc 
%f2sav =    'RNN_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 
f2sav = 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';

cfg = getParams(f2sav); if strcmp(cfg.brainROI, 'vvs') sub2exc = [18 22]; elseif strcmp(cfg.brainROI, 'pfc')sub2exc = [1]; end
paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);

clear hbL
for layi = 1:size(nnFit{1}, 1)
        
    clear nnH
    for subji = 1:length(nnFit)
        if strcmp(cfg.period(1), 'M')
           nnH(subji, : ,:) = nnFit{subji, 1}(layi,1:37);       
        else
            nnH(subji, : ,:) = nnFit{subji, 1}(layi,1:21);       
        end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
            allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(layi, :, :) = 0;
    end
    
end






%% compute clusters in each permutation BANDS for all layers
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%f2sav =  'RNN_pfc_M123_[8-8-56]_13-29_1_1_0_0_.1_5_1_1000p.mat'; 
%f2sav = 'CORrt_vvs_E123_[1-8]_39-54_1_0_0_0_.1_5_1_1000p.mat';
f2sav = 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1_1000p.mat';



cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);
cd ([paths.results.DNNs])
load(f2sav);

if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end

f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{13}(1:end-5))));
if ndims (nnFitPerm) == 3
    nnFitPerm1(:, :, 1,:) = nnFitPerm;
    nnFitPerm = nnFitPerm1;

end

for layi = 1:size(nnFitPerm, 3)
    for permi = 1:nPerm
        
        dataP = squeeze(nnFitPerm(permi, :,layi, :));
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



%% MULTI - ITEM TRIALS

clear
%...__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials;1:Trials)__timeRes__win-width__mf_FST
f2sav = 'BLnext2_vvs_MALL_[6]_3-54_0_0_1_0_.1_5_1.mat'; 
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI);
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
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts.oneListIds_c, sessi, paths);
        
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    end
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');







%%  plot all layers MULTI-ITEM
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials;1:Trials)__timeRes__win__mf_FST
clear 
f2sav = 'BLnext2_pfc_MALL_[3]_3-54_0_0_1_0_.1_5_1.mat'; 

cfg = getParams(f2sav);


if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2]
end
 

paths = load_paths_WM(cfg.brainROI);
load([paths.results.DNNs f2sav]);


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
        [networkRDMs ids2rem]       = createNetworkRDMs(cfg, cfg_contrasts.oneListIds_c, sessi, paths);
        
        neuralRDMs(ids2rem, :, : ,:) = []; 
        neuralRDMs(:,ids2rem, : ,:) = []; 
        
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
        
    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav], 'nnFit');
    

end







%% STIMULUS REPRESENTATION IN DNNs
% % ALEXNET 
% % % 
clear, clc
f2sav = 'Alex_pfc_E123_[1-8]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM(cfg.brainROI);
[ACT] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


%% Plot all layers Alexnet one line horizontal

figure(); set(gcf, 'Position', [100 100 1500 500]);
myCmap = brewermap([], 'Spectral') 
colormap(myCmap)

for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [-.2 .6])
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

%% plot MDS one Line Vertical
cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);
%plot(1,1,'.','color',cols(11,:), 'Markersize', 2000) % check the color
% also nice palette here: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5a6dcde-4a80-11e4-9553-005056977bd0/a64ed616-e9ce-4b1d-8d18-5118cc03f8d2/images/screenshot.png

figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
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
 legend
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Representational consistency all layers / time points Alexnet


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  


allMS = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS(9,:) = nan; allMS(:,9) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 


% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0 .9], 'xtick',  (1:8) +.5, 'ytick', (1:8) +.5, ...
                'xticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'}, ...
                'yticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer Alexnet

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
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
plot (CCI, 'Linewidth', 3); axis square
set(gca, 'FontSize', 25, 'xlim', [0 9], 'ylim', [0 .75])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%  permutations
clearvars -except ACT CCI

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
    
    for layi = 1:size(ACT, 1)
        d2p = squeeze(ACT(layi, :,:)); 
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
lay = 8; 
obsT = CCI(lay); 
permT = CCIP(:, lay);


figure
histogram(permT); hold on; 
scatter(obsT,0, 'filled','r');



%% % RNN 
% % % 
clear, clc
f2sav = 'RNN_pfc_E123_[1-56]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM(cfg.brainROI);
[ACT] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
%[ACT] = load_rnn_eco(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet

%% correlate imagenet and ecoset activations 

for layi = 1:56

    rdmECO = squeeze(ACT_ECO(layi, :, :)); 
    rdmECO = vectorizeRDM(rdmECO);
    rdmIMA = squeeze(ACT_IMA(layi, :, :)); 
    rdmIMA = vectorizeRDM(rdmIMA);
    allRs(layi, :) = corr(rdmECO', rdmIMA', 'type', 's');

end 

%% 
figure()
plot(allRs, 'LineWidth', 3)
xlabel('Layer/Time')
ylabel('Rho')
set(gca, 'FontSize', 14)

%%

figure(); 
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';
lw = 3;
plot (allRs(1:8), 'Linewidth', lw, 'Color', cols(1, :)); hold on; 
plot (allRs(9:16), 'Linewidth', lw, 'Color', cols(2, :)); hold on; 
plot (allRs(17:24), 'Linewidth', lw, 'Color', cols(3, :)); hold on; 
plot (allRs(25:32), 'Linewidth', lw, 'Color', cols(4, :)); hold on; 
plot (allRs(33:40), 'Linewidth', lw, 'Color', cols(5, :)); hold on; 
plot (allRs(41:48), 'Linewidth', lw, 'Color', cols(6, :)); hold on; 
plot (allRs(49:56), 'Linewidth', lw, 'Color', cols(7, :)); hold on; 
set(gca, 'FontSize', 16, 'xlim', [1 8])
xlabel('Time')
ylabel('Rho')

%% RNN all RDMS

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
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
   d2p = squeeze(ACT(layi, :,:)); 
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

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
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

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
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
set(gca, 'FontSize', 25, 'xlim', [0 9], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for the last time point in each layer 

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
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


%figure()
%plot (CCI, 'Linewidth', 2)
%set(gca, 'FontSize', 20)

figure
lw = 4;
plot (CCI([8:8:56]), 'Linewidth', lw); hold on; 
set(gca, 'FontSize', 25, 'xlim', [.5 7.5], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);



%% compute CCI for each layer / timepoint FOR ECO AND IMAGE

clc
clearvars -except ACT_ECO ACT_IMA

%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(ACT_ECO, 1)
    d2p_ECO = squeeze(ACT_ECO(layi, :,:)); 
    d2p_ECO = 1- d2p_ECO;
    rdmMDS_ECO = d2p_ECO; 
    rdmMDS_ECO(find(eye(size(rdmMDS_ECO)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS_ECO(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS_ECO(M == 0), 'all', 'omitnan');
    CCI_ECO(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    %CCI_ECO(layi) = (mAcross - mWithin) ;
    
    d2p_IMA = squeeze(ACT_IMA(layi, :,:)); 
    d2p_IMA = 1- d2p_IMA;
    rdmMDS_IMA = d2p_IMA; 
    rdmMDS_IMA(find(eye(size(rdmMDS_IMA)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS_IMA(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS_IMA(M == 0), 'all', 'omitnan');
    CCI_IMA(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    %CCI_IMA(layi) = (mAcross - mWithin) ;
    
    
end

%% 

figure()
%plot (CCI_ECO, 'Linewidth', 2); hold on; 
%plot (CCI_IMA, 'Linewidth', 2); hold on; 
plot(CCI_ECO-CCI_IMA, 'Linewidth', 2); hold on; 
set(gca, 'FontSize', 20)




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





%% plot only last time point one line

figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 .8])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
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

 myCmap = brewermap([], '*Spectral') 
colormap (myCmap)

exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% plot only last time point one line vertical

figure(); set(gcf, 'Position', [100 100 700 700]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
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





%% plot MDS only last time point one line
cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 56];
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
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

%% MDS all layers (black and yellow circles plot)

act = ACT; 
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

act = ACT; 
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

act = ACT; 
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


act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 

% % % matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse');
set(gca, 'FontSize', 15, 'clim', [0 1], 'xtick',  (1:7) +.5, 'ytick', (1:7) +.5, 'xticklabels', {'1', '2', '3', '4', '5', '6', '7', '8'}, ...
                'yticklabels', {'1', '2', '3', '4', '5', '6', '7', '8'})

set(gca, 'FontSize', 25)
colormap(myCmap)
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



%% Changes in representational consistency 

act = ACT; 
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
set(gca, 'ylim', [0.8 1], 'xlim', [0 8], 'xtick', [1:7], 'xticklabels', {'1-2' '2-3' '3-4' '4-5' '5-6' '6-7' '7-8'}, 'FontSize', 25)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% 
% % CORNET 
% % % 
clear, clc
f2sav = 'CORrt_pfc_E123_[1-8]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM(cfg.brainROI);
[ACT] = load_CORrt_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet

%% Plot all layers CORNET one line horizontal

figure(); set(gcf, 'Position', [100 100 1500 500]);
myCmap = brewermap([], '*Spectral') 
colormap(myCmap)

for layi = [2 4 6 8]
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 .8])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot MDS one Line Vertical
cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);
%plot(1,1,'.','color',cols(11,:), 'Markersize', 2000) % check the color
% also nice palette here: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5a6dcde-4a80-11e4-9553-005056977bd0/a64ed616-e9ce-4b1d-8d18-5118cc03f8d2/images/screenshot.png

figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = [2 4 6 8]
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
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
 
 %legend
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Representational consistency all layers / time points CORNET


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  


allMS = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS(9,:) = nan; allMS(:,9) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 


% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0 .9], 'xtick',  (1:8) +.5, 'ytick', (1:8) +.5, ...
                'xticklabels', {'V1 in', 'V1 out', 'V2 in', 'V2 out', 'V4 in', 'V4 out', 'IT in', 'IT out'}, ...
                'yticklabels', {'V1 in', 'V1 out', 'V2 in', 'V2 out', 'V4 in', 'V4 out', 'IT in', 'IT out'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% Representational consistency all layers / ONLY 4 TIME points (OUTPUT) CORNET


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act1 = act1([2 4 6 8]);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 4, []); act3 = act2(:,all(~isnan(act2)));  


allMS = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS(5,:) = nan; allMS(:,5) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 


% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0 .9], 'xtick',  (1:8) +.5, 'ytick', (1:8) +.5, ...
                'xticklabels', {'V1', 'V2', 'V4', 'IT'}, ...
                'yticklabels', {'V1', 'V2', 'V4', 'IT'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer Cornet

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
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
plot (CCI([2 4 6 8]), 'Linewidth', 3); axis square
set(gca, 'FontSize', 25, 'xlim', [0 5], 'ylim', [0 .5], 'xtick', [1 2 3 4], 'xticklabels', {'V1' 'V2' 'V4' 'IT'})
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%  permutations
clearvars -except ACT CCI

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
    
    for layi = 1:size(ACT, 1)
        d2p = squeeze(ACT(layi, :,:)); 
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
lay = 2; 
obsT = CCI(lay); 
permT = CCIP(:, lay);


figure
histogram(permT); hold on; 
scatter(obsT,0, 'filled','r');



%% MDS all layers (pink and blue triangle-circles plot) only for first and last time points


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  
allM = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS = tril(allM, -1); allM(allM==0) = nan; 

clear c1 c2 c3 cols
c1 = (1:4)/4'; % sorted by time point
c2 = repmat(zeros(1), 4, 1)';
c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);

d2p = 1- allMS;
[rdmMDS] = cmdscale(allM);

figure()
plot(rdmMDS([1 2],1),rdmMDS([1 2],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([3 4],1),rdmMDS([3 4],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([5 6],1),rdmMDS([5 6],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([7 8],1),rdmMDS([7 8],2),'k', 'linewidth', 2);hold on; 

fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 8);

for i = 1:8
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gca,'FontSize', 30, 'xlim', [-.35 .55], 'ylim', [-.25 .5]); 
end



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);



















