%% VVS PFC
clear 

paths = load_paths_WM('vvs'); 

%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
load ([paths.trial_level 'RNN_pfc_M123_[8-8-56]_13-29_0_0_0_1_.1_5_1.mat']); 
pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
load ([paths.trial_level 'RNN_vvs_M123_[8-8-56]_9-12_0_0_0_1_.1_5_1.mat']); 
vvs_fits = nnFit([7 9 13 18 19 20 21 23 27 28]); 
load ([paths.electrodes_path 'pfc_elec']); 
pfc = pfc_fits; 
vvs = vvs_fits; 


%% bands
pfc_fits = pfc; 
vvs_fits = vvs; 

tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp([2 3  5  9 10 11 12 14 15 16],1) > 3); 
%sub2exc = []; 
pfc_fits(s2e_pfc) = [];vvs_fits(s2e_pfc) = [];


tt1 = 3:8; % taken from pfc band data (13-29Hz)
tt2 = 21:35;



nLays = 7;
clear cR

for subji = 1:length(pfc_fits)
   
    x = pfc_fits{subji};
    y = vvs_fits{subji};
    
%     figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:nLays % 7:7 %
        x1 = squeeze(mean(x(layi, :, tt1), 3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(y(layj, :, tt2), 3, 'omitnan'))';
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
%             subplot(7 ,7, ploti)
%             % % scatterHistDiff(x1, y1);
%             scatter(x1, y1, 50, '.'); hold on; 
%             p = polyfit(x1,y1,1); 
%             f = polyval(p,x1); 
%             plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
%             
%             ploti = ploti+1;
            
            
        end
    end    
    
   % exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
   % close all
end

cR = atanh(cR)
[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure()
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22) %, 'clim', [-4 4]

%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


%exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)








%% check imagesc

M = [ 1 2 3 ; 4 5 6; 7 8 9];
figure
imagesc(M)
figure
contourf(myresizem(M, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar





%% frequency resolved

%% VVS PFC
clear 

paths = load_paths_WM('vvs'); 
%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
load ([paths.trial_level 'RNN_pfc_M123_[8-8-56]_3-54_0_0_0_1_.1_5_1.mat']); 
pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
load ([paths.trial_level 'RNN_vvs_M123_[8-8-56]_3-54_0_0_0_1_.1_5_1.mat']); 
vvs_fits = nnFit([7 9 13 18 19 20 21 23 27 28]); 
load ([paths.electrodes_path 'pfc_elec']); 
pfc = pfc_fits; 
vvs = vvs_fits; 



%% frequnecy resolved
pfc_fits = pfc; 
vvs_fits = vvs; 

tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp([2 3  5  9 10 11 12 14 15 16],1) > 2); 
%sub2exc = []; 
pfc_fits(s2e_pfc) = [];vvs_fits(s2e_pfc) = [];


ff1 = 18:28; 
tt1 = 21:26; 
ff2 = 30:38;
tt2 = 41:60; 


nLays = 7;
clear cR

for subji = 1:length(pfc_fits)
   
    x = pfc_fits{subji};
    y = vvs_fits{subji};
    
%     figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:7 %5:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'), 3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'), 3, 'omitnan'))'; 
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
%             subplot(7 ,7, ploti)
%             % % scatterHistDiff(x1, y1);
%             scatter(x1, y1, 50, '.'); hold on; 
%             p = polyfit(x1,y1,1); 
%             f = polyval(p,x1); 
%             plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
%             
%             ploti = ploti+1;
            
            
        end
    end    
    
   % exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
   % close all
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure(); %set(gcf, 'Position', [1000 1000 500 200])
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); hold on; colorbar; axis equal
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22, 'clim', [-4 4])
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)



%% VVS Hippocampus
clear 

paths = load_paths; 

%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
load ([paths.neural_rdms_path 'HIP_8-16-24-32-40-48-56_01-54_1_1_0.1_5_1_.mat']); 
hip_fits = nnFit(1:16);
load ([paths.neural_rdms_path 'VVS_8-16-24-32-40-48-56_01-54_1_1_0.1_5_1_.mat']); 
vvs_fits = nnFit([1 5 10 12:15 17 18 21 23:28]);
load ([paths.electrodes_path 'hipp_elec']); 
hip = hip_fits; 
vvs = vvs_fits; 


%% frequency resolved
hip_fits = hip; 
vvs_fits = vvs; 

tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_hip = ~(tmp([1:16],1) > 2); 
%sub2exc = []; 
hip_fits(s2e_hip) = [];vvs_fits(s2e_hip) = [];

ff1 = 18:28; 
tt1 = 21:31; 
ff2 = 39:54;
tt2 = 51:60; 

nLays = 7;
clear cR

for subji = 1:length(hip_fits)
   
    x = hip_fits{subji};
    y = vvs_fits{subji};
    
%     figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:3 %5:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi+4, :, ff1, tt1), 4, 'omitnan'), 3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'), 3, 'omitnan'))'; 
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
%             subplot(7 ,7, ploti)
%             % % scatterHistDiff(x1, y1);
%             scatter(x1, y1, 50, '.'); hold on; 
%             p = polyfit(x1,y1,1); 
%             f = polyval(p,x1); 
%             plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
%             
%             ploti = ploti+1;
            
            
        end
    end    
    
   % exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
   % close all
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure(); set(gcf, 'Position', [1000 1000 500 200])
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); hold on; colorbar; %axis equal
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[5:7]})
set(gca, 'FontSize', 22, 'clim', [-4 4])
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)







%%
clear 
%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
load PFC_8-16-24-32-40-48-56_13-29_0_1_0.01_50_1_.mat
pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
load VVS_8-16-24-32-40-48-56_30-38_0_1_0.01_50_1_.mat
vvs_fits = nnFit([7 9 13 18 19 20 21 23 27 28]); 
load pfc_elec
pfc = pfc_fits; 
vvs = vvs_fits; 











%%