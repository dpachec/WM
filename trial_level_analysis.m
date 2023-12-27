%% 
%% Trial level analysis PFC

clear
f2sav = 'BLNETi_pfc_M123_[56]_3-54_0_0_1_1_.1_5_1'; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);
load([paths.results.clusters 'clustinfo_PFC_px2.mat']);


for subji = 1:length(nnFit)


    nnF = nnFit{subji, 1}; 
    for triali = 1:size(nnF, 2)
        nnFT = squeeze(nnF(1, triali, :, 1:40)); 
        % % % check that cluster is correct
%         times = 1:400; freqs = 1:520; 
%         contourf(times, freqs, myresizem(nnFT, 10), 100, 'linecolor', 'none'); hold on; colorbar
%         h = zeros(52,40);
%         h(clustinfo.PixelIdxList{2}) = 1; 
%         contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        fitTR(triali, :) = mean(nnFT(clustinfo.PixelIdxList{2}), 'all');
    end

    ids = nnFit{subji, 2}; 
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids = double(string(cellfun(@(x) x(9), ids, 'UniformOutput', false)));

    fTRC(subji, :) = mean(fitTR(ids==1)); 
    fTRI(subji, :) = mean(fitTR(ids==0)); 




end

fTRC(sub2exc) = []; 
fTRI(sub2exc) = []; 

%% 


data.data = [fTRC fTRI];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.1 .2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
%[h p ci t] = ttest (data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)


%% VVS

clear

%f2sav = 'BLNETi_pfc_M123_[56]_3-54_0_0_1_1_.1_5_1'; 
f2sav = 'BLNETi_vvs_M123_[32 40 48]_3-54_0_0_1_1_.1_5_1'; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);
load([paths.results.clusters 'all_clustinfo_VVS.mat']);



for subji = 1:length(nnFit)


    nnF = nnFit{subji, 1}; 
    for triali = 1:size(nnF, 2)
        %nnFT = squeeze(nnF(1, triali, :, 1:40)); 
        %fitTR(triali, :) = mean(nnFT(allClustInfo{4}.PixelIdxList{14}), 'all');
        
        %nnFT = squeeze(nnF(2, triali, :, 1:40)); 
        %fitTR(triali, :) = mean(nnFT(allClustInfo{5}.PixelIdxList{25}), 'all');

        nnFT = squeeze(nnF(3, triali, :, 1:40)); 
        fitTR(triali, :) = mean(nnFT(allClustInfo{6}.PixelIdxList{17}), 'all');
    end

    ids = nnFit{subji, 2}; 
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids = double(string(cellfun(@(x) x(8), ids, 'UniformOutput', false)));

    fTRC(subji, :) = mean(fitTR(ids==1)); 
    fTRI(subji, :) = mean(fitTR(ids==0)); 

end


fTRC(sub2exc) = []; 
fTRI(sub2exc) = []; 


data.data = [fTRC fTRI];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.15 .2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
%[h p ci t] = ttest (data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)





















%%




























%%