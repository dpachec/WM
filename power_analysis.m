%% first load traces
%%
clear
region = 'vvs';
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);



for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   

    %ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    ids = cell2mat(cellfun(@(x) strcmp(x(1), '7'), cfg_contrasts.oneListIds_c, 'un', 0));
    cfg_contrasts.oneListTraces   = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 

    
    cfg = []; 
    cfg.timeRes = .01;  
    cfg.period = 'M'; 
    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 

    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    
    allPow{sessi} = cfg_contrasts; 
   
end

mkdir ([paths.results.power]);
save([paths.results.power 'pow2S_' region '_' num2str(cfg.timeRes*1000) 'ms'], 'allPow');


%% 
clear
paths = load_paths_WM('pfc');
load([paths.results.power 'pow2S_pfc_10ms'], 'allPow');

%% plot example trial 

size(allPow{1}.oneListPow)

d2p = squeeze(allPow{1}.oneListPow(1, 1,:,:));

figure
imagesc(d2p)


%% compute mean over trials and over electrodes for each subject

clear powSI powMI
for subji = 1:length(allPow)


    
    ids = cellfun(@(x) strsplit(x), allPow{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));

    
    powSI(subji, :, :,:) = squeeze(mean(mean(allPow{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    powMI(subji, :, :,:) = squeeze(mean(mean(allPow{subji}.oneListPow(ids == 4,: ,:,:), 2))); 
    


end


%% plot one example subject


d2p = squeeze(powSI(4,:,:));
figure
imagesc(d2p)


%% plot mean across subjects


mSI = squeeze(mean(powSI)); 
mMI = squeeze(mean(powMI)); 


%times = 1:550; 
times = -1:0.01:4.49
freqs = 1:54; 

figure;
myCmap = colormap(brewermap([],'*Spectral'));
contourf(times, freqs, mSI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 29 54], 'yticklabels', {'1' '30' '150'})
set(gca, 'FontSize', 14)
title('Single-item trials')
colormap(myCmap)

figure()
contourf(times, freqs, mMI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 29 54], 'yticklabels', {'1' '30' '150'})
set(gca, 'FontSize', 14)
title('Multi-item trials')
colormap(myCmap)




%% statistics 


[h p ci ts] = ttest(powSI, powMI)
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

max_clust_obs = max(abs(allTOBs))


figure()
times = 1:550; 
freqs = 1:54; 
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'clim', [-4 4])


%% plot all together

times = -1:0.01:4.49
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, mSI, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);%set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single-item trials')
nexttile
contourf(times, freqs, mMI, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); %set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Multi-item trials')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); %set(gca, 'clim', [-3 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',20, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});


%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = powSI(:,1:54,1:550); 
    c2B = powMI(:,1:54,1:550); 
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :, :);
           c1B(subji, :, :) = c2B(subji, :, :);
           c2B(subji, :, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1B, c2B); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = min(allSTs);
    else
        allSTs = 0; 
    end
    max_clust_sum_perm(permi,:) = allSTs(id); 

end

%%
clear p mcsR mcsP

mcsR = max_clust_obs; 
%mcsR  = 608.294938960190
mcsP = abs(max_clust_sum_perm);

allAb = mcsP(mcsP > mcsR);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm





%% 
figure
histogram(max_clust_sum_perm); hold on; 
scatter(max_clust_obs,0, 200, 'filled','r');
set(gca, 'FontSize', 14);
xlabel('T')
exportgraphics(gcf, [paths.results.power 'myP.png'], 'Resolution',150)






%% direct comparison all trials


%% 
clear
paths = load_paths_WM('pfc');
load([paths.results.power 'pow2S_vvs_10ms'], 'allPow');
allPow_VVS = allPow; 
load([paths.results.power 'pow2S_pfc_10ms'], 'allPow');
allPow_PFC = allPow; 


%%

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];


allT_VVS = allPow_VVS(vvs_ids)'; 
allT_PFC = allPow_PFC(pfc_ids)'; 



p_vvs = struct2cell(cat(1, allT_VVS{:}))';
p_vvs = p_vvs(:, 4);
p_vvs2 = cellfun(@(x) squeeze(mean(mean(x))), p_vvs, 'un', 0);
p_vvs2 = cat(3, p_vvs2{:}); 
p_vvs3 = permute(p_vvs2, [3 1 2]);



p_pfc = struct2cell(cat(1, allT_PFC{:}))';
p_pfc = p_pfc(:, 4);
p_pfc2 = cellfun(@(x) squeeze(mean(mean(x))), p_pfc, 'un', 0);
p_pfc2 = cat(3, p_pfc2{:}); 
p_pfc3 = permute(p_pfc2, [3 1 2]);


%% comparison VVS PFC for all trials


m_VVS = squeeze(mean(p_vvs3)); 
m_PFC = squeeze(mean(p_pfc3)); 


[h p ci ts] = ttest(p_vvs3, p_pfc3)
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

max_clust_obs = max(abs(allTOBs))

times = -1:0.01:4.49
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, m_VVS, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);%set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('All trials VVS')
nexttile
contourf(times, freqs, m_PFC, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); %set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('All trials PFC')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); %set(gca, 'clim', [-3 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',20, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});

exportgraphics(gcf, 'myP.png', 'Resolution',150)

%% compute separately for single and multi item trials 


clear powSI_VVS powMI_VVS powSI_PFC powMI_PFC
for subji = 1:10


    
    ids = cellfun(@(x) strsplit(x), allT_VVS{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));
    powSI_VVS(subji, :, :,:) = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    powMI_VVS(subji, :, :,:) = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids == 4,: ,:,:), 2))); 

    ids = cellfun(@(x) strsplit(x), allT_PFC{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));
    powSI_PFC(subji, :, :,:) = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    powMI_PFC(subji, :, :,:) = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids == 4,: ,:,:), 2))); 

    
end


%% plot mean across subjects SINGLE ITEM


mSI_VVS = squeeze(mean(powSI_VVS)); 
mSI_PFC = squeeze(mean(powSI_PFC)); 


[h p ci ts] = ttest(powSI_VVS, powSI_PFC)
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

max_clust_obs = max(abs(allTOBs))

times = -1:0.01:4.49
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, mSI_VVS, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);%set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single-item trials VVS')
nexttile
contourf(times, freqs, mSI_PFC, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); %set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single-item trials PFC')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); %set(gca, 'clim', [-3 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',20, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});


exportgraphics(gcf, 'myP.png', 'Resolution',150)


%% plot mean across subjects MULTI-ITEM


mMI_VVS = squeeze(mean(powMI_VVS)); 
mMI_PFC = squeeze(mean(powMI_PFC)); 


[h p ci ts] = ttest(powMI_VVS, powMI_PFC)
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

max_clust_obs = max(abs(allTOBs))

times = -1:0.01:4.49
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, mMI_VVS, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);%set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Multi-item trials VVS')
nexttile
contourf(times, freqs, mMI_PFC, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); %set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Multi-item trials PFC')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); %set(gca, 'clim', [-3 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',20, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});


exportgraphics(gcf, 'myP.png', 'Resolution',150)

%% 
clear pvvs
for subji = 1:10


    pvvs{subji} = allT_VVS{subji}.oneListPow; 



end





















%%