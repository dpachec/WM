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
load([paths.results.power 'pow2S_vvs_10ms'], 'allPow');

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
myCmap = colormap(brewermap([],'*Spectral'));

%times = 1:550; 
times = -1:0.01:4.49
freqs = 1:54; 

figure;
contourf(times, freqs, mSI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 29 54], 'yticklabels', {'1' '30' '150'})
set(gca, 'FontSize', 14)
title('Single-item trials')
colormap(myCmap)

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





%%























%%