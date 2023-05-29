%% first load traces
%%
clear
region = 'pfc';
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

    %cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    %cfg_contrasts = normalize_WM(cfg_contrasts, 0, 'sess', [51 100]);
    
    allPow{sessi} = cfg_contrasts; 
   
end

mkdir ([paths.results.power]);
save([paths.results.power 'pow2S_' region '_' num2str(cfg.timeRes*1000) 'ms'], 'allPow');


%% 
clear
paths = load_paths_WM('pfc');
load([paths.results.power 'pow2S_pfc_10ms'], 'allPow');

%% DIRECT COMPARISON all trials
 
clear
paths = load_paths_WM('pfc');
load([paths.results.power 'pow2S_vvs_10ms'], 'allPow');
allPow_VVS = allPow; 
load([paths.results.power 'pow2S_pfc_10ms'], 'allPow');
allPow_PFC = allPow; 



%% plot example trial 

size(allPow{1}.oneListPow)

d2p = squeeze(allPow{1}.oneListPow(1, 1,:,:));

figure
imagesc(d2p)


%% compute mean over trials and over electrodes for each subject

clearvars -except allPow_VVS allPow_PFC

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];


allT_VVS = allPow_VVS(vvs_ids)'; 
allT_PFC = allPow_PFC(pfc_ids)'; 

clear powSI_VVS powMI_VVS powSI_PFC powMI_PFC powMI_PFC_TR powMI_VVS_TR
for subji = 1:length(allT_PFC)
    ids = cellfun(@(x) strsplit(x), allT_PFC{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));
    powSI_PFC(subji, :, :,:) = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    powMI_PFC(subji, :, :,:) = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids == 4,: ,:,:), 2))); 

    powSI_VVS(subji, :, :,:) = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    powMI_VVS(subji, :, :,:) = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids == 4,: ,:,:), 2))); 


end


%% plot one example subject


d2p = squeeze(powSI_PFC(4,:,:));
figure
imagesc(d2p)


%% plot mean across subjects

powSI = powSI_PFC(:, 3:54,:);
powMI = powMI_PFC(:, 3:54,:);

mSI = squeeze(mean(powSI)); 
mMI = squeeze(mean(powMI)); 


%times = 1:550; 
times = -1:0.01:4.49;
freqs = 1:52; 

figure;
myCmap = colormap(brewermap([],'*Spectral'));
contourf(times, freqs, mSI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 27 52], 'yticklabels', {'3' '30' '150'})
%set(gca, 'clim', [-3 3], 'ytick', [1 29 54], 'yticklabels', {'1' '30' '150'})
set(gca, 'FontSize', 14)
title('Single-item trials')
colormap(myCmap)

figure()
contourf(times, freqs, mMI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 27 52], 'yticklabels', {'3' '30' '150'})
%set(gca, 'clim', [-3 3], 'ytick', [1 29 54], 'yticklabels', {'1' '30' '150'})
set(gca, 'FontSize', 14)
title('Multi-item trials')
colormap(myCmap)




%% statistics 


[h p ci ts] = ttest(powSI, powMI);
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

[max_clust_obs id] = max(abs(allTOBs));

%h = zeros(52, 550);
%h(clustinfo.PixelIdxList{id})= 1; 


figure()
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'clim', [-4 4])


%% plot all together

times = -1:0.01:4.49;
freqs = 1:52;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, mSI, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single-item trials')
nexttile
contourf(times, freqs, mMI, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.1 .1])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Multi-item trials')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-5 5])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 27 52], 'yticklabels', {'3', '30', '150'}, 'xlim', [-.5 3.5]);


exportgraphics(gcf, ['myP.png'], 'Resolution',300)



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



%% trial level correlations between regions MI trials
load clust_vvs_pfc_52.mat
%load clusts_vvs_pfc.mat

clearvars -except allPow_VVS allPow_PFC clust_vvs clust_pfc

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

allT_VVS = allPow_VVS(vvs_ids)'; 
allT_PFC = allPow_PFC(pfc_ids)'; 

% remove intersect
inter2u = intersect(clust_vvs, clust_pfc);
%inter2u = 9500:11000;
clust_vvs = setdiff(clust_vvs, inter2u);
clust_pfc = setdiff(clust_pfc, inter2u);
% % just2check
%intersect(clust_pfc, inter2u)
% figure()
% h= zeros(54, 550); 
% h(clust_vvs) = 1;
% contourf(h);

clear powSI_VVS powMI_VVS powSI_PFC powMI_PFC powMI_PFC_TR powMI_VVS_TR
for subji = 1:length(allT_PFC)
    ids = cellfun(@(x) strsplit(x), allT_PFC{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));

    %pfc2cmp = squeeze(mean(allT_PFC{subji}.oneListPow(ids == 4,: ,:,:), 2)); 
    %vvs2cmp = squeeze(mean(allT_VVS{subji}.oneListPow(ids == 4,: ,:,:), 2)); 

    pfc2cmp = squeeze(mean(allT_PFC{subji}.oneListPow(:,:, 3:54,:), 2)); 
    vvs2cmp = squeeze(mean(allT_VVS{subji}.oneListPow(:,:, 3:54,:), 2)); 

    % z-score the data separately for each condition
    pfc2cmpSI = pfc2cmp(ids~=4,:,:); 
    pfc2cmpSIN = normalize_within_cond_WM(pfc2cmpSI);
    pfc2cmpMI = pfc2cmp(ids==4,:,:); 
    pfc2cmpMIN = normalize_within_cond_WM(pfc2cmpMI);
    vvs2cmpSI = vvs2cmp(ids~=4,:,:); 
    vvs2cmpSIN = normalize_within_cond_WM(vvs2cmpSI);
    vvs2cmpMI = vvs2cmp(ids==4,:,:); 
    vvs2cmpMIN = normalize_within_cond_WM(vvs2cmpMI);
    
    pfc2cmp(ids~= 4, :, :) = pfc2cmpSIN;
    pfc2cmp(ids== 4, :, :) = pfc2cmpMIN;
    vvs2cmp(ids~= 4, :, :) = vvs2cmpSIN;
    vvs2cmp(ids== 4, :, :) = vvs2cmpMIN;

    clear valTrVVS valTrPFC
    for triali = 1:size(vvs2cmp, 1)
        x = vvs2cmp(triali, clust_vvs);
        y = pfc2cmp(triali, clust_pfc);
        valTrVVS(triali,:) = mean(x);
        valTrPFC(triali,:) = mean(y);
    end

    powMI_PFC_TR{subji,:} = valTrVVS;
    powMI_VVS_TR{subji,:} = valTrPFC;

    figure()
    x1 = valTrVVS; y1 = valTrPFC;
    scatter(x1, y1, 50, '.'); hold on; 
    p = polyfit(x1,y1,1); 
    f = polyval(p,x1); 
    plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 


    allRho(subji, :) = corr(valTrVVS, valTrPFC, 'type','s');

end

[h p ci ts] = ttest(allRho);
t = ts.tstat; 

disp (['t = ' num2str(t) '  ' ' p = ' num2str(p)]);

%% plot one bar
data.data = [allRho]; 


figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.15 .5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   



%% RESHAPE DATA FOR LME 


pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];
allT_VVS = allPow_VVS(vvs_ids)'; 
allT_PFC = allPow_PFC(pfc_ids)'; 

clear pow_PFC_TR pow_VVS_TR ids2u
for subji = 1:length(allT_PFC)
    ids = cellfun(@(x) strsplit(x), allT_PFC{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));
    pfc2cmp = squeeze(mean(allT_PFC{subji}.oneListPow(:,:, 3:54,:), 2)); 
    vvs2cmp = squeeze(mean(allT_VVS{subji}.oneListPow(:,:, 3:54,:), 2)); 

    % % % % z-score the data separately for each condition
    pfc2cmpSI = pfc2cmp(ids~=4,:,:); 
    pfc2cmpSIN = normalize_within_cond_WM(pfc2cmpSI);
    pfc2cmpMI = pfc2cmp(ids==4,:,:); 
    pfc2cmpMIN = normalize_within_cond_WM(pfc2cmpMI);
    vvs2cmpSI = vvs2cmp(ids~=4,:,:); 
    vvs2cmpSIN = normalize_within_cond_WM(vvs2cmpSI);
    vvs2cmpMI = vvs2cmp(ids==4,:,:); 
    vvs2cmpMIN = normalize_within_cond_WM(vvs2cmpMI);
    
    pfc2cmp(ids~= 4, :, :) = pfc2cmpSIN;
    pfc2cmp(ids== 4, :, :) = pfc2cmpMIN;
    vvs2cmp(ids~= 4, :, :) = vvs2cmpSIN;
    vvs2cmp(ids== 4, :, :) = vvs2cmpMIN;


    clear valTrVVS valTrPFC
    for triali = 1:size(vvs2cmp, 1)
        x = vvs2cmp(triali, clust_vvs);
        y = pfc2cmp(triali, clust_pfc);
        valTrVVS(triali,:) = mean(x);
        valTrPFC(triali,:) = mean(y);
    end

    pow_PFC_TR{subji,1} = valTrVVS;
    pow_VVS_TR{subji,1} = valTrPFC;
    ids2u{subji,1} = ids; 
    


end



%%
clc
clear tbl 

lenSj = cellfun(@length, pow_VVS_TR);
subj2u = [ones(1, lenSj(1)), (ones(1, lenSj(2)) +1), (ones(1, lenSj(3)) +2), (ones(1, lenSj(4)) +3),(ones(1, lenSj(5)) +4), ...
    (ones(1, lenSj(6)) +5),(ones(1, lenSj(7)) +6),(ones(1, lenSj(8)) +7),(ones(1, lenSj(9)) +8),(ones(1, lenSj(10)) +9)]';

bb1 = [cat(1, pow_VVS_TR{:})];
bb2 = [cat(1, pow_PFC_TR{:})];

Y = cat(1, bb1, bb2)
idsModel = [cat(1, ids2u{:})];
Trial_Type = idsModel; 
Trial_Type(idsModel<4 ) = 1; Trial_Type(idsModel==4 ) = 2; 

tbl = [bb1 , bb2, Trial_Type, subj2u];
tbl2 = table(tbl(:,1), tbl(:,2), tbl(:,3), tbl(:,4),'VariableNames',{'theta_VVS','theta_PFC','Trial_Type', 'Subj'});




%% fit model
clc

%lme = fitlme(tbl2,'theta_VVS ~ 1+ theta_PFC');
%lme = fitlme(tbl2,'theta_VVS ~ theta_PFC'); % first two are the same (included by default)
%lme = fitlme(tbl2,'theta_VVS~ theta_PFC+(1|Subj)'); 
%lme = fitlme(tbl2,'theta_VVS~theta_PFC+(1+theta_PFC|Subj)');
%lme = fitlme(tbl2,'theta_VVS~theta_PFC+(1|Subj)+(1|Trial_Type)'); 
lme = fitlme(tbl2,'theta_VVS~theta_PFC+Trial_Type+(1|Subj)');
%lme = fitlme(tbl2,'theta_VVS~theta_PFC+Trial_Type+(1+theta_PFC|Subj)');
%lme = fitlme(tbl2,'theta_VVS~theta_PFC+(1+theta_PFC|Subj)+(1+theta_PFC|Trial_Type)');
%lme = fitlme(tbl2,'theta_PFC~theta_VVS+(1+theta_VVS|Subj)+(1+theta_VVS|Trial_Type)');

lme
%lme.Coefficients

%%
compare(lme2, lme3)


%% 
mean(Y(Trial_Type == 1))
mean(Y(Trial_Type == 2))
figure()

histogram(Y(Trial_Type == 1)); hold on; 
histogram(Y(Trial_Type == 2)); 

x = Y(Trial_Type == 1); 
y = Y(Trial_Type == 2);
[h p ci ts] = ttest2(x,y)






















%% DIRECT COMPARISON all trials
 
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


%% 

clear powDIFF_VVS powDIFF_PFC  
for subji = 1:length(allT_VVS)


    
    ids = cellfun(@(x) strsplit(x), allT_VVS{subji}.oneListIds_c, 'un', 0);
    ids = cell2mat(cellfun(@(x) double(string(x(2))), ids, 'un', 0));

    t1 = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    t2 = squeeze(mean(mean(allT_VVS{subji}.oneListPow(ids == 4,: ,:,:), 2))); 
    powDIFF_VVS(subji, :, :)  = t1-t2; 
    
    t1 = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids ~= 4,: ,:,:), 2))); 
    t2 = squeeze(mean(mean(allT_PFC{subji}.oneListPow(ids == 4,: ,:,:), 2))); 
    powDIFF_PFC(subji, :, :)  = t1-t2; 


end




%% comparison VVS PFC for all trials


m_Diff= powDIFF_PFC - powDIFF_VVS; 
m_VVS = squeeze(mean(powDIFF_VVS)); 
m_PFC = squeeze(mean(powDIFF_PFC)); 


[h p ci ts] = ttest(m_Diff);
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allTOBs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));%
end

max_clust_obs = max(abs(allTOBs));

times = -1:0.01:4.49
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, m_VVS, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.2 .2])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single minus multi-item VVS')
nexttile
contourf(times, freqs, m_PFC, 40, 'linecolor', 'none'); hold on; colorbar
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.2 .2])
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
title('Single minus multi-item PFC')
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