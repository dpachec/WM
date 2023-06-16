%% first load traces
%% power computation
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

    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    allPow{sessi} = cfg_contrasts; 
   
end

mkdir ([paths.results.power]);
save([paths.results.power 'pow2S_' region '_' num2str(cfg.timeRes*1000) 'ms'], 'allPow');


%% DIRECT COMPARISON all trials
 
clear
paths = load_paths_WM('pfc', []);
load([paths.results.power 'pow2S_vvs_10ms'], 'allPow');
allPow_VVS = allPow; 
load([paths.results.power 'pow2S_pfc_10ms'], 'allPow');
allPow_PFC = allPow; 


%% compute mean over trials and over electrodes for each subject

clearvars -except allPow_VVS allPow_PFC

pfc_ids = [2 3  5  9 10 11 12 14 15 16]; %ids of subject with electrodes in both regions
vvs_ids = [7 9 13 18 19 20 21 23 27 28]; % ids of subjects with electrodes in both regions


allT_VVS = allPow_VVS(vvs_ids)'; %% select VVS or PFC
allT_PFC = allPow_PFC(pfc_ids)'; %% select VVS or PFC

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

powSI = powSI_VVS(:, 3:54,:);
powMI = powMI_VVS(:, 3:54,:);

mSI = squeeze(mean(powSI)); 
mMI = squeeze(mean(powMI)); 


%times = 1:550; 
times = -1:0.01:4.49;
freqs = 1:52; 

figure;
myCmap = colormap(brewermap([],'*Spectral'));
contourf(times, freqs, mSI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 27 52], 'yticklabels', {'3' '30' '150'})
set(gca, 'FontSize', 14)
title('Single-item trials')
colormap(myCmap)

figure()
contourf(times, freqs, mMI, 100, 'linecolor', 'none'); hold on; colorbar
set(gca, 'clim', [-0.08 0.08], 'ytick', [1 27 52], 'yticklabels', {'3' '30' '150'})
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

h = zeros(52, 550);
h(clustinfo.PixelIdxList{id})= 1; 


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
colormap(brewermap([],'Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 27 52], 'yticklabels', {'3', '30', '150'}, 'xlim', [-.5 3.5]);


exportgraphics(gcf, ['myP.png'], 'Resolution',300)



%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = powSI(:,1:52,101:450); 
    c2B = powMI(:,1:52,101:450); 
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
mcsP = abs(max_clust_sum_perm);

allAb = mcsP(mcsP > mcsR);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm





%% 
figure
histogram(max_clust_sum_perm); hold on; 
scatter(max_clust_obs,0, 200, 'filled','r');
set(gca, 'FontSize', 14);
xlabel('T')


%% trial level correlations between regions MI trials

cd D:\_WM\analysis
load clust_vvs_pfc_52.mat %contains the two clusters identified previously
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

    clear valTrVVS valTrPFC
    for triali = 1:size(vvs2cmp, 1)
        x = vvs2cmp(triali, clust_vvs);
        y = pfc2cmp(triali, clust_pfc);
        valTrVVS(triali,:) = mean(x);
        valTrPFC(triali,:) = mean(y);
    end

    powMI_PFC_TR{subji,:} = valTrVVS;
    powMI_VVS_TR{subji,:} = valTrPFC;

    x1 = valTrVVS(ids<4); %% only for single or multi-item trials 
    y1 = valTrPFC(ids<4); %% only for single or multi-item trials 

    figure()
    scatter(x1, y1, 50, '.'); hold on; 
    p = polyfit(x1,y1,1); 
    f = polyval(p,x1); 
    plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 


    allRho(subji, :) = corr(x1, y1, 'type','s');

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

%export_fig(2, '_2.png','-transparent', '-r300');
%close all;   



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


tbl2.Trial_Type = nominal(tbl2.Trial_Type)

%% fit model
clc

lme = fitlme(tbl2,'theta_VVS~theta_PFC+Trial_Type+(1|Subj)'); % random intercept model
lme


%% 'caseorder' | 'fitted' | 'lagged' | 'probability' | 'symmetry'
%plotResiduals(lme)
%plotResiduals(lme,'histogram')
%plotResiduals(lme,'caseorder')
plotResiduals(lme,'fitted')
%plotResiduals(lme,'lagged')
%plotResiduals(lme,'probability')
%plotResiduals(lme,'symmetry')

%%
F = fitted(lme);
R = residuals(lme);

plot(F,R,'bx')
xlabel('Fitted Values')
ylabel('Residuals')


%%
F = fitted(lme);
R = response(lme);
figure();
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')
%%
figure()
gscatter(F,R,tbl2.Trial_Type)

%%
figure
ypred = predict(lme)
plot(tbl2.theta_VVS,tbl2.theta_PFC,'o',ypred,tbl2.theta_PFC, 'x')
legend('Data','Predictions')







%% plot two clusters ONLY, CARTOON FIGURE 1F
h = zeros(52, 550);
h1 = h; h2 = h; 
h1(clust_vvs)= 1;
h2(clust_pfc)= 2;
h = h1+h2; 
hc = h; hc(hc>0) = 1;

times = -1:0.01:4.49;
freqs = 1:52;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 400 800])
nexttile
imagesc(times, freqs, flipud(h)); hold on;
plot([.5 .5],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))
%set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 27 52], 'yticklabels', {'3', '30', '150'}, 'xlim', [-.5 3.5]);
set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 22 52], 'yticklabels', {'150', '30', '3'}, 'xlim', [0.25 1.6], 'ylim', [22 52]);
exportgraphics(gcf, ['myP.png'], 'Resolution',300)















%%