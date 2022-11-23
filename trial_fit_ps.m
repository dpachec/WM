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














%%