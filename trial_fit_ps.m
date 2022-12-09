%% LOAD ALL IN PFC AND VVS
%% ENCODNG
clear, clc

region = 'vvs'; 
paths = load_paths_WM(region);

%contrasts = { 'DISC_EE' 'DIDC_EE'};
contrasts = { 'ALL_EE' };

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

ALLE_VVS = ALL_EE; 
ids_VVS = all_IDs; 

clearvars -except ALLE_VVS ids_VVS contrasts

region = 'pfc'; 
paths = load_paths_WM(region);


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

ALLE_PFC = ALL_EE; 
ids_PFC = all_IDs; 

ALLE_PFC = ALLE_PFC([2 3  5  9 10 11 12 14 15 16]);
ALLE_VVS = ALLE_VVS([7 9 13 18 19 20 21 23 27 28]); 
ids_PFC = ids_PFC([2 3  5  9 10 11 12 14 15 16]);
ids_VVS = ids_VVS([7 9 13 18 19 20 21 23 27 28]); 

clearvars -except ALLE_PFC ALLE_VVS  ids_VVS ids_PFC

%% analysis for one frequency band

clearvars -except ALLE_PFC ALLE_VVS  ids_VVS ids_PFC

timeL = 6:15; 

for subji = 1:10

    ids = ids_VVS{subji}; 
    c1 = ids(: ,1); 
    c1 = double(string(cellfun(@(x) x{1}(7), c1, 'un', 0)));
    c2 = ids(: ,2); 
    c2 = double(string(cellfun(@(x) x{1}(7), c2, 'un', 0)));

    for triali = 1:length(ids)
        tr2cmp = c1(triali); 
        
        % % % extract diagoanl for DISC and DIDC comparisons
        cc1VVS = ALLE_VVS{subji}(c2 == tr2cmp, timeL, timeL); 
        for i = 1:size(cc1VVS, 1)
            ccc1VVS(i,:) = diag(squeeze(cc1VVS(i,:,:)));
        end
        cc2VVS = ALLE_VVS{subji}(c2 ~= tr2cmp, timeL, timeL); 
        for i = 1:size(cc2VVS, 1)
            ccc2VVS(i,:) = diag(squeeze(cc2VVS(i,:,:)));
        end
        
        cc1PFC = ALLE_PFC{subji}(c2 == tr2cmp, timeL, timeL); 
        for i = 1:size(cc1PFC, 1)
            ccc1PFC(i,:) = diag(squeeze(cc1PFC(i,:,:)));
        end
        cc2PFC= ALLE_PFC{subji}(c2 ~= tr2cmp, timeL, timeL); 
        for i = 1:size(cc2PFC, 1)
            ccc2PFC(i,:) = diag(squeeze(cc2PFC(i,:,:)));
        end

        corrTR_VVS(triali,:) = mean(ccc1VVS, 'all', 'omitnan') - mean(ccc2VVS, 'all', 'omitnan') ;
        corrTR_PFC(triali,:) = mean(ccc1PFC, 'all', 'omitnan') - mean(ccc2PFC, 'all', 'omitnan') ;
        
    end


    
    rhoALL(subji, :) = corr(corrTR_PFC, corrTR_VVS, 'type', 's');


end


disp('done')

[h p ci ts] = ttest(rhoALL)




%% ANAYLSIS FOR ALL FREQUENCY BANDS

clear , clc
paths = load_paths_WM('pfc');
main_path_pfc = paths.currSessEE; 
paths = load_paths_WM('vvs');
main_path_vvs = paths.currSessEE;  


fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

for bandi = 1:length(fnames)

    flist = dir([main_path_vvs fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_vvs fnames{bandi} '/' f2load])
    ALLE_VVS_B{bandi,:}  = ALL_EE([7 9 13 18 19 20 21 23 27 28]); 
    ids_VVS_B{bandi,:}  = all_IDs([7 9 13 18 19 20 21 23 27 28]); 
    
    flist = dir([main_path_pfc fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_pfc fnames{bandi} '/' f2load])
    ALLE_PFC_B{bandi,:}  = ALL_EE([2 3  5  9 10 11 12 14 15 16]);
    ids_PFC_B{bandi,:}  = all_IDs([2 3  5  9 10 11 12 14 15 16]);
    


end

%% analysis for all frequency bands

clearvars -except ALLE_PFC_B ALLE_VVS_B ids_VVS_B ids_PFC_B

timeL = 6:15; 

tic

for bandi = 1:length(ALLE_VVS_B)

    ids_VVS = ids_VVS_B{bandi};
    ids_PFC = ids_PFC_B{bandi}; 
    ALLE_VVS = ALLE_VVS_B{bandi}; 
    ALLE_PFC = ALLE_PFC_B{bandi}; 

    for subji = 1:10
    
        ids = ids_VVS{subji}; 
        c1 = ids(: ,1); 
        c1 = double(string(cellfun(@(x) x{1}(7), c1, 'un', 0)));
        c2 = ids(: ,2); 
        c2 = double(string(cellfun(@(x) x{1}(7), c2, 'un', 0)));
    
        for triali = 1:length(ids)
            tr2cmp = c1(triali); 
            
            % % % extract diagoanl for DISC and DIDC comparisons
            cc1VVS = ALLE_VVS{subji}(c2 == tr2cmp, timeL, timeL); 
            for i = 1:size(cc1VVS, 1)
                ccc1VVS(i,:) = diag(squeeze(cc1VVS(i,:,:)));
            end
            cc2VVS = ALLE_VVS{subji}(c2 ~= tr2cmp, timeL, timeL); 
            for i = 1:size(cc2VVS, 1)
                ccc2VVS(i,:) = diag(squeeze(cc2VVS(i,:,:)));
            end
            
            cc1PFC = ALLE_PFC{subji}(c2 == tr2cmp, timeL, timeL); 
            for i = 1:size(cc1PFC, 1)
                ccc1PFC(i,:) = diag(squeeze(cc1PFC(i,:,:)));
            end
            cc2PFC= ALLE_PFC{subji}(c2 ~= tr2cmp, timeL, timeL); 
            for i = 1:size(cc2PFC, 1)
                ccc2PFC(i,:) = diag(squeeze(cc2PFC(i,:,:)));
            end
    
            corrTR_VVS(triali,:) = mean(ccc1VVS, 'all', 'omitnan') - mean(ccc2VVS, 'all', 'omitnan') ;
            corrTR_PFC(triali,:) = mean(ccc1PFC, 'all', 'omitnan') - mean(ccc2PFC, 'all', 'omitnan') ;
            
        end
    
    
        
        rhoALL(subji, :) = corr(corrTR_PFC, corrTR_VVS, 'type', 's');
    
    end


    disp('done')
    
    [h(bandi) p(bandi) ci ts] = ttest(rhoALL);
    t(bandi) = ts.tstat ;
    rhoALL_B(:, bandi) = rhoALL; 



end


toc


save ('rhoALL_B', 'rhoALL_B', 't', 'p', 'h')

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