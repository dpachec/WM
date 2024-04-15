%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 1 - SLIDE 22 - PFC ENC-MAINT WITHIN AND BETWEEN CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate clusters schematics 
%% PFC MAINT
load ([paths.results.clusters 'clustinfo_PFC_px2']);

h = zeros(52, 40); 
h(clustinfo.PixelIdxList{2}) = 1; 
figure()
imagesc(flipud(h)); hold on; 
plot([5 5],get(gca,'ylim'), 'w:','lineWidth', 2);
plot([5+5 5+5],get(gca,'ylim'), 'w:','lineWidth', 2);
set(gca, 'xTick', [], 'yTick', [])
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% PFC ENCODING
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);

h = zeros(52, 15); 
h(clustinfo.PixelIdxList{8}) = 1; 
figure(); set(gcf, 'Position', [100 100 200 400])
imagesc(flipud(h)); hold on; 
plot([5 5],get(gca,'ylim'), 'w:','lineWidth', 2);
%plot([5+8 5+8],get(gca,'ylim'), 'w:','lineWidth', 2);
set(gca, 'xTick', [], 'yTick', [])
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%% LOAD MEAN RDMS IN PFC DURING ENCODING AND MAINTENANCE (HAS TO BE ORDERED)
clear, clc

% % % PFC ENCODING
f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        %meanRDM_PFC_E{subji, 1} = rdmS; 
        %meanRDM_PFC_E{subji, 2} = oneListIds(:, 3); 

        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_PFC_E(subji,:,:) = rdm; 
        
    end
end



% % % PFC MAINTENANCE
f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        %meanRDM_PFC_M{subji, 1} = rdmS; 
        %meanRDM_PFC_M{subji, 2} = oneListIds(:, 3); 

        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_PFC_M(subji,:,:) = rdm; 
        
    end
end


clearvars -except meanRDM_PFC_E meanRDM_PFC_M 

clearvars -except meanRDM_PFC_E meanRDM_PFC_M 


%% Correlate each region during encoding and maintenance


for subji = 1:16

    mPFC_E = squeeze(meanRDM_PFC_E(subji, :, :)); 
    mPFC_M = squeeze(meanRDM_PFC_M(subji, :, :)); 

    mE = vectorizeRDM_WM(mPFC_E); 
    mM = vectorizeRDM_WM(mPFC_M); 

    isnE = isnan(mE); 
    isnM = isnan(mM); 

    mE(isnE|isnM) = []; 
    mM(isnE|isnM) = []; 

    rhoPFC(subji, :) = corr(mE, mM, 'type', 's'); 


end

sub2excPFC = [1]; 

rhoPFC1 = rhoPFC; 
rhoPFC1(sub2excPFC) = []; 

[h p ci ts] = ttest(atanh(rhoPFC1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);



%% plot one bar
clear data
data.data = atanh(rhoPFC1); 

figure(); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 50, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Mean and variance of item correlations in the two regions during the two time-periods
% FIRST LOAD WITHOUT ORDERING
clear, clc

% % % PFC ENCODING
f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM_PFC_E{subji, 1} = rdmS; 
        meanRDM_PFC_E{subji, 2} = oneListIds(:, 3); 

        
    end
end



% % % PFC MAINTENANCE
f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM_PFC_M{subji, 1} = rdmS; 
        meanRDM_PFC_M{subji, 2} = oneListIds(:, 3); 

        
    end
end


clearvars -except meanRDM_PFC_E meanRDM_PFC_M 



%% COMPUTE METRICS FOR ALL 4 TIME PERIODS

allPeriods = {meanRDM_PFC_E; meanRDM_PFC_M};

for periodi = 1:length(allPeriods)

    meanRDM = allPeriods{periodi}; 
    
    for subji = 1:size(meanRDM, 1)
        mRDM = squeeze(meanRDM{subji,1}); 
        ind = floor(meanRDM{subji, 2}/100);     
        ids = meanRDM{subji, 2}; 
        CM = generate_catModel_WM(ids); 
    
        CCI = compute_CCI_WM(mRDM, ind); 
        
        m1_CCI{periodi}(subji, :) = CCI; 
        m2_WCSM{periodi}(subji, :) = mean(averageCAT_WM(mRDM, ind));
        m3_WCSV{periodi}(subji, :) = var(averageCAT_WM(mRDM, ind));
        m4_WISM{periodi}(subji, :) = mean(vectorizeRDM(mRDM), 'omitnan'); 
        m5_WISV{periodi}(subji, :) = var(vectorizeRDM_WM(mRDM), 'omitnan'); 
        m6_WISMWC{periodi}(subji, :) = mean(mRDM(CM==1), 'omitnan'); 
        m7_WISMBC{periodi}(subji, :) = mean(mRDM(CM==0), 'omitnan'); 
        m8_WISVWC{periodi}(subji, :) = var(mRDM(CM==1), 'omitnan'); 
        m9_WISVBC{periodi}(subji, :) = var(mRDM(CM==0), 'omitnan');      
        

    end
end



%% plot two bar
clear data

%data.data = [m4_WISM{1} m4_WISM{2}]; % mean PFC
%data.data = [m5_WISV{1} m5_WISV{2}]; % VAR PFC
%data.data = [m6_WISMWC{1}, m6_WISMWC{2}]; 
%data.data = [m7_WISMBC{1}, m7_WISMBC{2}]; 
%data.data = [m8_WISVWC{1}, m8_WISVWC{2}]; 
data.data = [m9_WISVBC{1} m9_WISVBC{2}]; ylim = [-.0 .05];  % VAR PFC

%data.data = [m5_WISV_E m5_WISV_M];
%data.data = [m6_WISMWC_E, m6_WISMWC_M]; 
%data.data = [m7_WISMBC_E, m7_WISMBC_M]; 
%data.data = [m8_WISVWC_E, m8_WISVWC_M]; 
%data.data = [m9_WISVBC_E, m9_WISVBC_M]; 


sub2exc = 1; 
data.data(sub2exc, :) = []; 
figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2], data.data, 50, 'k'); hold on;
hb = plot(data.data'); hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);box on; 
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2); set(gca, 'LineWidth', 2);
set(gca, 'ylim', ylim); 

[h p ci ts] = ttest (data.data(:,1), data.data(:,2));
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)



exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Correlate with BL-NET after excluding within and between category correlations
% FIRST LOAD THE RDM in the PFC cluster during MAINTENANCE
clear, clc

f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end

%% OR LOAD THE RDM in the PFC cluster during ENCODING

clear, clc

f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end



%%

clear rhoALL 
subj_ch_fr = 7;

witBet = 0; 
AlexOrBLNET = 'Alex'; 

for subji = 1:16

    if strcmp(AlexOrBLNET, 'Alex')
        f2sav = 'Alex_pfc_E123_[7]_3-54_0_0_0_1_.1_5_1.mat'; 
        cfg = getParams(f2sav);
        paths = load_paths_WM(cfg.brainROI, cfg.net2load);
        [ACT] = load_alex_activ(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    elseif strcmp(AlexOrBLNET, 'BLNET')
        f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
        cfg = getParams(f2sav);
        paths = load_paths_WM(cfg.brainROI, cfg.net2load);
        [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    end

    ACT = squeeze(ACT); 

    ids = meanRDM{subji, 2}; 
    ids2 = char(string(ids)); 
    clear ids3
    for i = 1:length(ids2)
        idT = ids2(i,:); 
        ids3(i,:) = double(string(idT([1 3])));
    end
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);

    ACTH = ACT(ids4, ids4); 
    CM = kron(eye(6), ones(10));
    CM = CM(ids4, ids4); 
    CM(CM==0) = 2000; 
    CM = tril(CM, -1);
    CM(CM==0) = nan; 
    CM(CM==2000) = 0; 
    
    BLDist = ACTH(CM==witBet);
    %BLDist = vectorizeRDM_WM(ACTH);
   
    ACT = squeeze(meanRDM{subji, 1}); 
    neuDist = ACT(CM==witBet); 
    %neuDist = vectorizeRDM_WM(ACT); 

    rhoALL(subji, :) = corr(neuDist, BLDist, 'type', 's'); 
    
end

sub2exc = 1; 
rhoALL(sub2exc) = []; 
[h p ci ts] = ttest(rhoALL);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


%% plot one bar
clear data
data.data = [rhoALL]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 24, 'linew',1, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

set(gca, 'LineWidth', 2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% 6x6 RDM analysis


%% Correlate with BL-NET after excluding within and between category correlations
% FIRST LOAD THE RDM in the PFC cluster during MAINTENANCE
clear, clc

f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end

%% OR LOAD THE RDM in the PFC cluster during ENCODING

clear, clc

f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end


%% Compute correlation neural Distances and BL-NET distances
nPerm = 200; 
tril2u = 0; 
clear rhoALL rhoPerm rhoALLCAT neuralCatDist catDistBLNET rhoCATPerm catDistBLNETPERM rhoALLCATCM


subj_ch_fr = 7;
for subji = 2:16
    ids = meanRDM{subji, 2}; 
    
    ACT = squeeze(meanRDM{subji, 1}); 
    neuralDistances = vectorizeRDM(ACT); 
    ind = floor(ids/100); 
    ACT(logical(eye(size(ACT, 1)))) = nan; 
    neuralCatDistB = averageCAT_WM(ACT, ind); 
    
    neuralCatDist = neuralCatDistB(tril(true(length(neuralCatDistB)), tril2u));
   
    
    f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
    %f2sav = 'Alex_pfc_E123_[8]_3-54_0_0_0_1_.1_5_1.mat'; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    %[ACT] = load_alex_activ(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    ACT = squeeze(ACT); 

    ids = meanRDM{subji, 2}; 
    ids2 = char(string(ids)); 
    clear ids3
    for i = 1:length(ids2)
        idT = ids2(i,:); 
        ids3(i,:) = double(string(idT([1 3])));
    end
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);

    ACTH = ACT(ids4, ids4); 
    BLNETDistances = vectorizeRDM(ACTH); 
    BLNETDistancesB = averageCAT_WM(ACTH, ind); 
    catDistBLNET = BLNETDistancesB(tril(true(length(BLNETDistancesB)), tril2u));

    catDistCM = kron(eye(6), ones(1)); 
    catDistCM = catDistCM(tril(true(length(BLNETDistancesB)), tril2u));
    [rhoALLCATCM(subji,:)] = corr(neuralCatDist, catDistCM, 'type', 's'); 
    

    [rhoALL(subji,:)] = corr(neuralDistances', BLNETDistances', 'type', 's'); 
    [rhoALLCAT(subji,:)] = corr(neuralCatDist, catDistBLNET, 'type', 's'); 
    
    
    for permi = 1:nPerm
        
        ACTH = ACT(ids4, ids4); 
        idsperm = randperm(size(ACTH, 1)); 
        ACTH = ACTH(idsperm, idsperm); 
        BLNETDistancesPerm = vectorizeRDM(ACTH); 
        rhoPerm(permi, subji, :) = corr(neuralDistances', BLNETDistancesPerm', 'type', 's'); 

        BLNETDistancesBPERM = averageCAT_WM(ACTH, ind); 
        catDistBLNETPERM = BLNETDistancesBPERM(tril(true(length(BLNETDistancesBPERM)), tril2u));
        rhoCATPerm(permi, subji, :) = corr(neuralCatDist, catDistBLNETPERM, 'type', 's'); 

    end

end

sub2exc = 1; 
rhoPerm(:, sub2exc, :) = []; 
rhoALL(sub2exc, :) = []; 
rhoCATPerm(:, sub2exc, :) = []; 
rhoALLCAT(sub2exc, :) = []; 
rhoALLCATCM(sub2exc, :) = []; 


%% plot one subject one row
cols = repelem([0.5, 0.5, 0.5], nPerm, 1); 
figure()
scatter(rhoPerm, 1:15, 5, cols, 'o'); hold on;
scatter(rhoALL, 1:15, 500, 'r.');

clear rankInEverySubject
for subji = 1:15
    allPS = rhoPerm(:, subji); 
    allAB = length(allPS(allPS>rhoALL(subji))); 
    rankInEverySubject(subji, :) = (1-((allAB+1)/nPerm))*100
end
set(gca, 'ylim', [0 16])
text(repelem(0.1, 15), 1:15, num2str(rankInEverySubject, '%.1f'))
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot one subject one row CATEGORY
cols = repelem([0.5, 0.5, 0.5], nPerm, 1); 
figure()
scatter(rhoCATPerm, 1:15, 5, cols, 'o'); hold on;
scatter(rhoALLCAT, 1:15, 500, 'r.');
clear rankInEverySubject
for subji = 1:15
    allPS = rhoCATPerm(:, subji); 
    allAB = length(allPS(allPS>rhoALLCAT(subji))); 
    rankInEverySubject(subji, :) = (1-((allAB+1)/nPerm))*100
end
set(gca, 'ylim', [0 16])
text(repelem(0.8, 15), 1:15, num2str(rankInEverySubject, '%.1f'))
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

[h p ci ts] = ttest(rhoALLCAT)



%% plot average in six categories for encoding and maintenance
% FIRST LOAD MEAN RDMS DURING ENCODING AND MAINTENANCE (HAS TO BE ORDERED)
clear, clc

% % % PFC ENCODING
f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        %meanRDM_PFC_E{subji, 1} = rdmS; 
        %meanRDM_PFC_E{subji, 2} = oneListIds(:, 3); 

        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_PFC_E(subji,:,:) = rdm; 
        
    end
end



% % % PFC MAINTENANCE
f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        %meanRDM_PFC_M{subji, 1} = rdmS; 
        %meanRDM_PFC_M{subji, 2} = oneListIds(:, 3); 

        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_PFC_M(subji,:,:) = rdm; 
        
    end
end


clearvars -except meanRDM_PFC_E meanRDM_PFC_M 



%%
clc
ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
for subji = 1:16
    mPFC_E(subji, :, :) = averageCAT_WM(squeeze(meanRDM_PFC_E(subji, :, :)), ind);
    mPFC_M(subji, :, :) = averageCAT_WM(squeeze(meanRDM_PFC_M(subji, :, :)), ind);
end

%% compute BL-NET 6*6
ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
subj_ch_fr = 7; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
[ACT] = squeeze(load_BLNETi(cfg, subji, subj_ch_fr, paths));
ACT = averageCAT_WM(ACT, ind); 

%% Plot 3 RDMS (Encoding, Maintenance, BLNET)

d2p1 = squeeze(mean(mPFC_E)); 
d2p2 = squeeze(mean(mPFC_M)); 

tiledlayout(1, 3)
nexttile
imagesc(d2p1); axis square; colorbar
set(gca, 'clim', [-.02 0]); 
nexttile
imagesc(d2p2); axis square; colorbar
set(gca, 'clim', [-.02 0]); 
nexttile
imagesc(ACT); axis square; colorbar


exportgraphics(gcf, 'allM.png', 'Resolution', 300);


















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 2 - SLIDE 20 - HIGHER ORDER REPRESENTATIONS IN PFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First load behavioral data
clear

[all_events allTrlInfo ] = loadLogsWM;
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = allTrlInfo(~cellfun('isempty',allTrlInfo));
allTrlInfo = allTrlInfo';




%% correct item-cat and by pos

clearvars -except allTrlInfo
tic 

for vp=1:length(allTrlInfo)
    
    trialInfo = allTrlInfo{vp};
    
    cues = trialInfo(:,10);
    
    indall = find(cues == 4);
    ind123 = find(cues ~= 4);
    
    %all
    corr_all(1)=sum(trialInfo(indall,40))./numel(indall); %item
    corr_all(2)=sum(trialInfo(indall,41))./numel(indall); %category
    gData.corr_all(vp, :) = corr_all;
    
    %123
    corr_123(1)=sum(trialInfo(ind123,40))./numel(ind123); %item
    corr_123(2)=sum(trialInfo(ind123,41))./numel(ind123); %category
    gData.corr_123(vp, :) = corr_123;
    
    %all trials in each position 
    corr_all_pos(1,:)=nansum(trialInfo(indall,34:36))./numel(indall);
    corr_all_pos(2,:)=nansum(trialInfo(indall,37:39))./numel(indall);
    gData.corr_all_pos{vp} = corr_all_pos;
    %123 trials in each position 
    corr_123_pos(1,:)=nansum(trialInfo(ind123,34:36))./(numel(ind123)./3);
    corr_123_pos(2,:)=nansum(trialInfo(ind123,37:39))./(numel(ind123)./3);
    %den2use = sum(~isnan(trialInfo(ind123,37:39))); 
    %corr_123_pos(1,:)=nansum(trialInfo(ind123,34:36))./den2use;
    %corr_123_pos(2,:)=nansum(trialInfo(ind123,37:39))./den2use;
    gData.corr_123_pos{vp} = corr_123_pos;
    
end

toc


%% average subjects not sessions

corr_all = gData.corr_all;
corr_123 = gData.corr_123;
corr_all_pos = cell2mat (gData.corr_all_pos); 
corr_all_pos_it  = corr_all_pos(1,:); corr_all_pos_it1 = reshape(corr_all_pos_it', 3, [])';
corr_all_pos_cat = corr_all_pos(2,:); corr_all_pos_cat1 = reshape(corr_all_pos_cat', 3, [])';
corr_all_pos1(:, 1:3) = corr_all_pos_it1; 
corr_all_pos1(:, 4:6) = corr_all_pos_cat1;
corr_123_pos = cell2mat (gData.corr_123_pos); 
corr_123_pos_it  = corr_123_pos(1,:); corr_123_pos_it1 = reshape(corr_123_pos_it', 3, [])';
corr_123_pos_cat = corr_123_pos(2,:); corr_123_pos_cat1 = reshape(corr_123_pos_cat', 3, [])';
corr_123_pos1(:, 1:3) = corr_123_pos_it1; 
corr_123_pos1(:, 4:6) = corr_123_pos_cat1;


region = 'all';
sub2exc = [];
corr_all = average_xGM (corr_all, region); 
corr_123 = average_xGM (corr_123, region); 
corr_all_pos = average_xGM (corr_all_pos1, region); 
corr_123_pos = average_xGM (corr_123_pos1, region); 

m_corr_all = mean(corr_all) ; std_corr_all = std(corr_all); se_corr_all = std_corr_all / sqrt ( size(corr_all, 1));
m_corr_123 = mean(corr_123) ; std_corr_123 = std(corr_123); se_corr_123 = std_corr_123 / sqrt ( size(corr_123, 1));
m_corr_all_pos = mean(corr_all_pos) ; std_corr_all_pos = std(corr_all_pos); se_corr_all_pos = std_corr_all_pos / sqrt ( size(corr_all_pos, 1));
m_corr_123_pos = mean(corr_123_pos) ; std_corr_123_pos = std(corr_123_pos); se_corr_123_pos = std_corr_123_pos / sqrt ( size(corr_123_pos, 1));

allPerf = [corr_123 corr_all];
allPerfPos = [corr_123_pos corr_all_pos];




%% 3 Bar POSITIONS
clear data
meanPF= cat(3, corr_123_pos(:, 1:3), corr_all_pos(:, 1:3))
meanPF1 = squeeze(mean(meanPF, 3)); 
data.data = [meanPF1]; 


figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data);


h = bar (mean_S,'FaceColor','flat', 'linew', 2);hold on;
h.CData = [1 1 1; 1 1 1; 1 1 1]; 
hb = plot ([1:3], data.data, 'k'); hold on; 
set(hb, 'lineWidth', 1, 'Marker', 'o', 'MarkerSize',4);hold on;
set(hb,'linestyle','none', 'linew', 2);
set(gca,'XTick',[1 3],'XTickLabel',{'', ''}, 'FontSize', 22, 'xlim', [0.25 3.75], 'ylim', [.4 1.1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);hold on;
set(gca, 'LineWidth', 1);
box on

[h p ci ts] = ttest(data.data(:, 1), data.data(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(data.data(:, 1), data.data(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(data.data(:, 2), data.data(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

exportgraphics(gcf, 'myP.png','Resolution', '150');

%% 2 Bar SINGLE MULTIPLE
clear data
meanPF= cat(3, corr_123_pos(:, 1:3), corr_all_pos(:, 1:3))
meanPF1 = squeeze(mean(meanPF, 2)); 
data.data = [meanPF1]; 


figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data);


h = bar (mean_S,'FaceColor','flat', 'linew', 2);hold on;
h.CData = [1 1 1; 1 1 1]; 
hb = plot ([1:2], data.data, 'k'); hold on; 
set(hb, 'lineWidth', 1, 'Marker', 'o', 'MarkerSize',4);hold on;
%set(hb,'linestyle','none', 'linew', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 22, 'xlim', [0.25 2.75], 'ylim', [.4 1.1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);hold on;
set(gca, 'LineWidth', 1);
box on

[h p ci ts] = ttest(data.data(:, 1), data.data(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

exportgraphics(gcf, 'myP.png','Resolution', '150');


%% 6 Bar 
clear data
data.data = [corr_123_pos(:, 1:3) corr_all_pos(:, 1:3)]; 


figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data);


h = bar (mean_S,'FaceColor','flat', 'linew', 2);hold on;
h.CData = [1 1 1; 1 1 1; 1 1 1; .5 .5 .5; .5 .5 .5; .5 .5 .5]; 
hb = plot ([1:6], data.data, 'k'); hold on; 
set(hb, 'lineWidth', 1, 'Marker', 'o', 'MarkerSize',4);hold on;
set(hb,'linestyle','none', 'linew', 2);
set(gca,'XTick',[1 6],'XTickLabel',{'', ''}, 'FontSize', 22, 'xlim', [0.25 6.75], 'ylim', [.3 1.3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);hold on;
set(gca, 'LineWidth', 1);
box on

exportgraphics(gcf, 'myP.png','Resolution', '150');



%% LME
clc
allD = data.data(:); 
positions = [[ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]' ; [ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]'];
singORMult = [ones(1, 96) ones(1, 96)*2]';
subID = [1:32 1:32 1:32 1:32 1:32 1:32]';
d4LME = [allD positions singORMult subID];

tbl = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), ...
    'VariableNames',{'performance','position','trial_type', 'subID'});
lme = fitlme(tbl,'performance ~ position + trial_type + (1|subID)'); % random intercept model

lme

%% 
clc

[h p ci ts] = ttest(corr_123_pos(:, 1), corr_123_pos(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_123_pos(:, 1), corr_123_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_123_pos(:, 2), corr_123_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), corr_all_pos(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), corr_all_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 2), corr_all_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


[h p ci ts] = ttest(corr_123_pos(:, 1), mean([corr_123_pos(:, 2) corr_123_pos(:, 3)], 2) );
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), mean([corr_all_pos(:, 2) corr_all_pos(:, 3)], 2) );
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

singleItems = mean([corr_123_pos(:, 1) corr_123_pos(:, 2) corr_123_pos(:, 3)], 2) ; 
multiItems = mean([corr_all_pos(:, 1) corr_all_pos(:, 2) corr_all_pos(:, 3)], 2) ; 
[h p ci ts] = ttest(singleItems, multiItems);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);



%% ANOVA VERSION
clc
d4ANOVA = [data.data(:, 1) data.data(:, 2) data.data(:, 3) data.data(:, 4) data.data(:, 5) data.data(:, 6)]; 
T = array2table(d4ANOVA);
T.Properties.VariableNames = {'Pos1S' 'Pos2S' 'Pos3S' 'Pos1M' 'Pos2M' 'Pos3M'};
% create the within-subjects design
withinDesign = table([1 2 3 1 2 3]',[1 1 1 2 2 2]', 'VariableNames',{'Positions', 'SM'});
withinDesign.Positions = categorical(withinDesign.Positions);
withinDesign.SM = categorical(withinDesign.SM);
% create the repeated measures model and do the anova
rm = fitrm(T,'Pos1S-Pos3M ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Positions*SM'); % remove comma to see ranova's table
tbl = multcompare(rm, 'Positions', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
%tbl = multcompare(rm, 'SM', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
%tbl = multcompare(rm, 'Positions', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));



%%

f = 'https://www.mathworks.com/matlabcentral/answers/uploaded_files/696039/Q2.xls';
T = readtable(f);
% ingore extra data at end
T = T(1:2*2*17,:);
% organize the data in a 17x4 matrix (1 row per participant)
dv = array2table(reshape(T.Value, 17, [])); % dependent variable is 'Value'
% codes for column names: gr=grouped, ug=ungrouped, b1=beta1, b2=beta2
dv.Properties.VariableNames = {'grb1', 'grb2', 'ugb1', 'ugb2' };
% create the within-subjects design
withinDesign = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Crowding','Beta'});
withinDesign.Crowding = categorical(withinDesign.Crowding);
withinDesign.Beta = categorical(withinDesign.Beta);
% create repeated measures model
rm = fitrm(dv, 'grb1-ugb2 ~ 1', 'WithinDesign', withinDesign);
% perform anova (remove semicolon to view table generated by ranova)
AT = ranova(rm, 'WithinModel', 'Crowding*Beta');
% output an improved version of the anova table
disp(anovaTable(AT, 'Value'));








%% LOAD all conditions
clearvars

region = 'pfc'; 
paths = load_paths_WM(region, 'none');

contrasts = {
              
              'DISC_EE1' 'DIDC_EE1';
              'DISC_EE2' 'DIDC_EE2';
              'DISC_EE3' 'DIDC_EE3';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    if exist("all_IDs")
        idData{i,:} = all_IDs;
    end
end



%% ANOVA EE1 EE2 EE3
clc
%clear all; 

%sub2exc = [1]; 
sub2exc = [18 22]; 

%define conditions 
cond1 = 'DISC_EE1';
cond2 = 'DISC_EE2';
cond3 = 'DISC_EE3';
cond1X = 'DIDC_EE1';
cond2X = 'DIDC_EE2';
cond3X = 'DIDC_EE3';

cond1B = eval(cond1);
cond2B = eval(cond2);
cond3B = eval(cond3);
cond1BX = eval(cond1X);
cond2BX = eval(cond2X);
cond3BX = eval(cond3X);

cond1C = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond1B, 'un', 0)));
cond2C = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond2B, 'un', 0)));
cond3C = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond3B, 'un', 0)));
cond1CX = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond1BX, 'un', 0)));
cond2CX = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond2BX, 'un', 0)));
cond3CX = double(string(cellfun(@(x) mean(x(:, 6:13, 6:13), 'all', 'omitnan'), cond3BX, 'un', 0)));

cond1C = cond1C - cond1CX; 
cond2C = cond2C - cond2CX; 
cond3C = cond3C - cond3CX; 

cond1C(sub2exc) = []; 
cond2C(sub2exc) = []; 
cond3C(sub2exc) = []; 


%% 
nSub = size(cond1C, 1); 
d4ANOVA = [cond1C; cond2C ;cond3C]; 
d4ANOVA(:,2) = [ones(1,nSub) ones(1,nSub)*2 ones(1,nSub)*3];
d4ANOVA(:,3) = [1:nSub 1:nSub 1:nSub];

[p F] = RMAOV1(d4ANOVA);


%% ANOVA COMPARING EMS at each position (SAME ANAYLSIS FOR ENCODING IN THE TWO REGOINS)

%first load data

% VVS
%cd D:\_WM\analysis\pattern_similarity\vvs\100ms\EM21_EM22_EM23\3-8Hz

clear 

tic
region = 'pfc'; 
paths = load_paths_WM(region, 'none');

contrasts = {
              
              'DISC_EE1' 'DIDC_EE1';
              'DISC_EE2' 'DIDC_EE2';
              'DISC_EE3' 'DIDC_EE3';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    if exist("all_IDs")
        idData{i,:} = all_IDs;
    end
end

toc
%% ANOVA EE1 EE2 EE3
clc
%clear all; 

%sub2exc = [1]; 
%sub2exc = [18 22]; 
sub2exc = [1]; 

%define conditions 
cond1 = 'DISC_EE1';
cond2 = 'DISC_EE2';
cond3 = 'DISC_EE3';
cond1X = 'DIDC_EE1';
cond2X = 'DIDC_EE2';
cond3X = 'DIDC_EE3';

cond1B = eval(cond1);
cond2B = eval(cond2);
cond3B = eval(cond3);
cond1BX = eval(cond1X);
cond2BX = eval(cond2X);
cond3BX = eval(cond3X);

tP1 = 6:13; 
tP2 = 6:13; 
cond1C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond1B, 'un', 0)));
cond2C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond2B, 'un', 0)));
cond3C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond3B, 'un', 0)));
cond1CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond1BX, 'un', 0)));
cond2CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond2BX, 'un', 0)));
cond3CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond3BX, 'un', 0)));

cond1C = cond1C - cond1CX; 
cond2C = cond2C - cond2CX; 
cond3C = cond3C - cond3CX; 

cond1C(sub2exc) = []; 
cond2C(sub2exc) = []; 
cond3C(sub2exc) = []; 


%% 
nSub = size(cond1C, 1); 
d4ANOVA = [cond1C; cond2C ;cond3C]; 
d4ANOVA(:,2) = [ones(1,nSub) ones(1,nSub)*2 ones(1,nSub)*3];
d4ANOVA(:,3) = [1:nSub 1:nSub 1:nSub];

[p F] = RMAOV1(d4ANOVA);


%% WITHIN SUBJECTS DESING IN MATLAB

d4ANOVA = [[1:length(cond1C)]' cond1C cond2C cond3C]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'Pos1' 'Pos2' 'Pos3'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Model'});
withinDesign.Model = categorical(withinDesign.Model);
% create the repeated measures model and do the anova
rm = fitrm(T,'Pos1-Pos3~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));



%% ANOVA COMPARING EES at each position (SAME ANAYLSIS FOR PRIORITIZATION IN THE TWO REGOINS)

%first load data

% VVS
%cd D:\_WM\analysis\pattern_similarity\pfc\100ms\EM21_EM22_EM23\3-8Hz

clear 

tic
region = 'pfc'; 
paths = load_paths_WM(region, 'none');

contrasts = {
              
              'DISC_EM21' 'DIDC_EM21';
              'DISC_EM22' 'DIDC_EM22';
              'DISC_EM23' 'DIDC_EM23';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    if exist("all_IDs")
        idData{i,:} = all_IDs;
    end
end

toc
%% ANOVA EM21 EM22 EM23
clc
%clear all; 

%sub2exc = [1]; 
%sub2exc = [18 22]; 
%sub2exc = [18]; 
sub2exc = []; 

%define conditions 
cond1 = 'DISC_EM21';
cond2 = 'DISC_EM22';
cond3 = 'DISC_EM23';
cond1X = 'DIDC_EM21';
cond2X = 'DIDC_EM22';
cond3X = 'DIDC_EM23';

cond1B = eval(cond1);
cond2B = eval(cond2);
cond3B = eval(cond3);
cond1BX = eval(cond1X);
cond2BX = eval(cond2X);
cond3BX = eval(cond3X);

tP1 = 6:10; 
tP2 = 6:10; 
cond1C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond1B, 'un', 0)));
cond2C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond2B, 'un', 0)));
cond3C = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond3B, 'un', 0)));
cond1CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond1BX, 'un', 0)));
cond2CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond2BX, 'un', 0)));
cond3CX = double(string(cellfun(@(x) mean(x(:, tP1, tP2), 'all', 'omitnan'), cond3BX, 'un', 0)));

cond1C = cond1C - cond1CX; 
cond2C = cond2C - cond2CX; 
cond3C = cond3C - cond3CX; 

cond1C(sub2exc) = []; 
cond2C(sub2exc) = []; 
cond3C(sub2exc) = []; 


%% 
nSub = size(cond1C, 1); 
d4ANOVA = [cond1C; cond2C ;cond3C]; 
d4ANOVA(:,2) = [ones(1,nSub) ones(1,nSub)*2 ones(1,nSub)*3];
d4ANOVA(:,3) = [1:nSub 1:nSub 1:nSub];

[p F] = RMAOV1(d4ANOVA);


%% WITHIN SUBJECTS DESING IN MATLAB

d4ANOVA = [[1:length(cond1C)]' cond1C cond2C cond3C]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'Pos1' 'Pos2' 'Pos3'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Model'});
withinDesign.Model = categorical(withinDesign.Model);
% create the repeated measures model and do the anova
rm = fitrm(T,'Pos1-Pos3~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 6 - SLIDE 21 - NOISE AND SIMULATION ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Process and plot mean RDM in the PFC cluster (MAINTENANCE) 
clear, clc

f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end

%% Process and plot mean RDM in the PFC cluster (ENCODING)

clear, clc

f2load = 'pfc_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   
load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; % % USE HERE 7 or 8 for theta or beta cluster

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end


%% Plot category model under different levels of noise

noiseLevels = [0 .1 .5 1 3 5 10]

for noisei = 1:length(noiseLevels)

    mu = noiseLevels(noisei); 
    ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
    
    CM = kron(eye(6), ones(10));
    CM = CM + normrnd(0, mu,[60,60]);
    CM = normalize(CM, 'range', [-1, 1]);
    %CM = mat2gray(CM);
    CM(logical(eye(size(CM)))) = 1; 
    CM = tril(CM,-1)+tril(CM)'; 
    
    figure()
    tiledlayout(1, 2)
    nexttile
    imagesc(CM); ; axis square ;
    set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])

    d2p = 1- CM;
    [rdmMDS] = cmdscale(d2p);
    %Encoding
            
    cols = brewermap(6, 'Accent'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
    cols = repelem(cols, 10, 1);
    nexttile
    scatter(rdmMDS(:,1),rdmMDS(:,2),550, cols, '.'); axis square
    set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
    box on
    
    
    %sgtitle(['Noise: ' num2str(mu) ' CCI: ' num2str(CCI)])

    exportgraphics(gcf, ['allM' num2str(noisei, '%02.f') '.png'], 'Resolution', 300);
    close all


end


%% COUNT significant correlations for within and between categories
nPerm = 1000; 
noiseLevels = [0 .1 .5 1 3 5 10]

clearvars -except meanRDM noiseLevels nPerm
for subji = 1:16
    ids = meanRDM{subji, 2}; 
    ids2 = char(string(ids)); 
    clear ids3
    for i = 1:length(ids2)
        idT = ids2(i,:); 
        ids3(i,:) = double(string(idT([1 3])));
    end
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
    ind = floor(meanRDM{subji, 2}/100); 
    nTrials = length(ind); 

    f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    subj_ch_fr = subji; 
    [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    ACT = squeeze(ACT); 
    realDataBLNET = ACT(ids4, ids4); 
    realDataBLNET = vectorizeRDM_WM(realDataBLNET); 

     for noisei = 1:length(noiseLevels)
        mu = noiseLevels(noisei); 

        realDataSubj = meanRDM{subji, 1}; 
        realDataSubj = vectorizeRDM_WM(realDataSubj); 
    
        neuralBLNETCORR(subji, :) = corr(realDataBLNET, realDataSubj, 'type', 's'); 
    
    
        for permi = 1:nPerm
            CM = kron(eye(6), ones(10));
            CM = CM + normrnd(0, mu,[60,60]);
            CM = normalize(CM, 'range', [-1 1]);
            %CM(logical(eye(size(CM)))) = nan; 
            %CM = tril(CM,-1)+tril(CM)'; 
            
            CM = CM(ids4, ids4); 
            catModNoise = vectorizeRDM_WM(CM); 
        
            [rhoX pX] = corr(realDataSubj, catModNoise, 'type', 's');  
            rhoPermALL(subji, permi, noisei, :) = rhoX; 
            %pPermALL(subji, permi, noisei, :) = pX; 
            rhoPermALLBLNET(subji, permi,noisei, :) = corr(realDataBLNET, catModNoise, 'type', 's'); 
    
        end
     end


end

sub2exc = [1]; 

neuralBLNETCORR(sub2exc) = [];
rhoPermALL(sub2exc,:,:) = [];
rhoPermALLBLNET(sub2exc,:,:) = [];

[h p ci ts] = ttest(neuralBLNETCORR); 
t4 = ts.tstat; 

[hPerm pPerm ci ts] = ttest(rhoPermALL); 
tPermALL = squeeze(ts.tstat); pPermALL = squeeze(pPerm); 

[hPerm pPerm ci ts] = ttest(rhoPermALLBLNET); 
tPermALLBLNET = squeeze(ts.tstat); 

%% count % of significant p's

clear pH pH1
for noisei = 1:7
    pH = pPermALL(:,  noisei);
    pH1(:, noisei) = sum(pH<0.05, 2);
end

%% plot histograms for analysis 1 and 2

figure(); set(gcf, 'Position', [100 100 1000 1000])
tiledlayout(7, 2)
for tilei = 1:7
    if tilei == 1
        nexttile
        scatter(tPermALL(1, 1), 0, 550, '.'); hold on; 
        scatter(t4, 0, 550, '.', 'r'); box on; 
        %set(gca, 'xlim', [-5 7])
        nexttile
        scatter(tPermALLBLNET(1,1), 0, 550, '.'); hold on; 
        scatter(t4, 0, 550, '.', 'r'); box on; 
        %set(gca, 'xlim', [-5 7])
    else
        nexttile
        histogram(tPermALL(:, tilei), 10); hold on
        scatter(t4, 0, 550, '.', 'r'); box on; 
        %set(gca, 'xlim', [-5 7])
        nexttile
        histogram(tPermALLBLNET(:, tilei), 10); hold on
        scatter(t4, 0, 550, '.', 'r'); box on; 
        %set(gca, 'xlim', [-5 7])
    end
    
end
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% ANALYSIS OF VARIANCE across subjects for EACH PAIRWISE COMPARISON
%% First load mean RDM in the PFC cluster (MAINTENANCE) 
clear, clc

f2load = 'pfc_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]); 
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 


t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{2}, 3); 
nTimes = size(allNeuralRDMS{2}, 4); 

for subji = 1:nSubjs
    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
        oneListIds = double(string(cat(1, oneListIds{:})));
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        meanRDM{subji, 1} = rdmS; 
        meanRDM{subji, 2} = oneListIds(:, 3); 
        
    end
end






%% COMPUTE VARIANCE FOR PERMUTATIONS

nPerm = 100; 
tril2u = 0; 

clear neuralDistances neuralDistancesPerm 
for subji = 2:16
    clear ACT
    ids = meanRDM{subji, 2}; 
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (ids, idsF2);
    ACT(x3, x3) = squeeze(meanRDM{subji, 1}); 
    if size(ACT, 1) < 60
        ACT(60, 1:60) = 0; 
        ACT(1:60, 60) = 0; 
    end
    ACT(ACT ==0)= nan; 
    %ACT(logical(eye(size(ACT, 1)))) = nan; 
    %ACT = mat2gray(ACT); 
    ind = repelem(1:6, 10); 
    
    
    CM = kron(eye(6), ones(10));
    CM(CM==0) = 2000; 
    CM = tril(CM, -1);
    CM(CM==0) = nan; 
    CM(CM==2000) = 0; 
    
        
    neuralDistances(subji, :) = vectorizeRDM(ACT); 
    %neuralDistances(subji, :) = ACT(CM==tril2u); 

    for permi = 1:nPerm
        
        ACTH = ACT; 
        idsperm = randperm(size(ACTH, 1)); 
        ACTH = ACTH(idsperm, idsperm); 
        neuralDistancesPerm(permi, subji,:) = vectorizeRDM(ACTH); 
        %neuralDistancesPerm(permi, subji,:) = ACTH(CM==tril2u); 
        
    end

end

sub2exc = [1]; 
%sub2exc = []; 
neuralDistances(sub2exc, :) = []; 
neuralDistancesPerm(:, sub2exc, :) = []; 




%%
nIts = size(neuralDistances, 2);
var0 = var(neuralDistances,[],  1, 'omitnan'); 
var1 = squeeze(var(neuralDistancesPerm, [], 2, 'omitnan')); 
var1 = mean(var1, 1); 
varBelowMean = sum(double(var0 < var1)); 

perc = (varBelowMean * 100 ) / nIts



%% quantify variance with respect to shuffled variance for every item

nIts = size(neuralDistances, 2); 
h = zeros(nIts, 1); 
thre2u =  nPerm*0.05; 
var1 = squeeze(var(neuralDistancesPerm, [], 2, 'omitnan')); 

for itemi = 1:nIts

    varDistr = var1(:, itemi); 
    varObs = var0(itemi); 

    allAB = sum(varObs<varDistr);

    if allAB <= thre2u
        h(itemi, :) = 1; 
    end
end

perc = (sum(h) * 100 ) / nIts


%% plot one subject one row
nIts = 20; %1770

d2p1 = var(neuralDistances,[],  1, 'omitnan'); 
d2p2 = squeeze(var(neuralDistancesPerm, [], 2, 'omitnan'));

cols = repelem([0.5, 0.5, 0.5], nPerm, 1); 
figure(); set(gcf, 'Position', [100 100 2000 200])
scatter(1:nIts, d2p2(:,1:nIts), 200, 'r.'); hold on
scatter(1:nIts, d2p1(:,1:nIts), 600, 'k.'); 
set(gca, 'xlim', [0 21], fontsize= 16)
box on
exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% PFC LATERAL ELECTRODES




























%%

























%%