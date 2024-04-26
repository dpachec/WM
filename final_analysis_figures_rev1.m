%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 1 comment 1 - SLIDE 8 
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

%% VVS ENCODING
load ([paths.results.clusters 'CAT_VVS_ENC_px1']);

h = zeros(52, 15); 
h(clustinfo.PixelIdxList{1}) = 1; 
figure(); set(gcf, 'Position', [100 100 200 400])
imagesc(flipud(h)); hold on; 
plot([5 5],get(gca,'ylim'), 'w:','lineWidth', 2);
set(gca, 'xTick', [], 'yTick', [])
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% VVS MAINT
load ([paths.results.clusters 'all_clustinfo_VVS']);
h = zeros(52, 40); 

%c2u = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);
c2u=intersect(intersect(allClustInfo{4}.PixelIdxList{14},allClustInfo{5}.PixelIdxList{25},'stable'), allClustInfo{6}.PixelIdxList{17},'stable')

h(allClustInfo{4}.PixelIdxList{14}) = 1; 
h(allClustInfo{5}.PixelIdxList{25}) = 2; 
h(allClustInfo{6}.PixelIdxList{17}) = 3; 
h(c2u) = 4; 


figure()
imagesc(flipud(h)); hold on; 
plot([5 5],get(gca,'ylim'), 'w:','lineWidth', 2);
plot([5+5 5+5],get(gca,'ylim'), 'w:','lineWidth', 2);
set(gca, 'xTick', [], 'yTick', [])
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%% LOAD MEAN RDMS IN THE TWO REGIONS DURING ENCODING AND MAINTENANCE (HAS TO BE ORDERED)
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

% % % VVS ENCODING
f2load = 'vvs_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'CAT_VVS_ENC_px1.mat']);
idsClust = unique([clustinfo.PixelIdxList{1}]);

for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    idsF1 = oneListIds(:, 3);
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (idsF1, idsF2);
    
    %rdmS         = squeeze(mean(mean(neuralRDMs(:, :, 1:6,12:24), 3),4)); %Theta cluster time period (checked in the out_real sigMH_real : Exact time period)
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    
    %meanRDM_VVS_E{subji, 1} = rdmS; 
    %meanRDM_VVS_E{subji, 2} = oneListIds(:, 3); 

    rdm = nan(60); 
    rdm(x3,x3) = rdmS; 
    meanRDM_VVS_E(subji,:,:) = rdm; 

    
end



% % % % VVS MAINTENANCE
clearvars -except meanRDM_PFC_E meanRDM_PFC_M meanRDM_VVS_E


f2load = 'vvs_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'all_clustinfo_VVS']);
idsClust = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);


for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    idsF1 = oneListIds(:, 3);
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (idsF1, idsF2);
    
    %rdmS         = squeeze(mean(mean(neuralRDMs(:, :, 1:6,12:24), 3),4)); %Theta cluster time period (checked in the out_real sigMH_real : Exact time period)
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    
    %meanRDM_VVS_M{subji, 1} = rdmS; 
    %meanRDM_VVS_M{subji, 2} = oneListIds(:, 3); 

    rdm = nan(60); 
    rdm(x3,x3) = rdmS; 
    meanRDM_VVS_M(subji,:,:) = rdm; 

    
end


clearvars -except meanRDM_PFC_E meanRDM_PFC_M meanRDM_VVS_E meanRDM_VVS_M

%% Correlate each region during encoding and maintenance

% % % % VVS 

for subji = 1:28
    
    mVVS_E = squeeze(meanRDM_VVS_E(subji, :, :)); 
    mVVS_M = squeeze(meanRDM_VVS_M(subji, :, :)); 

    mE = vectorizeRDM_WM(mVVS_E); 
    mM = vectorizeRDM_WM(mVVS_M); 

    isnE = isnan(mE); 
    isnM = isnan(mM); 

    mE(isnE|isnM) = []; 
    mM(isnE|isnM) = []; 

    rhoVVS(subji, :) = corr(mE, mM, 'type', 's'); 


end


% % % % PFC

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


sub2excVVS = [18 2]; 
sub2excPFC = [1]; 

rhoVVS1 = rhoVVS; 
rhoPFC1 = rhoPFC; 
rhoVVS1(sub2excVVS) = []; 
rhoPFC1(sub2excPFC) = []; 



[h p ci ts] = ttest(atanh(rhoVVS1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(atanh(rhoPFC1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

rhoVVS2= rhoVVS([7 9 13 18 19 20 21 23 27 28]);
rhoPFC2= rhoPFC([2 3  5  9 10 11 12 14 15 16]);

diff = atanh(rhoVVS2)-atanh(rhoPFC2); 
[h p ci ts] = ttest(diff);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

%% plot one bar
clear data

%data.data = atanh(rhoVVS1); 
data.data = atanh(rhoPFC1); 

figure(); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 50, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .175])
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

% % % VVS ENCODING
f2load = 'vvs_E123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'CAT_VVS_ENC_px1.mat']);
idsClust = unique([clustinfo.PixelIdxList{1}]);

for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    
    %rdmS         = squeeze(mean(mean(neuralRDMs(:, :, 1:6,12:24), 3),4)); %Theta cluster time period (checked in the out_real sigMH_real : Exact time period)
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    
    meanRDM_VVS_E{subji, 1} = rdmS; 
    meanRDM_VVS_E{subji, 2} = oneListIds(:, 3); 

    
end



% % % % VVS MAINTENANCE
clearvars -except meanRDM_PFC_E meanRDM_PFC_M meanRDM_VVS_E

f2load = 'vvs_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'all_clustinfo_VVS']);
idsClust = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);


for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    
    meanRDM_VVS_M{subji, 1} = rdmS; 
    meanRDM_VVS_M{subji, 2} = oneListIds(:, 3); 

    
    
end


clearvars -except meanRDM_PFC_E meanRDM_PFC_M meanRDM_VVS_E meanRDM_VVS_M


%% COMPUTE METRICS FOR ALL 4 TIME PERIODS

allPeriods = {meanRDM_PFC_E; meanRDM_PFC_M; meanRDM_VVS_E; meanRDM_VVS_M};

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



% data.data = [m4_WISM{1} m4_WISM{2}]; % mean PFC
% sub2exc = 1; 
% data.data(sub2exc, :) = []; 

% data.data = [m4_WISM{3} m4_WISM{4}]; % mean VVS
% sub2exc = [18 22]; 
% data.data(sub2exc, :) = []; 

% % % Direct contrast MEAN
diffPFC = [m4_WISM{1}-m4_WISM{2}];
diffVVS = [m4_WISM{3}-m4_WISM{4}];
diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
data.data = [diffPFC diffVVS]; 
data.data = abs(data.data)



%data.data = [m5_WISV{1} m5_WISV{2}]; % VAR PFC
%sub2exc = 1; 
%data.data(sub2exc, :) = []; 

% data.data = [m5_WISV{3} m5_WISV{4}]; % VAR VVS
% sub2exc = [18 22]; 
% data.data(sub2exc, :) = []; 


% data.data = [m7_WISMBC{1} m7_WISMBC{2}]; % MEAN between cat PFC
% sub2exc = 1; 
% data.data(sub2exc, :) = []; 


 % data.data = [m7_WISMBC{3} m7_WISMBC{4}]; % MEAN between cat VVS
 % sub2exc = 1; 
 % data.data(sub2exc, :) = []; 

% % Direct contrast MEAN BETWEEN CATEGORY
% diffPFC = [m7_WISMBC{1}-m7_WISMBC{2}];
% diffVVS = [m7_WISMBC{3}-m7_WISMBC{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 





% % % Direct contrast MEAN WITHIN CATEGORY
% diffPFC = [m6_WISMWC{1}-m6_WISMWC{2}];
% diffVVS = [m6_WISMWC{3}-m6_WISMWC{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 



% % % Direct contrast VARIANCE 
% diffPFC = [m5_WISV{1}-m5_WISV{2}];
% diffVVS = [m5_WISV{3}-m5_WISV{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 
% data.data = abs(data.data)

% diffPFC = [m8_WISVWC{1}-m8_WISVWC{2}];
% diffVVS = [m8_WISVWC{3}-m8_WISVWC{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 
% 
% diffPFC = [m9_WISVBC{1}-m9_WISVBC{2}];
% diffVVS = [m9_WISVBC{3}-m9_WISVBC{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 






figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2], data.data, 50, 'k'); hold on;
hb = plot(data.data'); hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);box on; 
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2); set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-.025 .055]); 

[h p ci ts] = ttest (data.data(:,1), data.data(:,2));
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)



exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% LOAD MEAN RDMS IN THE TWO REGIONS DURING MAINTENANCE
clear

f2load = 'vvs_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'all_clustinfo_VVS']);
idsClust = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);


for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    idsF1 = oneListIds(:, 3);
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (idsF1, idsF2);
    
    %rdmS         = squeeze(mean(mean(neuralRDMs(:, :, 1:6,12:24), 3),4)); %Theta cluster time period (checked in the out_real sigMH_real : Exact time period)
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    
    rdm = nan(60); 
    rdm(x3,x3) = rdmS; 
    meanRDM_VVS(subji,:,:) = rdm; 

    
end


clearvars -except meanRDM_VVS
clc
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
        
        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_PFC(subji,:,:) = rdm; 

        
    end
end



%% select only subjects with electrodes in both regions

pfc_rdms = meanRDM_PFC([2 3  5  9 10 11 12 14 15 16],:,:);
vvs_rdms = meanRDM_VVS([7 9 13 18 19 20 21 23 27 28],:,:);


for subji = 1:10
    
    pfc_rdmH = pfc_rdms(subji, : ,: ); 
    vvs_rdmH = vvs_rdms(subji, : ,: ); 

    pfc1 = vectorizeRDM_WM(pfc_rdmH); 
    vvs1 = vectorizeRDM_WM(vvs_rdmH); 

    isnPFC = isnan(pfc1); 
    isnVVS = isnan(vvs1); 

    pfc1(isnPFC|isnVVS) = []; 
    vvs1(isnPFC|isnVVS) = []; 

    [rhoALL(subji,:)] = corr(pfc1', vvs1', 'type', 's'); 
    
end


[h p ci ts] = ttest(rhoALL);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% plot one bar
clear data
data.data = [rhoALL]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.075 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth',2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth',2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM


%% plot mean RDM in each region

sub2excVVS = [18 22]; 
sub2excPFC = [1]; 

meanRDM_VVS(sub2excVVS, :,:) = []; 
meanRDM_PFC(sub2excPFC, :,:) = []; 


d2p1 = squeeze(mean(meanRDM_VVS, 'omitnan')); 
d2p1(logical(eye(size(d2p1, 1)))) = nan; 

d2p2 = squeeze(mean(meanRDM_PFC, 'omitnan')); 
d2p2(logical(eye(size(d2p2, 1)))) = nan; 

tiledlayout(1, 2)
nexttile
imagesc(d2p1); axis square; colorbar
set(gca, 'xTick', [], 'yTick', [], 'Fontsize', 20)
%set (gca, 'clim', [.35 .65])
nexttile
imagesc(d2p2); axis square; colorbar
set(gca, 'xTick', [], 'yTick', [], 'Fontsize', 20)
exportgraphics(gcf, 'myP.png', 'Resolution', 300)





%% RDM schematics: within category correlations

CM = kron(eye(6), ones(10));
CM(CM==0) = 2000; 
CM = tril(CM, -1);
CM(CM==0) = nan; 
CM(CM==2000) = 0; 

imagesc(CM); axis square; 
set(gca, 'xTick', [], 'yTick', [])

exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM

%% RDM schematics: between category correlations

CM = kron(eye(6), ones(10));
CM(CM==0) = 2000; 
CM = tril(CM, -1);
CM(CM==0) = nan; 
CM(CM==1) = 0;
CM(CM==2000) = 1; 

imagesc(CM); axis square; 
set(gca, 'xTick', [], 'yTick', [])

exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM



%% BL-NET FIT TO WITHIN ITEM / or BETWEEN CATEGORY CORRELATIONS PFC

clear, clc

subj_ch_fr = 7;
triu2u = 0; %select Within or between category correlations

% % load PFC RDMs
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


% compute the fits
for subji = 1:16
    
    f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat';
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
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
    
    BLDist = ACTH(CM==triu2u);
    %BLDist = vectorizeRDM_WM(ACTH);
   
    ACT = squeeze(meanRDM{subji, 1}); 
    neuDist = ACT(CM==triu2u); 
    %neuDist = vectorizeRDM_WM(ACT); 

    rhoALLPFC(subji, :) = corr(neuDist, BLDist, 'type', 's'); 
    
end


sub2exc = 1; 
rhoALLPFC(sub2exc) = []; 
[h p ci ts] = ttest(rhoALLPFC);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);



%% plot one bar
clear data
data.data = [rhoALLPFC]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.15 .2])
plot(get(gca,'xlim'), [0 0],'k','lineWidth',2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth',2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM

%% BL-NET FIT TO WITHIN ITEM / or BETWEEN CATEGORY CORRELATIONS VVS

clearvars -except rhoALLPFC
clc
subj_ch_fr = 17;
triu2u = 1; %Within or between category correlations

f2load = 'vvs_M123_[]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'all_clustinfo_VVS']);
idsClust = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);


for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    idsF1 = oneListIds(:, 3);
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (idsF1, idsF2);
    
    %rdmS         = squeeze(mean(mean(neuralRDMs(:, :, 1:6,12:24), 3),4)); %Theta cluster time period (checked in the out_real sigMH_real : Exact time period)
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
    meanRDM{subji, 1} = rdmS; 
    meanRDM{subji, 2} = oneListIds(:, 3); 
    
end

allF2sav = {'BLNETi_pfc_E123_[32]_3-54_0_0_0_1_.1_5_1.mat'; 'BLNETi_pfc_E123_[40]_3-54_0_0_0_1_.1_5_1.mat'; 'BLNETi_pfc_E123_[48]_3-54_0_0_0_1_.1_5_1.mat'}; 
% compute the fits
for subji = 1:size(meanRDM, 1)

    for layi = 1:3
        
        f2sav = allF2sav{layi}; 
        
        cfg = getParams(f2sav);
        paths = load_paths_WM(cfg.brainROI, cfg.net2load);
        [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
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
        
        BLDist = ACTH(CM==triu2u);
        %BLDist = vectorizeRDM_WM(ACTH);
       
        ACT = squeeze(meanRDM{subji, 1}); 
        neuDist = ACT(CM==triu2u); 
        %neuDist = vectorizeRDM_WM(ACT); 
    
        rhoALL(subji, layi, :) = corr(neuDist, BLDist, 'type', 's'); 
    end
    
end

rhoALLVVS = squeeze(mean(rhoALL, 2)); 

sub2exc = [18 22]; 
rhoALLVVS(sub2exc) = []; 
[h p ci ts] = ttest(rhoALLVVS);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% plot one bar
clear data
data.data = [rhoALLVVS]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.15 .3])
plot(get(gca,'xlim'), [0 0],'k','lineWidth',2);
box on

[h p ci ts] = ttest (data.data(:,1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

set(gca, 'LineWidth',2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM


%% contrasts VVS and PFC for 10 subjects

rhoALLPFC = rhoALLPFC([2 3  5  9 10 11 12 14 15 16],:,:);
rhoALLVVS = rhoALLVVS([7 9 13 18 19 20 21 23 27 28],:,:);

[h p ci ts] = ttest (rhoALLPFC, rhoALLVVS);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 1 comment 2 - ITEM SPECIFIC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%

M = zeros(110); 
M(1:35, 1:35) = 1; M(36:65, 36:65) = 1; M(66:86, 66:86) = 1; M(87:end, 87:end) = 1; 
imagesc(M); axis square; 
set(gca, 'xtick', [], 'ytick', [],   'xticklabels', [],  'yticklabels', []); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%%

M = zeros(110); 
%M(1:35, 1:35) = rand(35); M(36:65, 36:65) = rand(30);  M(66:86, 66:86) = rand(21); M(87:end, 87:end) = rand(24); 
M(1:65, 1:65) = rand(65);M(66:110, 66:110) = rand(45); 
M(M==0) = nan;
imagesc(M); axis square; 
set(gca, 'xtick', [], 'ytick', [],   'xticklabels', [],  'yticklabels', []); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% Schematic (item model)

CM = kron(eye(6), ones(10));
CM(CM==0) = 2000; 
CM = tril(CM, -1);
CM(CM==0) = nan; 
CM(CM==2000) = 0; 

imagesc(CM); axis square; 
set(gca, 'xTick', [], 'yTick', [])

exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM




%% compute CCI in PFC cluster
clear, clc
paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
load ([paths.results.CCI 'pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1_CCI.mat']);
clustinfo = allClustInfo{7}; 
for subji = 1:size(CCI, 1)
    cciSubj = squeeze(CCI(subji, :, 1:40)); 
    cciClust(subji, :) = mean(cciSubj(clustinfo.PixelIdxList{7}), 'all');
end
sub2exc = [1]
cciClust(sub2exc) = []; 
[h p ci ts] = ttest (cciClust);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 2)]);


%% compute CCI in VVS cluster
clear, clc
paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'all_clustinfo_VVS']);
idsClust=intersect(intersect(allClustInfo{4}.PixelIdxList{14},allClustInfo{5}.PixelIdxList{25},'stable'), allClustInfo{6}.PixelIdxList{17},'stable');
load ([paths.results.CCI 'vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1_CCI.mat']);

for subji = 1:size(CCI, 1)
    cciSubj = squeeze(CCI(subji, :, 1:40)); 
    cciClust(subji, :) = mean(cciSubj(idsClust), 'all');
end
sub2exc = [18 22]
cciClust(sub2exc) = []; 
[h p ci ts] = ttest (cciClust);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


%% count how many items repetitions
clear , clc
f2sav = 'ITM_vvs_M123_[1]_3-54_0_0_1_0_.1_5_1'

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


for subji = 1:length(nnFit)

    ids = nnFit{subji, 2}; 
    oneListIds = cellfun(@(x) strsplit(x, ' '), ids, 'un', 0);
    oneListIds = double(string(cat(1, oneListIds{:})));
    idsF1 = oneListIds(:, 3);
    [v, w] = unique(idsF1, 'stable' );
    duplicate_indices = setdiff( 1:numel(idsF1), w ); 
    nRep(subji, :) = length(duplicate_indices); 

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 1 comment 3 - SLIDE 10 - Correct vs incorrect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  COMPARE CORRECT VS INCORRECT
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'Alex_vvs_M123_[1-8]_3-54_0_0_1_1_.1_5_1';
f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';

minSubCrit = 5; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,8, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1800 1300])
else
    set(gcf, 'Position', [100 100 670 1300])
end

for layi = 1:size(nnFit{2}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:40));
             ids = cellfun(@(x) strsplit(x, ' '), nnFit{subji, 2}, 'un', 0); %str2num does not work
             ids = double(string(cat(1, ids{:})));%str2num does not work
           elseif strcmp(cfg.period(1), 'E') % for clarity
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:15));
             ids= str2num(cell2mat(nnFit{subji, 2}));
           end
            idsCC = ids(:,20)==1; 
            idsIC = ids(:,20)==0; 
            idsCI = ids(:,19)==1; 
            idsII = ids(:,19)==0; 

            nTR(subji,:) = sum(idsIC==1); 
            nnH1(subji, : ,:) = squeeze(mean(avTR(:, idsCC,:,:), 2)); 
            nnH2(subji, : ,:) = squeeze(mean(avTR(:, idsIC,:,:), 2)); 

       end
    end
    
    sub2exc2 = find(nTR<minSubCrit); 
    sub2exc = union(sub2exc, sub2exc2);
    nnH1(sub2exc, :, :) = []; 
    nnH1 = squeeze(nnH1);

    nnH2(nnH2==inf) = nan; 
    nnH2(sub2exc, :, :) = []; 
    nnH2 = squeeze(nnH2);

    nnH1L(layi, :, :,:) = nnH1; 
    nnH2L(layi, :, :,:) = nnH2; 

    [h p ci ts] = ttest(nnH1, nnH2);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH1 - nnH2, 'omitnan'));
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    clear allTObs
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs = 0;
    end

    [max2u id] = max(allTObs);
    tObs(layi) = allTObs(id); 
    
    freqs = 1:520; 
    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2)*10; 
        h = zeros(52, 40); 
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 1
            h(clustinfo.PixelIdxList{11}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 3
            h(clustinfo.PixelIdxList{13}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 4
            h(clustinfo.PixelIdxList{13}) = 1;
        end

    else
        times = 1:150; 
    end
    myCmap = colormap(brewermap([],'*spectral'));
    colormap(myCmap)
    contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 400], 'clim', [-5 5], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
    end
end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 
close all 


%% check category and item model in these clusters
clear, clc
% % first load neural RDMs in VVS
f2load = 'vvs_M123_[]_3-54_0_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'allClustInfo_CORR-INC_VVS.mat']);
ids7clust = [11 19 13 13 13 15 22]; 


%%
clear areThereExemp
for subji = 1:nSubjs

    if ~isempty(allNeuralRDMS{subji, 1})
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        
        CM = squeeze(load_CATMODEL_activ(ids)); 
        %CM = squeeze(load_ITMODEL_activ(ids)); 
        CM = vectorizeRDM_WM(CM); 
        nanIds = isnan(CM); 
        CM(nanIds) = []; 
        areThereExemp(subji,:) = ~isempty(find(CM));
        
        for layi = 1:7
            idsClust = ids7clust(layi); 
            rdm         = mean(neuralRDM1(:, :,idsClust ), 3, 'omitnan'); 
            rdm = vectorizeRDM_WM(rdm); 
            rdm(nanIds) = []; 
            allR(subji, layi) = corr(CM, rdm, 'type', 's');
        end
    end
end

%%
sub2exc = [18 22]; 
allR(sub2exc, :) = []; 


[h p ci t] = ttest (allR);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);






%%