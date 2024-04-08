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



%data.data = [m4_WISM{1} m4_WISM{2}]; % mean PFC
%sub2exc = 1; 
%data.data(sub2exc, :) = []; 

data.data = [m4_WISM{3} m4_WISM{4}]; % mean VVS
sub2exc = [18 22]; 
data.data(sub2exc, :) = []; 

%data.data = [m5_WISV{1} m5_WISV{2}]; % VAR PFC
%sub2exc = 1; 
%data.data(sub2exc, :) = []; 

% data.data = [m5_WISV{3} m5_WISV{4}]; % VAR VVS
% sub2exc = [18 22]; 
% data.data(sub2exc, :) = []; 


% % % Direct contrast MEAN
% diffPFC = [m4_WISM{1}-m4_WISM{2}];
% diffVVS = [m4_WISM{3}-m4_WISM{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 
% data.data = abs(data.data)

% % % Direct contrast MEAN WITHIN CATEGORY
% diffPFC = [m6_WISMWC{1}-m6_WISMWC{2}];
% diffVVS = [m6_WISMWC{3}-m6_WISMWC{4}];
% diffPFC = diffPFC([2 3  5  9 10 11 12 14 15 16]);
% diffVVS = diffVVS([7 9 13 18 19 20 21 23 27 28]);
% data.data = [diffPFC diffVVS]; 

% % % Direct contrast MEAN BETWEEN CATEGORY
 % diffPFC = [m7_WISMBC{1}-m7_WISMBC{2}];
 % diffVVS = [m7_WISMBC{3}-m7_WISMBC{4}];
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




%data.data = [m5_WISV_E m5_WISV_M];
%data.data = [m6_WISMWC_E, m6_WISMWC_M]; 
%data.data = [m7_WISMBC_E, m7_WISMBC_M]; 
%data.data = [m8_WISVWC_E, m8_WISVWC_M]; 
%data.data = [m9_WISVBC_E, m9_WISVBC_M]; 



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

    rhoALL(subji, :) = corr(neuDist, BLDist, 'type', 's'); 
    
end

sub2exc = 1; 
rhoALL(sub2exc) = []; 
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
set(gca, 'ylim', [-.15 .2])
plot(get(gca,'xlim'), [0 0],'k','lineWidth',2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth',2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM

%% BL-NET FIT TO WITHIN ITEM / or BETWEEN CATEGORY CORRELATIONS VVS

clear, clc
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

rhoALLM = squeeze(mean(rhoALL, 2)); 

sub2exc = [2 18]; 
rhoALLM(sub2exc) = []; 
[h p ci ts] = ttest(rhoALLM);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% plot one bar
clear data
data.data = [rhoALLM]; 

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

[h p ci t] = ttest (data.data(:,1));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

set(gca, 'LineWidth',2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);%% plot average RDM



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 1 comment 2 - SLIDE 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 2 comment 1 - SLIDE 14 - Correct vs incorrect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot all layers BANDS > recover method for plotting the category model in individual bands

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc


f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-8_1_0_0_0_.1_5_1'


cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
    %set(gcf, 'Position', [100 100 1800 1000])
end

for layi = 1:size(nnFit{2}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:15));
             %nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:134));
           end
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
        %times = 1:134;
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    hb = h; hb(h==1) = 0; hb(h==0) = nan; 
    
    plot(times, t); hold on; %colorbar
    plot(times, hb, linewidth=2)
    
    
    if strcmp(cfg.period(1), 'M')
        
    else
        
    end
    
end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 




%% compute clusters in each permutation BANDS for all layers
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%f2sav =  'RNN_pfc_M123_[8-8-56]_13-29_1_1_0_0_.1_5_1_1000p.mat'; 
%f2sav = 'CORrt_vvs_E123_[1-8]_39-54_1_0_0_0_.1_5_1_1000p.mat';
%f2sav = 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1_1000p.mat';


f2sav = [f2sav '_1000p.mat'];

cd ([paths.results.DNNs])
load(f2sav);


if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end


f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{13}(1:end-5))));


for permi = 1:nPerm
    
 

    for layi = 1:7
        dataP = squeeze(nnFitPerm(permi, :, layi,:));
        dataP(sub2exc, :,:) = []; 
        [h p ci ts] = ttest(dataP);
        h = squeeze(h); t = squeeze(ts.tstat);
    
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        %sort all clusters 
        if exist('allSTs')
            [max2u id] = max(allSTs);
            max_clust_sum_perm(permi, layi, :) = allSTs(id); 
        else
            max_clust_sum_perm(permi, layi, :) = 0; 
        end
    end

end

%% compute p value bands for all layers FREQ RES
clc
clear p

for layi = 1:size(allTObs, 1)
    clear mcsR mcsP
    mcsR = allTObs(layi); 
    mcsP = squeeze(max_clust_sum_perm(:,layi));
    
    %allAb = mcsP(abs(mcsP) > abs(mcsR));
    allAb = mcsP(mcsP > mcsR);
    p(layi, :) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
    
end

p_ranked = p; p_ranked(isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));
p

%% p last layer only
p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))


cd (paths.github)




%% load files to plot all CATEGORY RESULTS  -- IN one PLOT ONLY 

clear, clc 

listF2sav = {       
    
                    % 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';


                    'CAT_pfc_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';


%                     'CAT_vvs_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';

                    % 'CAT_vvs_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';

                    % 'BLNETi_pfc_M123_[8-8-56]_3-8_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_9-12_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_13-29_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_30-38_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_39-54_1_0_0_0_.1_5_1'


             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi allFnnH
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI,  cfg.net2load);
    filelistSess = getFiles([paths.results.DNNs]);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.results.DNNs listF2sav{listi}]);   
        
        nnH = cat(1, nnFit{:,1});
        
    end

    allFnnH(listi, :, :) = nnH; 


end
cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    %sub2exc = [18 22];
    sub2exc = [18];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [];
end
%allFnnH(:, sub2exc, :) = []; 

%% plot


for freqi = 1:size(allFnnH, 1)
    if strcmp(cfg.period(1), 'M')
        nnH = allFnnH(freqi, :, 1:40);
    else
        nnH = allFnnH(freqi, :, 1:15);
    end
    nnH = squeeze(nnH);
    nnH(sub2exc, :, :) = []; 

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(freqi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(freqi, :, :) = 0;
    end
    
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    clear hbl
    hbL(freqi, :) = hb; 
    
    mART(freqi, :) = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    

end


if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 500 350])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    times = 1:40;
    x = (-.015:-.003:-.028)';
    %hbL([1:5], :) = nan;     
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    %colormap(jet(5));
    times = 1:15;
    plot(times, mART, 'Linewidth', 5); hold on; 
    x = (-.015:-.0053:-.039)';
    %hbL(4, 1:3) = nan; 

    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 5); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 12, 'ylim', [-.04 .15]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);
    %legend

end


set(gca, 'ColorOrder', myCmap)


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 2 comment 3 - SLIDE 14 - Correct vs incorrect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANOVA comparing fits in ALEXNET BL-NET (LAYER 7 vs 8) AND CATEGORY MODEL > IN SELECTED ROI
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

timeBins =6:12; 
freqBins = 21:29; 

f2sav1 = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitALEX = nnFit; 

f2sav3 = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav3);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav3 '.mat']);
nnFitCAT = nnFit; 

load ([paths.results.clusters 'clustinfo_PFC_px2']);
idsClust = clustinfo.PixelIdxList{2}; 


clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
    if ~isempty(nnFitBLNET{subji, 1})
        nnHBLNET(subji, : ,:) = mean(mean(atanh(nnFitBLNET{subji, 1}(7,freqBins,timeBins)), 2), 3);
        nnHALEX(subji, : ,:) = mean(mean(atanh(nnFitALEX{subji, 1}(8,freqBins,timeBins)), 2), 3);
        nnHCAT(subji, : ,:) = mean(mean(atanh(nnFitCAT{subji, 1}(1,freqBins,timeBins)), 2), 3);
    end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHALEX([18 22], :, :) = []; 
    nnHCAT([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHALEX([1], :, :) = []; 
    nnHCAT([1], :, :) = []; 
end


d4ANOVA = [nnHBLNET; nnHALEX ;nnHCAT]; 
d4ANOVA(:,2) = [ones(1,15) ones(1,15)*2 ones(1,15)*3];
d4ANOVA(:,3) = [1:15 1:15 1:15];
[p F] = RMAOV1(d4ANOVA);


%% plot three bar
clear data

data.data = [nnHBLNET nnHALEX nnHCAT]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter(1:3, data.data, 50, 'k'); hold on;
%hb = plot(data.data'); hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);box on; 
set(gca,'XTick',1:3,'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 4] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2); set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-.025 .055]); 



exportgraphics(gcf, 'allM.png', 'Resolution', 300);





%% 
[h p ci ts] = ttest(nnHBLNET, nnHALEX);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])


[h p ci ts] = ttest(nnHBLNET, nnHCAT);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])

[h p ci ts] = ttest(nnHALEX, nnHCAT);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])




%% WITHIN SUBJECTS DESING IN MATLAB

d4ANOVA = [[1:15]' nnHBLNET nnHALEX nnHCAT]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'nnHBLNET' 'nnHALEX' 'nnHCAT'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Model'});
withinDesign.Model = categorical(withinDesign.Model);
% create the repeated measures model and do the anova
rm = fitrm(T,'nnHBLNET-nnHCAT ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));




%% Ttest at every time-frequency point comparing fits in ALEXNET VS BL-NET (LAYER 7 vs 8)

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitALEX = nnFit; 

tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end


ax1 = nexttile;
clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
        if ~isempty(nnFitBLNET{subji, 1})
           if strcmp(cfg.period(1), 'M')
                nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(7,:,1:40));
           else
             nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(7,:,1:15));
           end
        end 
end
for subji = 1:length(nnFitALEX)
   if strcmp(cfg.period(1), 'M')
     nnHALEX(subji, : ,:) = atanh(nnFitALEX{subji, 1}(7,:,1:40));
   else
     nnHALEX(subji, : ,:) = atanh(nnFitALEX{subji, 1}(7,:,1:15));
   end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHALEX([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHALEX([1], :, :) = []; 
end


nnHBLNET = squeeze(nnHBLNET);
nnHALEX = squeeze(nnHALEX);

[h p ci ts] = ttest(nnHBLNET, nnHALEX);
h = squeeze(h); t = squeeze(ts.tstat); 
h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

freqs = 1:52; 
clustinfo = bwconncomp(h);
allClustInfo{1} = clustinfo; 

% store allTOBS
if ~isempty(clustinfo.PixelIdxList)
    for pixi = 1:length(clustinfo.PixelIdxList)
         allTObs(1, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
         allTObs1(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
    end
else
    allTObs(1, :, :) = 0;
end

if exist('allTObs1')
    [max2u id] = max(abs(allTObs1));
    tObs = allTObs1(id); 
end


if strcmp(cfg.period(1), 'M')
    times = 1:size(t, 2); 
else
    times = 1:15; 
end
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);

if strcmp(cfg.period(1), 'M')
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'xlim', [1 40], 'clim', [-3 3], 'FontSize', 10);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 8, 'clim', [-5 5]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    
end


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% Permutations 
nPerm = 1000; 
timePeriod = [6:15]; 

clear max_clust_sum_perm allTObs
for permi = 1:nPerm 

    junts = cat(1, nnHBLNET, nnHALEX);
    junts = junts (randperm(size(junts, 1)), :, timePeriod);
    nnHBLNETPerm = junts(1:floor(size(junts, 1)/2), :, :); 
    nnHALEXPerm = junts(floor(size(junts, 1)/2)+1:end, :, :); 
    
    [h p ci ts] = ttest(nnHBLNETPerm, nnHALEXPerm);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

    clustinfo = bwconncomp(h);
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(permi,:) = 0;
    end
    
    if exist('allTObs')
        [max2u id] = max(abs(allTObs));
        max_clust_sum_perm(permi,:) = allTObs(id); 
    end


end

allAB = max_clust_sum_perm(max_clust_sum_perm> tObs);
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

disp(['p = ' num2str(p)])


%% plot nicely ALEXNET layer 7


freqs = 1:520; 
times = 1:size(t, 2)*10; 
h = zeros(52, 40); 
h(clustinfo.PixelIdxList{5}) = 1; 

figure(); set(gcf, 'Position', [100 100 300 500])
myCmap = colormap(brewermap([],'*spectral'));
colormap(myCmap)
contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);


set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 400], 'clim', [-5 5], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
set(gca, 'xlim', [1 150])


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% plot nicely ALEXNET layer 8

freqs = 1:520; 
times = 1:size(t, 2)*10; 
h = zeros(52, 40); 
h(clustinfo.PixelIdxList{6}) = 1; 

figure(); set(gcf, 'Position', [100 100 300 500])
myCmap = colormap(brewermap([],'*spectral'));
colormap(myCmap)
contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);


set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 400], 'clim', [-5 5], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
set(gca, 'xlim', [1 150])


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% ANOVA comparing fits in ALEXNET BL-NET AND CATEGORY MODEL > IN SELECTED ROI VVS LAYER 5
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

timeBins = 26:36; 
freqBins = 7:11; % alpha range

f2sav1 = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'Alex_vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitALEX = nnFit; 

f2sav3 = 'CAT_vvs_M123_[1]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav3);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav3 '.mat']);
nnFitCAT = nnFit; 


clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
    if ~isempty(nnFitBLNET{subji, 1})
        nnHBLNET(subji, : ,:) = mean(mean(atanh(nnFitBLNET{subji, 1}(5,freqBins,timeBins)), 2), 3);
        nnHALEX(subji, : ,:) = mean(mean(atanh(nnFitALEX{subji, 1}(5,freqBins,timeBins)), 2), 3);
        nnHCAT(subji, : ,:) = mean(mean(atanh(nnFitCAT{subji, 1}(1,freqBins,timeBins)), 2), 3);
    end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHALEX([18 22], :, :) = []; 
    nnHCAT([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHALEX([1], :, :) = []; 
    nnHCAT([1], :, :) = []; 
end


 
d4ANOVA = [nnHBLNET; nnHALEX ;nnHCAT]; 
d4ANOVA(:,2) = [ones(1,26) ones(1,26)*2 ones(1,26)*3];
d4ANOVA(:,3) = [1:26 1:26 1:26];

[p F] = RMAOV1(d4ANOVA);


%% plot three bar
clear data

data.data = [nnHBLNET nnHALEX nnHCAT]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter(1:3, data.data, 50, 'k'); hold on;
%hb = plot(data.data'); hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);box on; 
set(gca,'XTick',1:3,'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 4] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2); set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-.025 .055]); 



exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% 
[h p ci ts] = ttest(nnHBLNET, nnHALEX);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])

[h p ci ts] = ttest(nnHBLNET, nnHCAT);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])

[h p ci ts] = ttest(nnHALEX, nnHCAT);
disp(['t= ' num2str(ts.tstat) ', p = ' num2str(p)])




%% WITHIN SUBJECTS DESING IN MATLAB

d4ANOVA = [[1:26]' nnHBLNET nnHALEX nnHCAT]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'nnHBLNET' 'nnHALEX' 'nnHCAT'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Model'});
withinDesign.Model = categorical(withinDesign.Model);
% create the repeated measures model and do the anova
rm = fitrm(T,'nnHBLNET-nnHCAT ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 2 comment 4 - SLIDE 18 - RANDOM SUBSAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD all conditions
clear
currFolder = pwd; 
cd D:\_WM\analysis\pattern_similarity\vvs\100ms\EM2\3-8Hz
paths = load_paths_WM('vvs', 'none');

contrasts = {
              
              'DISC_EM2' 'DIDC_EM2';
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

cd(pwd)

%% COMPUTE T-OBS FOR ALL SUBJECTS
clc
clearvars -except DISC_EM2 DIDC_EM2

%define conditions 
cond1 = 'DISC_EM2';
cond2 = 'DIDC_EM2';


sub2exc2 = [18 22]; 

cond1B = eval(cond1);
cond2B = eval(cond2);
 
cfg             =       [];
cfg.subj2exc    =       sub2exc2;% vvs;
cfg.clim        =       [-.01 .01];
cfg.climT       =       [-7 7]; %color scale for t-map
cfg.saveimg     =       1;
cfg.res         =       '100_perm'; %'100_perm'; '100_norm'
cfg.cut2        =       '1-4'; %1-1 1-4 4-4
cfg.cond1       =       cond1;
cfg.cond2       =       cond2;
cfg.lwd1        =       2; %baseline 
cfg.lwd2        =       2; %significant outline
cfg.remClust    =       0; 
cfg.plot1clust  =       0;  
cfg.clust2plot  =       [1];  %VVS > 3-8: 4; 9-12: 6; 13-29: 8; 30-75: 7; 75-150: [3 17]; vector of pixels to print
cfg.all_cond1   =       cond1B; 
cfg.all_cond2   =       cond2B; 
cfg.alpha       =       0.05; 


[out_real]   = plot_reinst_map_wm_no_plot(cfg);
 
mDiff = squeeze(mean(mean(out_real.meanReal_cond1(:, 1:8, 1:15) - out_real.meanReal_cond2(:, 1:8, 1:15), 2), 3));
[h p ci ts] = ttest(mDiff); 
tObs = ts.tstat;



%save('1000p16S.mat','alltObs', 'allPs', 'allTs')

%% RUN this same analysis by randomly excluding subjects in VVS
clc
clearvars -except DISC_EM2 DIDC_EM2

nPerm = 1000; 

%define conditions 
cond1 = 'DISC_EM2';
cond2 = 'DIDC_EM2';


for permi = 1:nPerm
    permi 

    ids = [1:17 19:21 23:28]; 
    ids2 = randperm(26, 11); 
    ids3 = ids(ids2); 
    sub2exc2 = [18 22 ids3]; 
    
    cond1B = eval(cond1);
    cond2B = eval(cond2);
     
    cfg             =       [];
    cfg.subj2exc    =       sub2exc2;% vvs;
    cfg.clim        =       [-.01 .01];
    cfg.climT       =       [-7 7]; %color scale for t-map
    cfg.saveimg     =       1;
    cfg.res         =       '100_perm'; %'100_perm'; '100_norm'
    cfg.cut2        =       '1-4'; %1-1 1-4 4-4
    cfg.cond1       =       cond1;
    cfg.cond2       =       cond2;
    cfg.lwd1        =       2; %baseline 
    cfg.lwd2        =       2; %significant outline
    cfg.remClust    =       0; 
    cfg.plot1clust  =       0;  
    cfg.clust2plot  =       [1];  %VVS > 3-8: 4; 9-12: 6; 13-29: 8; 30-75: 7; 75-150: [3 17]; vector of pixels to print
    cfg.all_cond1   =       cond1B; 
    cfg.all_cond2   =       cond2B; 
    cfg.alpha       =       0.05; 
    
    
    [out_real]   = plot_reinst_map_wm_no_plot(cfg);
     
    
    % %%perm
    % cfg_perm                    =       [];
    % cfg.runperm                 =       1;  
    % cfg_perm.n_perm             =       1000; 
    % cfg_perm.savePerm           =       1;
    % cfg_perm.out_real           =       out_real;
    % cfg_perm.pval               =       0.05;
    % cfg_perm.cond1              =       cond1;
    % cfg_perm.cond2              =       cond2;
    % cfg_perm.ids_all_cond       =       [];
     
    %allTs(permi, :) = mean(out_real.sigMT_real(1:8, 6:15), 'all');
    %tObs = squeeze(mean(mean(out_real.meanReal_cond1(:, 6:13, 6:15) - out_real.meanReal_cond2(:, 6:13, 6:15), 2), 3));
    mDiff = squeeze(mean(mean(out_real.meanReal_cond1(:, 1:8, 1:15) - out_real.meanReal_cond2(:, 1:8, 1:15), 2), 3));
    [h p ci ts] = ttest(mDiff); 
    allPs(permi, :) = p;
    allTs(permi, :) = ts.tstat;

    
    % [out_perm] = myPerm(cfg_perm);
    % if ~isempty (out_perm.max_clust_sum_real)
    %     max_clust_sum = out_perm.max_clust_sum;
    %     obs = max(out_real.all_clust_tsum_real(:,1));
    %     %allAb = max_clust_sum(abs(max_clust_sum) > obs);
    %     allAb = max_clust_sum((max_clust_sum) > obs);
    %     p(permi, :) = (1 - (cfg_perm.n_perm - (length (allAb) ) )  /cfg_perm.n_perm) + (1/cfg_perm.n_perm);
    %     %disp (['p = ' num2str(p)]);
    % end

end

%save('1000p16S.mat','alltObs', 'allPs', 'allTs')


%% count how many times there is a significant effect when VVS N = 15
numberOfTimes = sum(allPs < 0.05);
perc2report = (numberOfTimes/nPerm)*100

%% 
figure()
histogram(allTs); hold on
set(gca, fontsize=16)
plot ([1.96, 1.96], ylim, lineWidth=2)
exportgraphics(gcf, 'myP.png', 'Resolution', 300);


%% RANDOM SUBSAMPLE SUBJECTS DNN FITS
clear, clc

nPerm = 1000; 
lay2u = 6; 
f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';


cfg = getParams(f2sav);

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);

%load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
load ([paths.results.clusters 'all_clustinfo_VVS']);

% compute overlapping cluster in layers 4 5 and 6
%c2u = unique([allClustInfo{4}.PixelIdxList{14} ; allClustInfo{5}.PixelIdxList{25} ; allClustInfo{6}.PixelIdxList{17}]);
c2u = allClustInfo{6}.PixelIdxList{17};

for permi = 1:nPerm

    
    clear nnHClust_pfc7 nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6
    
    ids = [1:17 19:21 23:28]; 
    ids2 = randperm(26, 11); 
    ids3 = ids(ids2); 
    sub2exc2 = [18 22 ids3]; 


    clear nnH
    for subji = 1:length(nnFit)
        if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(lay2u,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(lay2u,:,1:15));
           end
        end
    end
    
    nnH(sub2exc2, :, :) = []; 


    for subji = 1:size(nnH, 1)
        nnHSubj = squeeze(nnH(subji, :, :)); 
        if strcmp(cfg.period(1), 'M')
            nnHClust(subji, :) = mean(nnHSubj(c2u)); 
        elseif strcmp(cfg.period(1), 'E')
            nnHClust_vvsE1(subji, :) = mean(nnHSubj(allClustInfo{1}.PixelIdxList{1})); 
            nnHClust_vvsE2(subji, :) = mean(nnHSubj(allClustInfo{2}.PixelIdxList{1})); 
            nnHClust_vvsE3(subji, :) = mean(nnHSubj(allClustInfo{3}.PixelIdxList{1})); 
            nnHClust_vvsE4(subji, :) = mean(nnHSubj(allClustInfo{4}.PixelIdxList{1})); 
            nnHClust_vvsE5(subji, :) = mean(nnHSubj(allClustInfo{5}.PixelIdxList{1})); 
            nnHClust_vvsE6(subji, :) = mean(nnHSubj(allClustInfo{6}.PixelIdxList{1})); 
            nnHClust_vvsE7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{1})); 
        end
    end
   


    
    [h p ci ts] = ttest(nnHClust);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %disp (['t = ' num2str(ts.tstat) ',' ' p = ' num2str(p)]);

    allPs(permi, :) = p; 
    allTs(permi, :) = ts.tstat; 

end

%% count how many times there is a significant effect when VVS N = 15
numberOfTimes = sum(allPs < 0.05);
perc2report = (numberOfTimes/nPerm)*100

%% 
figure()
histogram(allTs); hold on
set(gca, fontsize=16)
plot ([1.96, 1.96], ylim, lineWidth=2)
exportgraphics(gcf, 'myP.png', 'Resolution', 300);

%% RANDOM SUBSAMPLE SUBJECTS ANALYSIS DURING WHOLE ENCODING TIME PERIOD
clear, clc

nPerm = 1000; 
lay2u = 6; 
f2sav = 'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';


cfg = getParams(f2sav);

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


for permi = 1:nPerm

    
    clear nnHClust_pfc7 nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6
    
    ids = [1:17 19:21 23:28]; 
    ids2 = randperm(26, 11); 
    ids3 = ids(ids2); 
    sub2exc2 = [18 22 ids3]; 


    clear nnH
    for subji = 1:length(nnFit)
        if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(lay2u,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(lay2u,:,1:15));
           end
        end
    end
    
    nnH(sub2exc2, :, :) = []; 


    for subji = 1:size(nnH, 1)
        nnHSubj = squeeze(nnH(subji, :, :)); 
        if strcmp(cfg.period(1), 'M')
            nnHClust(subji, :) = mean(nnHSubj(c2u)); 
        elseif strcmp(cfg.period(1), 'E')
            nnHClust(subji, :) = mean(nnHSubj(:, 6:15), 'all'); 
        end
    end
   


    
    [h p ci ts] = ttest(nnHClust);
    h = squeeze(h); t = squeeze(ts.tstat); 

    allPs(permi, :) = p; 
    allTs(permi, :) = ts.tstat; 

end

%% count how many times there is a significant effect when VVS N = 15
numberOfTimes = sum(allPs < 0.05);
perc2report = (numberOfTimes/nPerm)*100

%% 
figure()
histogram(allTs); hold on
set(gca, fontsize=16)
plot ([1.96, 1.96], ylim, lineWidth=2)
exportgraphics(gcf, 'myP.png', 'Resolution', 300);






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 1 - SLIDE 19 - WITHIN AND BETWEEN CATEGORY CORRELATIONS
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

%%
clear rhoALL rhoNoiALL

nPerm = 1000;
witBet = 1; 
noiseLevels = [0 .1 .5 1 3 5 10];

subj_ch_fr = 7;
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
    
    BLDist = ACTH(CM==witBet);
    %BLDist = vectorizeRDM(ACTH);
   
    ACT = squeeze(meanRDM{subji, 1}); 
    neuDist = ACT(CM==witBet); 
    %neuDist = vectorizeRDM(ACT); 

    rhoALL(subji, :) = corr(neuDist, BLDist, 'type', 's'); 
    
    
    for noisei = 1:length(noiseLevels)
    
        mu = noiseLevels(noisei); 
    
        for permi = 1:nPerm
            CMn = kron(eye(6), ones(10));
            CMn = CMn + normrnd(0, mu,[60,60]);
            %CMn = normalize(CMn, 'range');
            %CMn(logical(eye(size(CMn)))) = 1; 
            %CMn = tril(CMn,-1)+tril(CMn)'; 
            %CCI = compute_CCI_WM(CMn);
            CMn = CMn(ids4, ids4); 
            CMn = CMn(CM==witBet);
            rhoNoiALL(subji, permi, noisei, :) = corr(neuDist, CMn, 'type', 's'); 
        end
    end
end

sub2exc = 1; 
rhoALL(sub2exc, :) = []; 
rhoNoiALL(sub2exc, :,:) = []; 

%% check number of significant correaltions 

[h p ci ts] = ttest(rhoALL)


[h1 p1 ci ts] = ttest(rhoNoiALL)
h1 = squeeze(h1); p1 = squeeze(p1); 
t1 = squeeze(ts.tstat); 
percHs = sum (h1)

%% 
figure(); set(gcf, 'Position', [100 100 1200 300])
tiledlayout (1, 6)
for noisei = 2:7
    nexttile
    histogram(t1(:, noisei)); hold on
    plot ([1.96, 1.96], ylim, lineWidth=2)

end

exportgraphics(gcf, 'myP.png', 'Resolution', 300);



%% LOAD MAINTENANCE DATA (PFC CLUSTER)
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



%%

clear rhoALL 
witBet = 1; 
subj_ch_fr = 7;

for subji = 1:16
        
    
    f2sav = 'Alex_pfc_E123_[7]_3-54_0_0_0_1_.1_5_1.mat'; 
    %f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
    
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    %[ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    [ACT] = load_alex_activ(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet


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






























%%