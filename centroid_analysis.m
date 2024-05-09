%% Correlate with BL-NET after excluding within and between category correlations
% Process and plot mean RDM in the PFC cluster (MAINTENANCE) 
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
subj_ch_fr = 7;

for subji = 1:16
        
    
    f2sav = 'Alex_pfc_E123_[7]_3-54_0_0_0_1_.1_5_1.mat'; 
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
    
    %BLDist = ACTH(CM==0);
    BLDist = vectorizeRDM_WM(ACTH);
   
    ACT = squeeze(meanRDM{subji, 1}); 
    %neuDist = ACT(CM==0); 
    neuDist = vectorizeRDM_WM(ACT); 

    rhoALL(subji, :) = corr(neuDist, BLDist, 'type', 's'); 
    
end

sub2exc = 1; 
rhoALL(sub2exc) = []; 
[h p ci ts] = ttest(rhoALL);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

%% correlate BL-NET dist with category model 
CM2T = vectorizeRDM_WM(CM); 
[rho1 p1] = corr(BLDist, CM2T, 'type', 's')


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
set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);







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



%% plot one bar
clear data
data.data = [rhoALLCAT]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 24, 'linew',1, 'xlim', [0 2] );
set(gca, 'ylim', [-.6 .85])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% 

[h p ci ts] = ttest(rhoALLCAT, rhoALLCATCM)




%% COMPUTE VARIANCE FOR PERMUTATIONS

nPerm = 1000; 
tril2u = 0; 

clear neuralDistances neuralDistancesPerm neuralCatDist
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

    
    neuralCatDistB = averageCAT_WM(ACT, ind); 
    neuralCatDist(subji,:) = neuralCatDistB(tril(true(length(neuralCatDistB)), 0));

    
    for permi = 1:nPerm
        
        ACTH = ACT; 
        idsperm = randperm(size(ACTH, 1)); 
        ACTH = ACTH(idsperm, idsperm); 
        neuralDistancesPerm(permi, subji,:) = vectorizeRDM(ACTH); 
        %neuralDistancesPerm(permi, subji,:) = ACTH(CM==tril2u); 
        
        neuralCatDistBP = averageCAT_WM(ACTH, ind); 
        neuralCatDistPerm(permi, subji,:) = neuralCatDistBP(tril(true(length(neuralCatDistBP)), 0));
        
    end

end

sub2exc = [1]; 
neuralDistances(sub2exc, :) = []; 
neuralCatDist(sub2exc, :) = []; 
neuralDistancesPerm(:, sub2exc, :) = []; 
neuralCatDistPerm(:, sub2exc, :) = []; 
%% plot one subject one row

%d2p1 = realV; 
%d2p2 = permV; 

d2p1 = std(neuralDistances,[],  1, 'omitnan'); 
d2p2 = squeeze(std(neuralDistancesPerm, [], 2, 'omitnan'));
%d2p3 = squeeze(mean(d2p2, 1, 'omitnan'));

cols = repelem([0.5, 0.5, 0.5], nPerm, 1); 
figure(); set(gcf, 'Position', [100 100 2000 200])
nIts = 20; 
scatter(1:nIts, d2p2(:,1:nIts), 50, 'r.'); hold on
scatter(1:nIts, d2p1(:,1:nIts), 200, 'k.')

%scatter(1:nIts, d2p3(:,1:nIts), 10, 'b'); hold on
%scatter(1:1770, permV, 1, 'k'); hold on
%scatter(1:1770, realV, 10, 'r')

exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% quantify variance with respect to shuffled variance for every item

nIts = size(neuralDistances, 2); 
h = zeros(nIts, 1); 
for itemi = 1:nIts

    varDistr = d2p2(:, itemi); 
    varObs = d2p1(itemi); 

    allAB = sum(varObs>varDistr);

    if allAB < 50
        h(itemi, :) = 1; 
    end
end

sum(h)



%% 
realV = var(neuralDistances,[],  1, 'omitnan');
permV = squeeze(var(neuralDistancesPerm,[],  2, 'omitnan'));

f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
subj_ch_fr = 7; 
subji = 1;
[ACT_BLNET] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
ACT_BLNET = squeeze(ACT_BLNET); 
%ACT_BLNET(logical(eye(size(ACT_BLNET, 1)))) = nan; 
%ACT_BLNET = mat2gray(ACT_BLNET); 

BLNETDistances = vectorizeRDM(ACT_BLNET)'; 



%% plot one subject one row

%d2p1 = realV; 
%d2p2 = permV; 

d2p1 = var(neuralDistances,[],  1, 'omitnan'); 
d2p2 = squeeze(mean(neuralDistancesPerm, 2, 'omitnan'));
d2p3 = squeeze(mean(BLNETDistances, 1, 'omitnan'));

cols = repelem([0.5, 0.5, 0.5], nPerm, 1); 
figure()
nIts = 100; 
scatter(1:nIts, d2p1(:,1:nIts), 10, 'k'); hold on
scatter(1:nIts, d2p2(:,1:nIts), 1, 'r')
scatter(1:nIts, d2p3(:,1:nIts), 10, 'b'); hold on
%scatter(1:1770, permV, 1, 'k'); hold on
%scatter(1:1770, realV, 10, 'r')

exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%[rho pVal] = corr(d2p1', d2p3', 'type', 's')
diff = abs(d2p1-d2p3); 
[h p ci ts] = ttest(diff)


%%

mRV = mean(realV); 
mPV = mean(permV, 2); 

histogram(mPV); hold on
scatter(mRV, 0, 2550, '.')
exportgraphics(gcf, 'myP.png', 'Resolution', 300);


%% 
realV = var(neuralCatDist,[],  1, 'omitnan');
permV = squeeze(var(neuralCatDistPerm,[],  2, 'omitnan'));

%%

mRV = mean(realV); 
mPV = mean(permV, 2); 

histogram(mPV); hold on
scatter(mRV, 0, 2550, '.')
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



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


%% compute 15*60 vector


for subji = 1:16
    
    rdmH = squeeze(meanRDM{subji, 1}); 
    ind = floor(meanRDM{subji, 2}/100); 
    d2p = 1- rdmH;
    [rdmMDS] = cmdscale(d2p);
    
    %compute centroids and distances to the centers
    clear cX
    nDims = size(rdmMDS, 2); 
    for i = 1:nDims
        cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
    end
    
    %%compute distance 2 centroid
    clear dist2Centr
    X = rdmMDS(:,1:nDims);
    nTrials = length(ind); 
    for itemi = 1:nTrials
        itCoord = X(itemi, : ); 
        centrCoord = cX(ind(itemi),:); 
        dist2Centr(:, itemi) =  sqrt (sum ( (itCoord - centrCoord) .^2 , 2, 'omitnan'));
    end
    
    %compute distance 2 centroid and to other category centroids
    clear allCat iO
    allCat = repmat(1:6, nTrials, 1);
    for triali = 1:nTrials
       h = allCat(triali,:) == ind(triali) ;
       iO(triali, :) = allCat(triali,~h) ;
    end
    iO = [ind iO];
    iO = iO(:);
    
    clear i2O i2oR
    for i = 1:nDims
       i2O(:,i) = cX(iO, i); 
    end
    X = rdmMDS(:,1:nDims);
    i2O = reshape (i2O, size(X, 1),[],  nDims);
    
    clear dis2c dis2
    for i = 1:6
        i2Otmp = squeeze(i2O(:,i,:));
        dis2c(:,i,:) =  sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
    end
    
    dist2catCen = dis2c(:,1);
    dist2otherCen = squeeze(mean(dis2c(:,2:6),2));
    
    sameMDiff = dist2catCen - dist2otherCen; 
    

    %%compute distance among category centroids
    D_cat_cent = pdist(cX, 'euclidean');
    all2allDist = pdist(X, 'euclidean');

    cX_all{subji} = cX; 
    dist2Centr_all{subji, :} = dist2Centr; 
    dist2OthCentr_all{subji, :} = dist2otherCen; 
    D_cat_cent_all(subji, :) = D_cat_cent; 
    sameMDiffAll{subji, :} = sameMDiff;
    all2allDistALL{subji, :} = all2allDist;
    

end

%% organize dist2Centr_all in a 60*1 line 

for subji = 1:16

    ids = meanRDM{subji, 2}; 
    idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
    [x1 x2 x3] = intersect (ids, idsF2);

    all2allCorr = zeros(1, 60); 
    dist2C(x3) = dist2Centr_all{subji}; 
    dist2C(dist2C==0) = nan; 
    dist2CALL(subji, :) = dist2C; 
    
    sdiffC = zeros(1, 60); 
    sdiffC(x3) = sameMDiffAll{subji}; 
    sdiffC(sdiffC==0) = nan; 
    sdiffCALL(subji, :) = sdiffC; 

end



%% plot dist2Centr_all
mDist = squeeze(mean(dist2CALL, 'omitnan'))
stdD = std(dist2CALL, 'omitnan'); 
se1 = stdD/sqrt(size(dist2CALL,1));

figure(); set(gcf, 'Position', [100 100 1000 300])
shadedErrorBar(1:60, mDist,se1, 'k', 1); hold on; 
errorbar(mDist, se1, 'k', LineWidth=1)
set(gca, 'xlim', [0 61])


exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% plot same minus
mDist = squeeze(mean(sdiffCALL, 'omitnan'))
stdD = std(sdiffCALL, 'omitnan'); 
se1 = stdD/sqrt(size(sdiffCALL,1));

figure(); set(gcf, 'Position', [100 100 1000 300])
shadedErrorBar(1:60, mDist,se1, 'k', 1); hold on; 
errorbar(mDist, se1, 'k', LineWidth=1)
set(gca, 'xlim', [0 61])

exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% permutations


nPerm = 100; 
clear dist2CSP
for permi = 1:nPerm

    for subji = 1:16
       
        dist2CSR = dist2CALL(subji, :); 
        ids4perm = randperm(60);
        dist2CSPerm = dist2CSR(ids4perm);
        dist2CSP(permi, subji,:) = dist2CSPerm; 

    end

end

%%  plot with respect to distributino

figure(); set(gcf, 'Position', [100 100 1000 300])
d2p = squeeze(mean(dist2CSP, 2, 'omitnan'));  hold on; 
plot(d2p', 'k')
plot(mDist, 'r', LineWidth=3);

exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% plot D_cat_cent 

dCatCentr = mean(D_cat_cent_all); 
stdD = std(D_cat_cent_all, 'omitnan'); 
se1 = stdD/sqrt(size(D_cat_cent_all,1));
shadedErrorBar(1:15, dCatCentr,se1, 'k', 1); hold on; 
set(gca, 'xlim', [0 16])
errorbar(dCatCentr, se1, 'k', LineWidth=1)

exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Compute distances for BL-NET

subj_ch_fr = 7;

for subji = 1:16

    clearvars -except meanRDM subji ACT nDims mDist dist2Centr_all_BLNET dist2Centr_all cX_all D_cat_cent_all_BLNET D_cat_cent_all sameMDiffAll_BLNET ...
                        sameMDiffAll all2allDistALL_BLNET all2allDistALL dist2OthCentr_all dist2OthCentr_all_BLNET
    ids = meanRDM{subji, 2}; 
    ids2 = char(string(ids)); 
    for i = 1:length(ids2)
        idT = ids2(i,:); 
        ids3(i,:) = double(string(idT([1 3])));
    end
    idx = find(~mod(ids3, 10)); 
    ids4 = ids3-10; ids4(idx) = ids3(idx);
    ind = floor(meanRDM{subji, 2}/100); 



    f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    subj_ch_fr = subji; 
    [ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
    ACT = squeeze(ACT); 
       
    ACTH = ACT(ids4, ids4); 
    d2p = 1- ACTH;
    [rdmMDS] = cmdscale(d2p);
    
    %compute centroids and distances to the centers
    clear cX
    nDims = size(rdmMDS, 2); 
    for i = 1:nDims
        cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
    end
    %%compute distance among category centroids
    D_cat_centACT = pdist(cX); 
    
    %%compute distance 2 centroid
    clear dist2CentrACT
    X = rdmMDS(:,1:nDims);
    nTrials = length(ind); 
    for itemi = 1:nTrials
        itCoord = X(itemi, : ); 
        centrCoord = cX(ind(itemi),:); 
        dist2CentrACT(:, itemi) =  sqrt (sum ( (itCoord - centrCoord) .^2 , 2, 'omitnan'));
    end
    

    %compute distance 2 centroid and to other category centroids
    clear allCat iO
    allCat = repmat(1:6, nTrials, 1);
    for triali = 1:nTrials
       h = allCat(triali,:) == ind(triali) ;
       iO(triali, :) = allCat(triali,~h) ;
    end
    iO = [ind iO];
    iO = iO(:);
    
    clear i2O i2oR
    for i = 1:nDims
       i2O(:,i) = cX(iO, i); 
    end
    X = rdmMDS(:,1:nDims);
    i2O = reshape (i2O, size(X, 1),[],  nDims);
    
    clear dis2c dis2
    for i = 1:6
        i2Otmp = squeeze(i2O(:,i,:));
        dis2c(:,i,:) =  sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
    end
    
    dist2catCen = dis2c(:,1);
    dist2otherCen = squeeze(mean(dis2c(:,2:6),2));
    
    
    sameMDiffBLNET = dist2catCen - dist2otherCen; 

    all2allDist_BLNET = pdist(X);
    

    cX_all_BLNET{subji} = cX; 
    dist2Centr_all_BLNET{subji, :} = dist2CentrACT; 
    dist2OthCentr_all_BLNET{subji, :} = dist2otherCen; 
    D_cat_cent_all_BLNET(subji, :) = D_cat_centACT; 

    sameMDiffAll_BLNET{subji, :} = sameMDiffBLNET;
    all2allDistALL_BLNET{subji, :} = all2allDist_BLNET;


end

%%

for subji =1:16

    x = dist2Centr_all{subji}'; 
    y = dist2Centr_all_BLNET{subji}'; 
    [rhoAll1(subji,:) pAll(subji, :)] = corr(x, y, 'type', 's'); 

    x = D_cat_cent_all(subji,:)'; 
    y = D_cat_cent_all_BLNET(subji,:)'; 
    [rhoAll2(subji,:) pAll2(subji, :)] = corr(x, y, 'type', 's'); 

    x = sameMDiffAll{subji}; 
    y = sameMDiffAll_BLNET{subji}; 
    [rhoAll3(subji,:) pAll3(subji, :)] = corr(x, y, 'type', 's'); 

    x = all2allDistALL{subji}'; 
    y = all2allDistALL_BLNET{subji}'; 
    [rhoAll4(subji,:) pAll4(subji, :)] = corr(x, y, 'type', 's'); 

    x = dist2OthCentr_all{subji}; 
    y = dist2OthCentr_all_BLNET{subji}; 
    [rhoAll5(subji,:) pAll5(subji, :)] = corr(x, y, 'type', 's'); 

    % figure()
    % plot(x, LineWidth=2); hold on
    % plot(y, LineWidth=2)
    % set(gca, 'xlim', [0 61], 'ylim', [.1 1.1])
    % title([num2str(rhoAll(subji,:)) ' // ' num2str(pAll(subji, :))])
    % exportgraphics(gcf, ['allM' num2str(subji, '%02.f') '.png'], 'Resolution', 300);
    % close all

end

rhoAll1(1) = []; 
rhoAll2(1) = []; 
rhoAll3(1) = []; 
rhoAll4(1) = []; 
[h p ci ts] = ttest(rhoAll1);
t1 = ts.tstat; 
[h p ci ts] = ttest(rhoAll2);
t2 = ts.tstat; 
[h p ci ts] = ttest(rhoAll3);
t3 = ts.tstat; 
[h p ci ts] = ttest(rhoAll4);
t4 = ts.tstat; 
[h p ci ts] = ttest(rhoAll5);
t4 = ts.tstat; 

%% plot one bar
clear data
data.data = [rhoAll4]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);



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

nPerm = 100;
witBet = 0; 
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
percHs = sum (h1)


%% Analysis at different noise levels WITH RDMS, no MDS

noiseLevels = [0 .1 .5 1 3 5 10 20 30 50]
%noiseLevels = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.2 1.4 1.6 1.8 2]

for noisei = 1:length(noiseLevels)

    mu = noiseLevels(noisei); 
    ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
    
    CM = kron(eye(6), ones(10));
    CM = CM + normrnd(0, mu,[60,60]);
    CM = normalize(CM, 'range', [-1, 1]);
    %CM = mat2gray(CM);
    CM(logical(eye(size(CM)))) = 1; 
    CM = tril(CM,-1)+tril(CM)'; 
    CCI = compute_CCI_WM(CM, ind);
    
    figure()
    tiledlayout(1, 4)
    nexttile
    imagesc(CM); ; axis square ;
    d2p = 1- CM;
    [rdmMDS] = cmdscale(d2p);
    
    %Encoding
    CCI = compute_CCI_WM(CM, ind); 
    m1_CCI_S = CCI; 
    m2_WCSM_S= averageCAT_WM(CM, ind);
    m3_WCSV_S= var(averageCAT_WM(CM, ind));
    m4_WISM_S= vectorizeRDM(CM); 
    m5_WISV_S= var(vectorizeRDM_WM(CM), 'omitnan'); 

    
    cols = brewermap(6, 'Accent'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
    cols = repelem(cols, 10, 1);
    nexttile
    scatter(rdmMDS(:,1),rdmMDS(:,2),550, cols, '.'); axis square
    set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
    nexttile
    plot(m2_WCSM_S); axis square
    %bar(m2_WCSM_S); axis square
    nexttile
    plot(m4_WISM_S); axis square
    %bar(m3_WCSV_S); axis square
    set(gca, 'ylim', [-0.00 0.003])

    
    sgtitle(['Noise: ' num2str(mu) ' CCI: ' num2str(CCI)])

    exportgraphics(gcf, ['allM' num2str(noisei, '%02.f') '.png'], 'Resolution', 300);
    close all


end

%% Produce 100 equidistant

nPerm = 500; 
noiseLevels = [0 .1 .5 1 3 5 10 20 30 50]

clear m2_WCSM_S
for permi = 1:nPerm
    for noisei = 2 % 1:length(noiseLevels)
    
        mu = noiseLevels(noisei); 
        ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
        
        CM = kron(eye(6), ones(10));
        CM = CM + normrnd(0, mu,[60,60]);
        CM = normalize(CM, 'range', [-1, 1]);
        %CM = mat2gray(CM);
        CM(logical(eye(size(CM)))) = 1; 
        CM = tril(CM,-1)+tril(CM)'; 
        CCI = compute_CCI_WM(CM, ind);
        
        m2_WCSM_S(permi, :)= averageCAT_WM(CM, ind);
        m4_WISM_S(permi, :)= vectorizeRDM(CM);

    end
end



%%
N = 500; 
a=-1; % lower boundary
b=1; % higher boundary
nbins = 6; % number of bin
x = m4_WISM_S(:,2); 
histogram(x)
edges = linspace(a,b,nbins+1); % edges of the bins
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(x,'Expected',E,'Edges',edges)

%%
N=1000; % sample size
a=-1; % lower boundary
b=1; % higher boundary
nbins = 10; % number of bin

x=unifrnd(a,b,N,1);
%x(x<.9) = rand(sum(x<.9),1);
histogram(x, nbins)


edges = linspace(a,b,nbins+1); % edges of the bins
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(x,'Expected',E,'Edges',edges)



%%
N=6; % sample size
a=-1; % lower boundary
b=1; % higher boundary
nbins = 6; % number of bin

x=unifrnd(a,b,N,1);
%x = [.2 .2 .2 .2 .2 .2]

edges = linspace(a,b,nbins+1); % edges of the bins
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(x,'Expected',E,'Edges',edges)

%%
h = histogram(x, edges)
chi = sum((h.Values - N/nbins).^2 / (N/nbins));
k = nbins-1; % degree of freedom
z = chi2cdf(chi, k)
p2 = chi2cdf(chi,k,'upper')

%% Test how normalization affects 

CM = kron(eye(6), ones(10));
CM = CM + normrnd(0, 1,[60,60]);
CMN = normalize(CM, 'range');
%CM(logical(eye(size(CM)))) = 1; 
CM = tril(CM,-1)+tril(CM)'; 
%CMN(logical(eye(size(CM)))) = 1; 
CMN = tril(CMN,-1)+tril(CMN)'; 

figure
tiledlayout(1, 2)
nexttile
imagesc(CM); axis square; colorbar
nexttile
imagesc(CMN); axis square; colorbar

CM1 = vectorizeRDM(CM); 
CMN1 = vectorizeRDM(CMN); 
CMod = kron(eye(6), ones(10));
CMod = vectorizeRDM(CMod); 

rho1 = corr(CM1, CMod, 'type', 's'); 
rho2 = corr(CMN1, CMod, 'type', 's'); 



%% compare real data against noise 

nPerm = 1000; 
noiseLevel = 50; 


clear ids5
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

    mu = noiseLevel; 
    

    for permi = 1:nPerm
        CM = kron(eye(6), ones(10));
        CM = CM + normrnd(0, mu,[60,60]);
        CM = normalize(CM, 'range');
        CM(logical(eye(size(CM)))) = 1; 
        CM = tril(CM,-1)+tril(CM)'; 
        
        CM = CM(ids4, ids4); 
        d2p = 1- CM;
        [rdmMDS] = cmdscale(d2p);
        
        %compute centroids and distances to the centers
        clear cX
        nDims = size(rdmMDS, 2); 
        for i = 1:nDims
            cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
        end
        %D_cat_centACT = pdist(cX); 
    
        all2allDist = pdist(rdmMDS);

        %%compute distance 2 centroid
        clear dist2CentrACTPerm
        X = rdmMDS(:,1:nDims);
        
        for itemi = 1:nTrials
            itCoord = X(itemi, : ); 
            centrCoord = cX(ind(itemi),:); 
            dist2CentrACTPerm(:, itemi) =  sqrt (sum ( (itCoord - centrCoord) .^2 , 2, 'omitnan'));
        end
    
        
        %realDataSubj = dist2Centr_all{subji}; 
        realDataSubj = all2allDistALL{subji}; 


        %rhoPerm(subji, permi, :) = corr(realDataSubj', dist2CentrACTPerm', 'type', 's'); 
        rhoPerm(subji, permi, :) = corr(realDataSubj', all2allDist', 'type', 's'); 

        

    end


end



[hPerm pPerm ci ts] = ttest(rhoPerm); 
tPerm = ts.tstat; 


%%


histogram(tPerm); hold on
scatter(t4, 0, 2550, '.')
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Analysis at different noise levels (CORRELATION AND NOT DISTANCE)

noiseLevels = [0 .1 .5 1 3 5 10 20 30 50];

for noisei = 1:length(noiseLevels)

    mu = noiseLevels(noisei); 

    CM = kron(eye(6), ones(10));
    CM = CM + normrnd(0, mu,[60,60]);
    CM = normalize(CM, 'range', [0 1]);
    %CM = 2 * mat2gray(CM) - 1;

    CM(logical(eye(size(CM)))) = nan; 
    CM = tril(CM,-1)+tril(CM)'; 
    CCI = compute_CCI_WM(CM);
    
    figure()
    tiledlayout(1, 2)
    nexttile
    imagesc(CM); axis square ;
    
    nexttile
    d2p = vectorizeRDM(CM); 
    plot(d2p); axis square ;
    
    sgtitle(['Noise: ' num2str(mu) ' CCI: ' num2str(CCI)])

    exportgraphics(gcf, ['allM' num2str(noisei, '%02.f') '.png'], 'Resolution', 300);
    close all


end



%% EVERY SUBJECT DURING ENCODING AND MAINTENANCE
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
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        
        %rdm = nan(60); 
        %rdm(x3,x3) = rdmS; 
        %meanRDM_ENC(subji,:,:) = rdmS; 

        meanRDM_ENC{subji,1}= rdmS; 
        meanRDM_ENC{subji,2}= idsF1; 
        
    end
end

clearvars -except meanRDM_ENC 

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
        
        %rdm = nan(60); 
        %rdm(x3,x3) = rdmS; 
        %meanRDM_MAINT(subji,:,:) = rdmS; 

        meanRDM_MAINT{subji,1}= rdmS; 
        meanRDM_MAINT{subji,2}= idsF1; 
        
    end
end




%% Compute metrics (no plotting, no MDS)

for subji = 1:16
    
    %Encoding
    mRDM_E = squeeze(meanRDM_ENC{subji,1}); 
    ind_E = floor(meanRDM_ENC{subji, 2}/100);     
    ids_E = meanRDM_ENC{subji, 2}; 
    CME = generate_catModel_WM(ids_E); 

    CCI = compute_CCI_WM(mRDM_E, ind_E); 
    m1_CCI_E(subji, :) = CCI; 
    m2_WCSM_E(subji, :) = mean(averageCAT_WM(mRDM_E, ind_E));
    m3_WCSV_E(subji, :) = var(averageCAT_WM(mRDM_E, ind_E));
    m4_WISM_E(subji, :) = mean(vectorizeRDM(mRDM_E), 'omitnan'); 
    m5_WISV_E(subji, :) = var(vectorizeRDM_WM(mRDM_E), 'omitnan'); 
    m6_WISMWC_E(subji, :) = mean(mRDM_E(CME==1), 'omitnan'); 
    m7_WISMBC_E(subji, :) = mean(mRDM_E(CME==0), 'omitnan'); 
    m8_WISVWC_E(subji, :) = var(mRDM_E(CME==1), 'omitnan'); 
    m9_WISVBC_E(subji, :) = var(mRDM_E(CME==0), 'omitnan'); 

    % Maintenance
    mRDM_M = squeeze(meanRDM_MAINT{subji,1}); 
    ind_M = floor(meanRDM_MAINT{subji, 2}/100); 
    ids_M = meanRDM_MAINT{subji, 2}; 
    CMM = generate_catModel_WM(ids_M); 

    CCI = compute_CCI_WM(mRDM_M, ind_M); 
    m1_CCI_M(subji, :) = CCI; 
    m2_WCSM_M(subji, :) = mean(averageCAT_WM(mRDM_M, ind_M));
    m3_WCSV_M(subji, :) = var(averageCAT_WM(mRDM_M, ind_M));
    m4_WISM_M(subji, :) = mean(vectorizeRDM(mRDM_M), 'omitnan'); 
    m5_WISV_M(subji, :) = var(vectorizeRDM_WM(mRDM_M), 'omitnan'); 
    m6_WISMWC_M(subji, :) = mean(mRDM_M(CMM==1), 'omitnan'); 
    m7_WISMBC_M(subji, :) = mean(mRDM_M(CMM==0), 'omitnan'); 
    m8_WISVWC_M(subji, :) = var(mRDM_M(CMM==1), 'omitnan'); 
    m9_WISVBC_M(subji, :) = var(mRDM_M(CMM==0), 'omitnan'); 
    
end




%% plot two bar
clear data

%data.data = [m2_WCSM_E m2_WCSM_M];
%data.data = [m3_WCSV_E m3_WCSV_M];
%data.data = [m4_WISM_E m4_WISM_M];
%data.data = [m5_WISV_E m5_WISV_M];
%data.data = [m6_WISMWC_E, m6_WISMWC_M]; 
%data.data = [m7_WISMBC_E, m7_WISMBC_M]; 
data.data = [m8_WISVWC_E, m8_WISVWC_M]; 
%data.data = [m9_WISVBC_E, m9_WISVBC_M]; 


sub2exc = 1; 
data.data(sub2exc, :) = []; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2], data.data, 50, 'k'); hold on;
hb = plot(data.data'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 3] );
%set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
box on; 
[h p ci t] = ttest (data.data(:,1), data.data(:,2));
res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
disp (res2title);
title(res2title)



exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Characterize Category model
catDistCM = kron(eye(6), ones(10)); 
catDistCM = catDistCM(tril(true(length(catDistCM)), -1));
catDistVar = var(catDistCM)
catDistMean = mean(catDistCM)


f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
subj_ch_fr = subji; 
[ACT] = load_BLNETi(cfg, subji, subj_ch_fr, paths);%load network if not loaded yet
ACT = squeeze(ACT); 
BLNETDistances = vectorizeRDM(ACT); 
blnetDistVar = var(BLNETDistances)
blnetDistMean = mean(BLNETDistances)


%% CORRELATE RDMS BETWEEN ENCODING AND MAINTENANCE (SORT THEM FIRST) 

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
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        
        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        meanRDM_ENC(subji,:,:) = rdm; 

        
    end
end

clearvars -except meanRDM_ENC 

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
        meanRDM_MAINT(subji,:,:) = rdm; 

        
    end
end

%% 

clear rhoALL
for subji = 1:16

    rdmE = squeeze(meanRDM_ENC(subji, :, :)); 
    %rdmE = mat2gray(rdmE); 
    
    rdmM = squeeze(meanRDM_MAINT(subji, :, :)); 
    %rdmM = mat2gray(rdmM); 

    rdmEV = vectorizeRDM(rdmE); 
    rdmMV = vectorizeRDM(rdmM); 

    nanIDs = isnan(rdmMV); 
    rdmEV(nanIDs) = []; 
    rdmMV(nanIDs) = []; 

    [rhoALL(subji,:)] = corr(rdmEV', rdmMV', 'type', 's'); 

end


%% plot one bar
clear data
data.data = [rhoALL]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);










%% plot MDS
d2p = 1- squeeze(mRDM);
[rdmMDS] = cmdscale(d2p);
nDims = 2; 

cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);
scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
 

%% compute centroids and distances to the centers
%centroid coordinates for each category

ind = [ones(1, 10) ones(1, 10) *2 ones(1, 10) *3 ones(1, 10) *4 ones(1, 10) *5 ones(1, 10) *6]';
clear cX
for i = 1:nDims
    cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
end

%% plot MDS + centroids 
scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); hold on; 
scatter(cX(:,1),cX(:,2),1350, cols([1 11 21 31 41 51], :), '.'); axis square
set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
exportgraphics(gcf, 'myP.png', 'Resolution', 300)

%% compute distance 2 centroid

X = rdmMDS(:,1:nDims);
cIds = repelem(1:6, 10);

for itemi = 1:60

    itCoord = X(itemi, : ); 
    centrCoord = cX(cIds(itemi),:); 
    dist2Centr(itemi, :) =  sqrt (sum ( (itCoord - centrCoord) .^2 , 2, 'omitnan'));

end


figure()
histogram(dist2Centr);
exportgraphics(gcf, 'myP.png', 'Resolution', 300)



%% compute distance 2 centroid and to other category centroids

nTrials = 60; 
clear allCat iO
allCat = repmat(1:6, nTrials, 1);
for triali = 1:nTrials
   h = allCat(triali,:) == ind(triali) ;
   iO(triali, :) = allCat(triali,~h) ;
end
iO = [ind iO];
iO = iO(:);

clear i2O i2oR
for i = 1:nDims
   i2O(:,i) = cX(iO, i); 
end
X = rdmMDS(:,1:nDims);
i2O = reshape (i2O, size(X, 1),[],  nDims);


clear dis2c dis2
for i = 1:6
    i2Otmp = squeeze(i2O(:,i,:));
    dis2c(:,i,:) =  sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
end

dist2catCen = dis2c(:,1);
md2_CorrCat = accumarray(ind, dist2catCen, [6 1], @(x)mean(x,'omitnan')); % all trials

dist2otherCen = squeeze(mean(dis2c(:,2:6),2));
md2_OtherCat = accumarray(ind, dist2otherCen, [6 1], @(x)mean(x,'omitnan')); % all trials

sameMDiff = dist2catCen - dist2otherCen; 




%% compute distance among category centroids

D_enc_c1 = pdist(cX_enc_c1, 'euclidean')
D_enc_c2 = pdist(cX_enc_c2, 'euclidean')

Z_enc_1 = squareform(D_enc_c1);
Z_enc_2 = squareform(D_enc_c2);
Z_prio = squareform(D_prio);

tiledlayout(1,3)
nexttile
imagesc(Z_enc_1); axis square; 
set(gca, 'clim', [0 .2])

nexttile
imagesc(Z_enc_2); axis square; 
set(gca, 'clim', [0 .2])

nexttile
imagesc(Z_prio); axis square
set(gca, 'clim', [0 .2])

exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% plot bars
clear data
data.data = [sameMDiff_enc_c1 sameMDiff_prio]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2], data.data, 'k'); hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
%set(gca, 'ylim', [0 .45])

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);
exportgraphics(gcf, 'allM.png', 'Resolution', 300);







%% Compute differences in fit between category model and BLNET for every subject

nPerm = 100; 
clear rhoALL rhoPerm rhoALLCAT neuralCatDist catDistBLNET rhoCATPerm catDistBLNETPERM
for subji = 1:16
    ids = meanRDM{subji, 2}; 
    
    ACT = squeeze(meanRDM{subji, 1}); 
    neuralDistances = vectorizeRDM(ACT); 
    ind = floor(ids/100); 
    
    
    f2sav = 'BLNETi_pfc_E123_[56]_3-54_0_0_0_1_.1_5_1.mat'; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    subj_ch_fr = subji; 
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
    BLNETDistances = vectorizeRDM(ACTH); 

    CMDistances = kron(eye(6), ones(10));
    CMDistances = CMDistances(ids4, ids4); 
    CMDistances = vectorizeRDM(CMDistances); 
    
    
    
    [rhoBLNET] = corr(neuralDistances, BLNETDistances, 'type', 's'); 
    [rhoCM] = corr(neuralDistances, CMDistances, 'type', 's'); 

    rhoALL(subji, :) = rhoBLNET - rhoCM; 

end

[h p ci ts] = ttest(rhoALL)




%% COMPUTE DIFFERENCES IN FITS FOR CATAEGORY MODEL DURING ENCODING AND MAINTENANCE 

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
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        
        %rdm = nan(60); 
        %rdm(x3,x3) = rdmS; 
        %meanRDM_ENC(subji,:,:) = rdm; 
        
        meanRDM_ENC{subji,:} = rdmS; 
        allCME{subji,:} = squeeze(load_CATMODEL_activ(ids));
        
        


        
    end
end

clearvars -except meanRDM_ENC allCME

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
        
        %rdm = nan(60); 
        %rdm(x3,x3) = rdmS; 
        %meanRDM_MAINT(subji,:,:) = rdm; 

        meanRDM_MAINT{subji,:} = rdmS; 
        allCMM{subji,:} = squeeze(load_CATMODEL_activ(ids));

        
    end
end


%% 

clear rhoE rhoM
for subji = 1:16

    %rdmE = squeeze(meanRDM_ENC(subji, :, :)); 
    %rdmM = squeeze(meanRDM_MAINT(subji, :, :)); 

    rdmE = meanRDM_ENC{subji}; 
    rdmM = meanRDM_MAINT{subji}; 
    
    rdmEV = vectorizeRDM(rdmE); 
    rdmMV = vectorizeRDM(rdmM); 

    cME = vectorizeRDM(allCME{subji}); 
    cMM = vectorizeRDM(allCMM{subji});


    nanIDs = isnan(rdmMV); 
    rdmEV(nanIDs) = []; 
    rdmMV(nanIDs) = []; 

    [rhoE(subji,:)] = corr(rdmEV', cME', 'type', 's'); 
    [rhoM(subji,:)] = corr(rdmMV', cMM', 'type', 's'); 

end

%% 
[h p ci ts] = ttest(rhoE)

diffD = rhoE-rhoM; 
[h p ci ts] = ttest(diffD)


%% plot one bar
clear data
data.data = [diffD]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 24, 'linew',1, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
box on

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 2);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);





%% COMPARE NUMBER OF SIGNIFICANT ITEMS WITH RESPECT TO CHANCE (SORT FIRST)

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
        idsF1 = oneListIds(:, 3);
        idsF2 = [101:110 201:210 301:310 401:410 501:510 601:610]';
        [x1 x2 x3] = intersect (idsF1, idsF2);
        
        nRows       = size(neuralRDMs, 1); 
        neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
        rdmS         = mean(neuralRDM1(:, :,idsClust ), 3); 
        
        rdm = nan(60); 
        rdm(x3,x3) = rdmS; 
        %rdm= mat2gray(rdm); 
        %rdm= normalize(rdm, 'range', [-1 1]); 
        mRDME(subji,:,:) = rdm; 
        
        
    end
end

clearvars -except mRDME 

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
        %rdm= mat2gray(rdm); 
        %rdm= normalize(rdm, 'range', [-1 1]);  
        mRDMM(subji,:,:) = rdm; 

    end
end

sub2exc = 1; 
mRDME(sub2exc, :, :) = []; 
mRDMM(sub2exc, :, :) = []; 




%%

clear mEH mMH mEHP d4TT allAB
for subji = 1:15
    mEH(subji, :) = vectorizeRDM(squeeze(mRDME(subji, :, :))); 
    mMH(subji, :) = vectorizeRDM(squeeze(mRDMM(subji, :, :))); 

end

[h p ci ts] = ttest(mEH); 
hE = squeeze(h); tE = squeeze(ts.tstat); 
[h p ci ts] = ttest(mMH); 
hM = squeeze(h); tM = squeeze(ts.tstat); 

tE2 = tE; 
tE2(hE==0) = nan; 

figure(); 
nIts = 500; 
scatter(1:nIts, tE(:,1:nIts), 10, 'k'); hold on
scatter(1:nIts, tE2(:,1:nIts), 300, 'r.')
sum(hE)
sum(hM)



%%
nPerm = 100; 

clear mEH
for subji = 1:15
    mRDMEH = squeeze(mRDME(subji, :, :)); 
    mRDMMH = squeeze(mRDMM(subji, :, :)); 
    for permi = 1:nPerm
        ids4perm = randperm(60); 
        mEH(subji, permi, :) = vectorizeRDM(mRDMEH(ids4perm, ids4perm)); 
        mEM(subji, permi, :) = vectorizeRDM(mRDMMH(ids4perm, ids4perm)); 
    end
end


[h p ci ts] = ttest(mEH); 
hEP = squeeze(h); tEP = squeeze(ts.tstat); 
totalHEP = sum(hEP, 2); 

[h p ci ts] = ttest(mEM); 
hMP = squeeze(h); tMP = squeeze(ts.tstat); 
totalMEP = sum(hMP, 2); 















%% plot mean RDM for each subject during ENC and MAINT

for subji = 1:16

    mRDM_E = squeeze(meanRDM_ENC{subji,1}); 
    mRDM_M = squeeze(meanRDM_MAINT{subji,1}); 

    d2pE = mRDM_E; 
    d2pM = mRDM_M; 
    
    tiledlayout(1,4)
    nexttile
    imagesc(d2pE); axis square
    set(gca, 'clim', [-.25 .25])

    nexttile
    imagesc(d2pM); axis square
    set(gca, 'clim', [-.25 .25])

    nexttile
    ind = floor(meanRDM_ENC{subji, 2}/100); 
    %mRDM_E1 = mat2gray(mRDM_E); 
    %mRDM_E = mat2gray(mRDM_E); 
    d2pE = 1- squeeze(mRDM_E);
    [rdmMDS] = cmdscale(d2pE);
    cols = brewermap(6, 'Accent'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
    cols = cols(ind, :);
    [id1 id2] = sort(cols(:,1)); 
    cols = cols (id2,:); 


    nDims = size(rdmMDS, 2); 
    %nDims = 2;     
    clear cX 
    for i = 1:nDims
        cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
    end
    CCI = compute_CCI_WM(mRDM_E, ind); 
    dist2C3E(subji, :) = CCI; 
    dist2CE(subji, :) = var(pdist(cX, 'euclidean')); 
    dist2C2E(subji, :) = mean(pdist(rdmMDS(:, 1:nDims), 'euclidean')); 
    %mRDM_E1 = normalize(mRDM_E, 'range'); %%THIS DOES NOT WORK SMALL VALUES GET ZEROS
    %mRDM_E1 = mat2gray(mRDM_E); 
    %mRDM_E1 = mRDM_E;
    
    wCatSimE(subji, :) = mean(averageCAT_WM(mRDM_E, ind));
    dist2C4E(subji, :) = var(vectorizeRDM(mRDM_E), 'omitnan'); 

    [idu1 idu2 idu3] = unique(cols, 'First', 'rows'); 
    scatter(rdmMDS(:,1),rdmMDS(:,2),50, cols, '.'); hold on; 
    scatter(cX(:,1),cX(:,2),350, cols(idu2, :), '.'); axis square
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1])
    %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
    title(num2str(CCI));

   
    nexttile
    ind = floor(meanRDM_MAINT{subji, 2}/100); 
    %mRDM_M1 = mat2gray(mRDM_M); 
    d2pM = 1- squeeze(mRDM_M);
    [rdmMDS] = cmdscale(d2pM);
    cols = brewermap(6, 'Accent'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
    cols = cols(ind, :);
    [id1 id2] = sort(cols(:,1)); 
    cols = cols (id2,:); 
    
    
    nDims = size(rdmMDS, 2); 
    %nDims = 2; 
    clear cX 
    for i = 1:nDims
        cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
    end
    CCI = compute_CCI_WM(mRDM_M, ind); 
    dist2C3M(subji, :) = CCI; 
    dist2CM(subji, :) = var(pdist(cX, 'euclidean')); 
    dist2C2M(subji, :) = mean(pdist(rdmMDS(:, 1:nDims), 'euclidean')); 
    %mRDM_M1 = normalize(mRDM_M, 'range'); %%THIS DOES NOT WORK SMALL VALUES GET ZEROS
    %mRDM_M1 = mat2gray(mRDM_M); 
    %xT = corr(vectorizeRDM(mRDM_M), vectorizeRDM(mRDM_M1))     % JUST TO CHECK THE CORRELATION IS 1

    wCatSimM(subji, :) = mean(averageCAT_WM(mRDM_M1, ind));
    dist2C4M(subji, :) = var(vectorizeRDM(mRDM_M1), 'omitnan'); 

    [idu1 idu2 idu3] = unique(cols, 'First', 'rows'); 
    scatter(rdmMDS(:,1),rdmMDS(:,2),50, cols, '.'); hold on; 
    scatter(cX(:,1),cX(:,2),350, cols(idu2, :), '.'); axis square
    %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1])
    title(num2str(CCI));
    
    exportgraphics(gcf, ['myP' num2str(subji, '%02.f') '.png'], 'Resolution', 300)

end
%CCI = compute_CCI_WM(mRDM)


%% compare real data against noise (CORRELATION AND NOT DISTANCE)
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


%% compare real data against noise (CORRELATION AND NOT DISTANCE)

nPerm = 100; 
noiseLevels = [0 .1 1 3 5 10];


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

[h p ci ts] = ttest(neuralBLNETCORR); 
t4 = ts.tstat; 

[hPerm pPerm ci ts] = ttest(rhoPermALL); 
tPermALL = squeeze(ts.tstat); pPermALL = squeeze(pPerm); 

[hPerm pPerm ci ts] = ttest(rhoPermALLBLNET); 
tPermALLBLNET = squeeze(ts.tstat); 

%% count % of significant p's

clear pH pH1
for noisei = 1:6
    pH = pPermALL(:, :, noisei);
    pH1(:, noisei) = sum(pH<0.05, 2);
end


%%

figure(); set(gcf, 'Position', [100 100 500 500])
tiledlayout(6, 2)
for tilei = 1:6
    if tilei == 1
        nexttile
        scatter(tPermALL(1, 1), 0, 550, '.'); hold on; 
        scatter(t4, 0, 550, '.', 'r'); 
        nexttile
        scatter(tPermALLBLNET(1,1), 0, 550, '.'); hold on; 
        scatter(t4, 0, 550, '.', 'r'); 
    else
        nexttile
        histogram(tPermALL(:, tilei)); hold on
        scatter(t4, 0, 550, '.', 'r'); 
        nexttile
        histogram(tPermALLBLNET(:, tilei)); hold on
        scatter(t4, 0, 550, '.', 'r'); 
    end
    
end
exportgraphics(gcf, 'allM.png', 'Resolution', 300);






%%



