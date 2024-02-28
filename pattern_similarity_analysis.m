%% load cfg_contrasts
%% 
clearvars
region              = 'vvs';
paths = load_paths_WM(region, 'none');
filelistSess = getFilesWM(paths.out_contrasts);


frequncies2test = [{3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = { '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 0; 
meanInFreq          = 0; 
avMeth              = 'none';  % average across image repetitions or not
TG                  = 1; %temporal generalization
contr2save          = {'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3' 'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE'  }; %
%contr2save          = {'SISC_EE' 'DISC_EE' 'DIDC_EE' 'SISC_EM2' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'SCSP_M2M2' 'SCDP_M2M2'}; %{};
%contr2save          = {'DISC_EM1' 'DIDC_EM1'}; %{};
%contr2save          = {'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3'}; %{};
bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
n2s                 = 10000;
%n2s                 = 5000;
loadSurr            = 0; 
zScType             = 'sess'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
aVTime              = 0; % Average or not in time the feature vectors
 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 7:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.out_contrasts filelistSess{sessi}]);   
    
    disp ([ 'fnames = ' fnames{:} newline 'win_width = ' num2str(win_width) newline 'mf = ' num2str(mf) newline ...
            'meanInTime = ' num2str(meanInTime) newline 'meanInFreq = ' num2str(meanInFreq) newline 'TG = ' num2str(TG) newline ...
            'bline = ' num2str(bline) newline 'acrossTrials = ' num2str(acrossTrials) newline 'batch_bin = ' num2str(batch_bin)  ]);
        
    
    cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline);
% %     % % check the normalization across trials works 
% %     ids = str2num(cell2mat(cfg_contrasts.oneListIds));
% %     d2p = squeeze(mean(cfg_contrasts.oneListPow(ids(:, 1) ==1, :, :, :)));
% %     imagesc(squeeze(d2p(1,:,:))); colorbar % % % because of the normalization across trials

 
    cfg_contrasts.contr2save = contr2save';
    cfg_contrasts.n2s = n2s;
    cfg_contrasts.loadSurr = loadSurr;
    cfg_contrasts.batch_bin = batch_bin;
    

    if strcmp(avMeth,'pow')
        cfg.avRep = 1; 
        cfg_contrasts = average_repetitions(cfg, cfg_contrasts);
    end
    
    [out_contrasts] = create_contrasts_WM (cfg_contrasts);
        
 
    for freqi = 1:length(frequncies2test)
        fname = fnames{freqi};
        mkdir ([paths.results.bands fname]);
        cd ([paths.results.bands fname]);
        f           = frequncies2test{freqi}; 
        % % the rsa_WM function saves the similarity matrices (TG=1) or diagonals (TG = 0)
        rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, sessi, TG, aVTime)
        cd ..
    end
 
end

cd .. 

 
%%process Folders and create one file per condition including all subjects

clearvars -except region
paths = load_paths_WM(region, 'none');
currentDir = pwd; 
mkdir(paths.results.bands)
cd (paths.results.bands)
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
 
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym2use frequncies2test currentDir
    
    direct = fold(foldi);
    cd (direct.name)
 
    processFoldersWM; 

    cd .. 
end

cd (currentDir);
disp ('done')

% % % % % % 

clearvars
region              = 'pfc';
paths = load_paths_WM(region, 'none');
filelistSess = getFilesWM(paths.out_contrasts);


frequncies2test = [{3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = { '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 0; 
meanInFreq          = 0; 
avMeth              = 'none';  % average across image repetitions or not
TG                  = 1; %temporal generalization
contr2save          = { 'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3' 'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE' 'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2'}; %
%contr2save          = {'SISC_EE' 'DISC_EE' 'DIDC_EE' 'SISC_EM2' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'SCSP_M2M2' 'SCDP_M2M2'}; %{};
%contr2save          = {'DISC_EM1' 'DIDC_EM1'}; %{};
%contr2save          = {'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3'}; %{};
bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
n2s                 = 10000;
%n2s                 = 5000;
loadSurr            = 0; 
zScType             = 'sess'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
aVTime              = 0; % Average or not in time the feature vectors
 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.out_contrasts filelistSess{sessi}]);   
    
    disp ([ 'fnames = ' fnames{:} newline 'win_width = ' num2str(win_width) newline 'mf = ' num2str(mf) newline ...
            'meanInTime = ' num2str(meanInTime) newline 'meanInFreq = ' num2str(meanInFreq) newline 'TG = ' num2str(TG) newline ...
            'bline = ' num2str(bline) newline 'acrossTrials = ' num2str(acrossTrials) newline 'batch_bin = ' num2str(batch_bin)  ]);
        
    
    cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline);
% %     % % check the normalization across trials works 
% %     ids = str2num(cell2mat(cfg_contrasts.oneListIds));
% %     d2p = squeeze(mean(cfg_contrasts.oneListPow(ids(:, 1) ==1, :, :, :)));
% %     imagesc(squeeze(d2p(1,:,:))); colorbar % % % because of the normalization across trials

 
    cfg_contrasts.contr2save = contr2save';
    cfg_contrasts.n2s = n2s;
    cfg_contrasts.loadSurr = loadSurr;
    cfg_contrasts.batch_bin = batch_bin;
    

    if strcmp(avMeth,'pow')
        cfg.avRep = 1; 
        cfg_contrasts = average_repetitions(cfg, cfg_contrasts);
    end
    
    [out_contrasts] = create_contrasts_WM (cfg_contrasts);
        
 
    for freqi = 1:length(frequncies2test)
        fname = fnames{freqi};
        mkdir ([paths.results.bands fname]);
        cd ([paths.results.bands fname]);
        f           = frequncies2test{freqi}; 
        % % the rsa_WM function saves the similarity matrices (TG=1) or diagonals (TG = 0)
        rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, sessi, TG, aVTime)
        cd ..
    end
 
end

cd .. 

 
%%process Folders and create one file per condition including all subjects

clearvars -except region
paths = load_paths_WM(region, 'none');
currentDir = pwd; 
mkdir(paths.results.bands)
cd (paths.results.bands)
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
 
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym2use frequncies2test currentDir
    
    direct = fold(foldi);
    cd (direct.name)
 
    processFoldersWM; 

    cd .. 
end

cd (currentDir);
disp ('done')

 
%% LOAD all conditions
clearvars

region = 'pfc'; 
paths = load_paths_WM(region, 'none');

contrasts = {
              %'SCSP_M2M2' 'SCDP_M2M2';
              'SCSP_EM1' 'SCDP_EM1';
              
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


%%PLOT
clc
%clear all; 
tic; 

%define conditions 
cond1 = 'SCSP_EM1';
cond2 = 'SCDP_EM1';

cond1B = eval(cond1);
cond2B = eval(cond2);
 
cfg             =       [];
cfg.subj2exc    =       [18 22];% vvs;
%cfg.subj2exc   =       [1]; % pfc
cfg.clim        =       [-.01 .01];
cfg.climT       =       [-7 7]; %color scale for t-map
cfg.saveimg     =       1;
cfg.res         =       '100_norm'; %'100_perm'; '100_norm'
cfg.cut2        =       '1-4'; %1-1 1-4 4-4
cfg.cond1       =       cond1;
cfg.cond2       =       cond2;
cfg.lwd1        =       2; %baseline 
cfg.lwd2        =       2; %significant outline
cfg.remClust    =       0; 
cfg.plot1clust  =       0;  
cfg.clust2plot  =       [];  %VVS > 3-8: 4; 9-12: 6; 13-29: 8; 30-75: 7; 75-150: [3 17]; vector of pixels to print
cfg.all_cond1   =       cond1B; 
cfg.all_cond2   =       cond2B; 
cfg.alpha       =       0.05; 



[out_real]   = plot_reinst_map_wm(cfg);
 
toc




%% PERMUTATIONS

%%perm
cfg_perm                    =       [];
cfg.runperm                 =       1;  
cfg_perm.n_perm             =       1000; 
cfg_perm.savePerm           =       1;
cfg_perm.out_real           =       out_real;
cfg_perm.pval               =       0.05;
cfg_perm.cond1              =       cond1;
cfg_perm.cond2              =       cond2;
cfg_perm.ids_all_cond       =       [];
 

[out_perm] = myPerm(cfg_perm);
if ~isempty (out_perm.max_clust_sum_real)
    max_clust_sum = out_perm.max_clust_sum;
    obs = max(out_real.all_clust_tsum_real(:,1));
    %allAb = max_clust_sum(abs(max_clust_sum) > obs);
    allAb = max_clust_sum((max_clust_sum) > obs);
    p = (1 - (cfg_perm.n_perm - (length (allAb) ) )  /cfg_perm.n_perm) + (1/cfg_perm.n_perm);
    disp (['p = ' num2str(p)]);
end


%%
n_perm = 1000; 
max_clust_sum = out_perm.max_clust_sum;
obs = max(out_perm.max_clust_sum_real);
allAb = max_clust_sum((max_clust_sum) > obs);
p = (1 - (n_perm - (length (allAb) ) )  /n_perm) + (1/n_perm);
disp (['p = ' num2str(p)]);



%% stats

allAb = max_clust_sum(max_clust_sum > 121.939054751820);
p =1 - (nPerm - (length (allAb)+1) )  /nPerm;

disp (['p = ' num2str(p)]);



%%
figure
histogram(out_perm.max_clust_sum)



%% analysis high temporal resolution: Fig2E

clearvars

curr_fold = pwd; 
cd D:\_WM\analysis\pattern_similarity\pfc\50010ms\EE\avRepet\bands\category\TG\3-150Hz

contrasts = {'DISC_EE' 'DIDC_EE'};
c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    if exist("all_IDs")
        idData{i,:} = all_IDs;
    end
end


%% extract only diagonal for each subject in each conditoon

%subj2exc    =       [18 22];% vvs;
subj2exc   =       [1]; % pfc

DISC_EE(subj2exc) = []; DIDC_EE(subj2exc) = []; 

cond1 = cell2mat(cellfun(@(x) diag(squeeze(mean(x, 1, 'omitnan')))', DISC_EE, 'un', 0));
cond2 = cell2mat(cellfun(@(x) diag(squeeze(mean(x, 1, 'omitnan')))', DIDC_EE, 'un', 0));



%%
mc1 = mean(cond1); 
mc2 = mean(cond2); 

[h p ci ts] = ttest(cond1-cond2); 
hb = h; hb(hb==0) = nan; hb(hb==1) = 0; 

figure()
%plot(mc1, 'r', 'linewidth' ,2); hold on; 
%plot(mc2, 'b', 'linewidth' ,2); 
plot(mc1-mc2, 'linewidth' ,2);hold on; 
plot(hb, 'linewidth' ,2);


exportgraphics(gcf, 'myIm.png', 'Resolution', 200)


%% shaded error version 

mc1 = mean(cond1); 
mc2 = mean(cond2); 

[h p ci ts] = ttest(cond1-cond2); 
hb = h; hb(hb==0) = nan; hb(hb==1) = 0; 

diff2U = cond1-cond2; 
mDiff = mean(diff2U); 
stdDiff = std(diff2U); 
seDiff = std(diff2U)/sqrt(size(diff2U, 1))
times = -.75:.01:.75

figure()
shadedErrorBar(times, mDiff, seDiff, 'b', 1); hold on; 
plot(times, hb, 'linewidth' ,2);


exportgraphics(gcf, 'myIm.png', 'Resolution', 200)




%% analysis high temporal resolution: Fig2E

cd D:\_WM\analysis\pattern_similarity\pfc\50010ms\EE\avRepet\bands\category\TG\3-150Hz
load all % this file contains the traces from VVS and PFC computed in the cells above

mc1PFC = mean(cond1_PFC); 
mc2PFC = mean(cond2_PFC); 
cPFC = cond1_PFC-cond2_PFC;
mcPFC = mean(cPFC)
stdPFC = std(cPFC); 
sePFC = stdPFC/sqrt(size(cPFC, 1))

mc1VVS = mean(cond1_VVS); 
mc2VVS = mean(cond2_VVS); 
cVVS = cond1_VVS - cond2_VVS; 
mcVVS = mean(cVVS); 
stdVVS = std(cVVS); 
seVVS = stdVVS/sqrt(size(cVVS, 1))


l2y = -.015; 
[hPFC p ci ts] = ttest(cPFC); 
hbPFC = hPFC; hbPFC(hbPFC==0) = nan; hbPFC(hbPFC==1) = l2y+.0015; 
[hVVS p ci ts] = ttest(cVVS); 
hbVVS = hVVS; hbVVS(hbVVS==0) = nan; hbVVS(hbVVS==1) = l2y +.0045; 


[hBoth p ci ts] = ttest2(cVVS, cPFC); clustinfo = bwconncomp(hBoth);
clustinfo = bwconncomp(hBoth);
clear allTObs

for pixi = 1:length(clustinfo.PixelIdxList)
     tObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
end

hbBoth = hBoth; hbBoth(hBoth==0) = nan; hbBoth(hbBoth==1) = l2y +.00775; 

times = -.75:.01:.75

figure()
%plot(mc1, 'r', 'linewidth' ,2); hold on; 
%plot(mc2, 'b', 'linewidth' ,2); 
shadedErrorBar(times, mcPFC, sePFC, 'k', 1); hold on; 
shadedErrorBar(times, mcVVS, seVVS, 'r', 1); hold on; 
plot(times, hbPFC, 'k', 'linewidth' ,10);
plot(times, hbVVS, 'r', 'linewidth' ,10);
plot(times, hbBoth, 'g', 'linewidth' ,10);
set(gca, 'FontSize', 20, 'xlim', [-.5, .75], 'ylim', [l2y 0.06])
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 2);


exportgraphics(gcf, 'myIm.png', 'Resolution', 200)

%% analysis high temporal resolution: Fig2E randomly subselecting same number of subjects 1000 times

cd D:\_WM\analysis\pattern_similarity\pfc\50010ms\EE\avRepet\bands\category\TG\3-150Hz
load all % this file contains the traces from VVS and PFC computed in the cells above

mc1PFC = mean(cond1_PFC); 
mc2PFC = mean(cond2_PFC); 
cPFC = cond1_PFC-cond2_PFC;
mcPFC = mean(cPFC)
stdPFC = std(cPFC); 
sePFC = stdPFC/sqrt(size(cPFC, 1));

ids2u = randsample(26, 16);

mc1VVS = mean(cond1_VVS(ids2u, :)); 
mc2VVS = mean(cond2_VVS(ids2u, :)); 
cVVS = cond1_VVS(ids2u, :) - cond2_VVS(ids2u, :); 
mcVVS = mean(cVVS); 
stdVVS = std(cVVS); 
seVVS = stdVVS/sqrt(size(cVVS, 1))


l2y = -.015; 
[hPFC p ci ts] = ttest(cPFC); 
hbPFC = hPFC; hbPFC(hbPFC==0) = nan; hbPFC(hbPFC==1) = l2y+.0015; 
[hVVS p ci ts] = ttest(cVVS); 
hbVVS = hVVS; hbVVS(hbVVS==0) = nan; hbVVS(hbVVS==1) = l2y +.0045; 


[hBoth p ci ts] = ttest2(cVVS, cPFC); 
hbBoth = hBoth; hbBoth(hBoth==0) = nan; hbBoth(hbBoth==1) = l2y +.00775; 

times = -.75:.01:.75

figure()
%plot(mc1, 'r', 'linewidth' ,2); hold on; 
%plot(mc2, 'b', 'linewidth' ,2); 
shadedErrorBar(times, mcPFC, sePFC, 'k', 1); hold on; 
shadedErrorBar(times, mcVVS, seVVS, 'r', 1); hold on; 
plot(times, hbPFC, 'k', 'linewidth' ,10);
plot(times, hbVVS, 'r', 'linewidth' ,10);
plot(times, hbBoth, 'g', 'linewidth' ,10);
set(gca, 'FontSize', 20, 'xlim', [-.5, .75], 'ylim', [l2y 0.06])
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 2);


exportgraphics(gcf, 'myIm.png', 'Resolution', 200)






%% analysis high temporal resolution: Fig2E randomly subselecting same number of subjects 1000 times

cd D:\_WM\analysis\pattern_similarity\pfc\50010ms\EE\avRepet\bands\category\TG\3-150Hz
load all % this file contains the traces from VVS and PFC computed in the cells above

nPerm = 1000; 

clear max_clust_sum_perm
for permi = 1:nPerm
    cPFC = cond1_PFC-cond2_PFC;
    ids2u = randsample(26, 16);
    cVVS = cond1_VVS(ids2u, :) - cond2_VVS(ids2u, :); 
        
    [hBoth p ci ts] = ttest2(cVVS, cPFC); 
    t = ts.tstat; 

    clustinfo = bwconncomp(hBoth);
    clear allTObs
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
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

%% 

histogram (max_clust_sum_perm); hold on; 
plot([tObs tObs], get(gca,'ylim'),'r','lineWidth', 3);




%% stats: condition label shuffling 

clearvars -except cond1_PFC cond1_VVS cond2_PFC cond2_VVS

nPerm = 1000; 


for permi =1:1000



    mc1PFC = mean(cond1_PFC); 
    mc2PFC = mean(cond2_PFC); 
    cPFC = cond1_PFC-cond2_PFC;
    
    mc1VVS = mean(cond1_VVS); 
    mc2VVS = mean(cond2_VVS); 
    cVVS = cond1_VVS - cond2_VVS; 

    junts = [ cVVS; cPFC];
    rndIndex = [zeros(1, size(cVVS, 1)), ones(1 ,size(cPFC, 1))]';
    %rndIndex = rndIndex(randperm(length(rndIndex)));
    junts1 = junts(rndIndex==0, :);
    junts2 = junts(rndIndex==1, :);
        
    [hPerm p ci ts] = ttest2(junts1, junts2); 
    tPerm = ts.tstat; 

    clustinfo = bwconncomp(hPerm);
    [numPixPermi(permi) maxi] = max([0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    if numPixPermi(permi) > 0
        max_clust_sum(permi,:) = sum (tPerm(clustinfo.PixelIdxList{maxi-1}));
    else
        %disp (['no significant cluster in permutation ' num2str(permi)]);
        max_clust_sum(permi,:) = 0; 
    end

    
end

%% plot histogram

figure()
histogram(max_clust_sum); hold on; 
scatter(121.939054751820, 0, 2000, '.', 'r' )



%% stats

allAb = max_clust_sum(max_clust_sum > 121.939054751820);
p =1 - (nPerm - (length (allAb)+1) )  /nPerm;

disp (['p = ' num2str(p)]);





















%%


 
 













%%
