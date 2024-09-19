%% load cfg_contrasts
%% 

clearvars 
region              = 'vvs';
paths = load_paths_WM(region, 'none');
filelistSess = getFilesWM(paths.out_contrasts);


frequncies2test = [{3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = { '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

% frequncies2test = [{13:29}]';
% fnames = { '13-29Hz'}'; fnames = fnames';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 0; 
meanInFreq          = 0; 
avMeth              = 'none';  % average across image repetitions or not
TG                  = 1; %temporal generalization
%contr2save          = { 'DISC_EM2' 'DIDC_EM2'}; %
%contr2save          = { 'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2'}; %
%contr2save          = { 'DISC_EM11' 'DIDC_EM11' 'DISC_EM12' 'DIDC_EM12' 'DISC_EM13' 'DIDC_EM13'}; %
%contr2save          = { 'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3' 'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE' 'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2'}; %
%contr2save          = {'DISC_EE' 'DIDC_EE' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'SCSP_M2M2' 'SCDP_M2M2'}; %{};
%contr2save          = {'DISC_EM1' 'DIDC_EM1'}; %{};
%contr2save          = {'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3'}; %{};
%contr2save          = {'DISC_EM21' 'DIDC_EM21' 'DISC_EM22' 'DIDC_EM22' 'DISC_EM23' 'DIDC_EM23'}; %{};
%contr2save          = {'SCSP_EM2' 'SCDP_EM2' 'DCSP_EM2' 'DCDP_EM2'}; %{};
contr2save          = {'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2' 'SCSP_EM2' 'SCDP_EM2' 'DCSP_EM2' 'DCDP_EM2'}; %{};
%contr2save          = {'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE'}; %{};

bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
n2s                 = 5000;
loadSurr            = 0; 
zScType             = 'sess'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
aVTime              = 0; % Average or not in time the feature vectors
 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3

    clearvars -except region paths filelistSess frequncies2test fnames win_width mf meanInTime meanInFreq avMeth TG ...
        contr2save bline acrossTrials batch_bin n2s loadSurr zScType aVTime sessi
 

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
    

    %cfg_contrasts = take_only_lateral_electrodes(cfg_contrasts);


    if size(cfg_contrasts.chanNames, 1) > 1


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
 

clearvars
region              = 'pfc';
paths = load_paths_WM(region, 'none');
filelistSess = getFilesWM(paths.out_contrasts);


frequncies2test = [{3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = { '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

%frequncies2test = [{13:29}]';
%fnames = { '13-29Hz'}'; fnames = fnames';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 0; 
meanInFreq          = 0; 
avMeth              = 'none';  % average across image repetitions or not
TG                  = 1; %temporal generalization
%contr2save          = { 'DISC_EM2' 'DIDC_EM2'}; %
%contr2save          = { 'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2'}; %
%contr2save          = { 'DISC_EM11' 'DIDC_EM11' 'DISC_EM12' 'DIDC_EM12' 'DISC_EM13' 'DIDC_EM13'}; %
%contr2save          = { 'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3' 'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE' 'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2'}; %
%contr2save          = {'DISC_EE' 'DIDC_EE' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'SCSP_M2M2' 'SCDP_M2M2'}; %{};
%contr2save          = {'DISC_EM1' 'DIDC_EM1'}; %{};
%contr2save          = {'DISC_EE1' 'DIDC_EE1' 'DISC_EE2' 'DIDC_EE2' 'DISC_EE3' 'DIDC_EE3'}; %{};
%contr2save          = {'DISC_EM21' 'DIDC_EM21' 'DISC_EM22' 'DIDC_EM22' 'DISC_EM23' 'DIDC_EM23'}; %{};
%contr2save          = {'SCSP_EM2' 'SCDP_EM2' 'DCSP_EM2' 'DCDP_EM2'}; %{};
contr2save          = {'SCSP_M2M2' 'SCDP_M2M2' 'DCSP_M2M2' 'DCDP_M2M2' 'SCSP_EM2' 'SCDP_EM2' 'DCSP_EM2' 'DCDP_EM2'}; %{};
%contr2save          = {'SCSP_EE' 'SCDP_EE' 'DCSP_EE' 'DCDP_EE'}; %{};

bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
n2s                 = 5000;
loadSurr            = 0; 
zScType             = 'sess'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
aVTime              = 0; % Average or not in time the feature vectors
 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3

    clearvars -except region paths filelistSess frequncies2test fnames win_width mf meanInTime meanInFreq avMeth TG ...
        contr2save bline acrossTrials batch_bin n2s loadSurr zScType aVTime sessi
 

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
    

    %cfg_contrasts = take_only_lateral_electrodes(cfg_contrasts);


    if size(cfg_contrasts.chanNames, 1) > 1


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

clear, clc 

region = 'pfc'; 
paths = load_paths_WM(region, 'none');

tic
contrasts = {
              %'DISC_M1A' 'DIDC_M1A';

              %'DISC_EE1' 'DIDC_EE1';
              %'DISC_EE2' 'DIDC_EE2';
              %'DISC_EE3' 'DIDC_EE3';
              
               'SCSP_M2M2' 'SCDP_M2M2' ; 
              % 'DCSP_M2M2' 'DCDP_M2M2';

              %'SCSP_EM2' 'SCDP_EM2' ; 
              %'DCSP_EM2' 'DCDP_EM2';

              %'SCSP_EE' 'SCDP_EE' ; 
              %'DCSP_EE' 'DCDP_EE';

              %'DISC_EM11' 'DIDC_EM11';
              %'DISC_EM12' 'DIDC_EM12';
              %'DISC_EM11' 'DIDC_EM11';

               % 'DISC_EM21' 'DIDC_EM21';
               % 'DISC_EM22' 'DIDC_EM22';
               % 'DISC_EM23' 'DIDC_EM23';

              %'DISC_EM2' 'DIDC_EM2';

              %'DISC_EE' 'DIDC_EE';
              %'DISC_EM2' 'DIDC_EM2';
              
              
              %'DISC_M2M2' 'DIDC_M2M2';

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

%% PLOT
clc

%define conditions 
cond1 = 'SCSP_M2M2'; 
cond2 = 'SCDP_M2M2'; 


cond1B = eval(cond1);
cond2B = eval(cond2);
 
cfg             =       [];
%cfg.subj2exc   =       []; % pfc LATERAL
%cfg.subj2exc    =       [18 22];% vvs;
%cfg.subj2exc    =       [18];% vvs;
%cfg.subj2exc   =       [1]; % pfc
cfg.subj2exc   =       []; % pfc LATERAL
cfg.clim        =       [-.01 .01];
cfg.climT       =       [-7 7]; %color scale for t-map
cfg.saveimg     =       1;
cfg.res         =       '100_perm'; %'100_perm'; '100_norm'
cfg.cut2        =       '4-4'; %1-1 1-4 4-4
cfg.lwd1        =       2; %baseline 
cfg.lwd2        =       2; %significant outline
cfg.remClust    =       0; 
cfg.plot1clust  =       1;  
cfg.clust2plot  =       [8];  %VVS EMS maintenance > 3-8: 4; 9-12: 6; 13-29: 8; 30-75: 7; 75-150: [3 17]; vector of pixels to print
cfg.plotTrend   =       0; 
cfg.cond1       =       cond1;
cfg.cond2       =       cond2;
cfg.all_cond1   =       cond1B; 
cfg.all_cond2   =       cond2B; 
cfg.alpha       =       0.05; 


[out_real]   = plot_reinst_map_wm(cfg);

%tObs = squeeze(mean(mean(out_real.meanReal_cond1(:, 1:8, 6:15) - out_real.meanReal_cond2(:, 1:8, 6:15), 2), 3));
%[h p ci ts] = ttest(tObs); 
 
close all 





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




%% LOAD all conditions
clearvars

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

%% ANOVA EM21 EM22 EM23
clc
%clear all; 

%sub2exc = [1]; 
%sub2exc = [18 22]; 
sub2exc = [18]; 

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

tP1 = 6:13; 
tP2 = 6:25; 
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

checkTimes(1, :) = times; 
checkTimes(2, :) = hbBoth; 

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

nPerm = 10000; 

clear max_clust_sum_perm allPs
for permi = 1:nPerm
    cPFC = cond1_PFC-cond2_PFC;
    ids2u = randsample(26, 15);
    cVVS = cond1_VVS(ids2u, :) - cond2_VVS(ids2u, :); 

    mcPFC = mean(cPFC(:, 91:143), 2); 
    mcVVS = mean(cVVS(:, 91:143), 2); 


    [h p ci ts] = ttest2(mcVVS, mcPFC); 
    allPs(permi, :) = p; 

        
    % [hBoth p ci ts] = ttest2(cVVS, cPFC); 
    % t = ts.tstat; 
    % 
    % clustinfo = bwconncomp(hBoth);
    % clear allTObs
    % if ~isempty(clustinfo.PixelIdxList)
    %     for pixi = 1:length(clustinfo.PixelIdxList)
    %          allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
    %     end
    % else
    %     allTObs(permi,:) = 0;
    % end
    % 
    % if exist('allTObs')
    %     [max2u id] = max(abs(allTObs));
    %     max_clust_sum_perm(permi,:) = allTObs(id); 
    % else
    %     max_clust_sum_perm(permi,:) = 0; 
    % end
    % 

end


%% count how many times there is a significant effect when VVS N = 15
numberOfTimes = sum(allPs < 0.05);
perc2report = (numberOfTimes/nPerm)*100

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



%% compare pfc M2M2 vs EM2
clearvars

cd D:\_WM\analysis\pattern_similarity\pfc\100ms\EM2-M2_together
region = 'pfc'; 
paths = load_paths_WM(region, 'none');

contrasts = {
              'DISC_M2M2' 'DIDC_M2M2';
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


%% take all means


c1 = cell2mat(cellfun(@mean, DISC_EM2, 'un', 0)); 
c2 = cell2mat(cellfun(@mean, DIDC_EM2, 'un', 0)); 
c3 = cell2mat(cellfun(@mean, DISC_M2M2, 'un', 0)); 
c4 = cell2mat(cellfun(@mean, DIDC_M2M2, 'un', 0)); 

c5 = c1-c2; 
c6 = c3-c4; 

c7 = c1-c3; 

%transfMetric = mean(mean(c7(:, 8:35, 8:35), 2, 'omitnan'), 3, 'omitnan'); 
transfMetric = mean(mean(c7(:, 8:15, 8:15), 2, 'omitnan'), 3, 'omitnan'); 


[h p ci ts ] = ttest(transfMetric)









%%


 
 













%%
