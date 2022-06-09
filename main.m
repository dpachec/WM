%% load cfg_contrasts
%% 

clearvars
region              = 'vvs';
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts_path);


frequncies2test = [{3:54} {3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = {'3-54Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-38Hz' '39-54Hz' }'; fnames = fnames';

%frequncies2test = [{13:29}]';
%fnames = {'13-29Hz'}'; fnames = fnames';

%frequncies2test = [{3:54}]';
%fnames = {'3-54Hz'}'; fnames = fnames';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 1; 
meanInFreq          = 0; 
takeElec            = 0; 
takeFreq            = 0;
TG                  = 1; %temporal generalization
contr2save          = {'SISC_EM2U' 'DISC_EM2U' 'DIDC_EM2U'}; %{};
%contr2save          = {'SISC_EE' 'DISC_EE' 'DIDC_EE' 'SISC_EM2' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'DISC_M2123V1' 'DIDC_M2123V1' 'DISC_M2123V2' 'DIDC_M2123V2' 'DISC_M2123CNCV1' ...
%                          'DIDC_M2123CNCV1' 'DISC_M2123CNCV2' 'DIDC_M2123CNCV2' 'DISC_M2123NC' 'DIDC_M2123NC' ...
%                          'DISC_M2A' 'DIDC_M2A' 'DISC_M2A123' 'DIDC_M2A123'}; 
bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
n2s                 = 1000000;
loadSurr            = 0; 
zScType             = 'allTrials'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
avMeth              = 'pow';  
 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.out_contrasts_path filelistSess{sessi}]);   
    
    disp ([ 'fnames = ' fnames{:} newline ...
            'win_width = ' num2str(win_width) newline ...
            'mf = ' num2str(mf) newline ...
            'meanInTime = ' num2str(meanInTime) newline ...
            'meanInFreq = ' num2str(meanInFreq) newline ...
            'TG = ' num2str(TG) newline ...
            'bline = ' num2str(bline) newline ...
            'acrossTrials = ' num2str(acrossTrials) newline ...
            'batch_bin = ' num2str(batch_bin);
            ]);
        
    
    cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline);
 
    cfg_contrasts.contr2save = contr2save;
    cfg_contrasts.n2s = n2s;
    cfg_contrasts.loadSurr = loadSurr;
    cfg_contrasts.batch_bin = batch_bin;
    

    if strcmp(avMeth,'pow')
        cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    [out_contrasts] = create_contrasts_WM (cfg_contrasts);
        
        
    if ~exist('idxCH')
       idxCH = []; 
    end
    if ~exist('idxF')
       idxF = []; 
    end
  
    for freqi = 1:length(frequncies2test)
        %create folder
        fname = fnames{freqi};
        mkdir ([paths.results_path fname]);
        cd ([paths.results_path fname]);
        f           = frequncies2test{freqi}; 
 
        rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, takeElec, takeFreq, idxCH, idxF, sessi, TG, 0, paths)
 
        cd ..
    end
 
end

cd .. 

%diary off

%disp('done');
%disp(string(datetime));
 
%% 
clear 
paths = load_paths_WM; 


fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
 
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym2use frequncies2test
    
    direct = fold(foldi);
    cd (direct.name)
 
    processFoldersWM; 

    cd .. 
end

disp ('done')

%% run permutations and plot

disp ('Starting...');
disp(string(datetime));

paths = load_paths_WM; 
cd (paths.results_path)
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
% contrasts = {   'DISC_M2123V1' 'DIDC_M2123V1'; ... 
%                 'DISC_M2123V2' 'DIDC_M2123V2' ; ...
%                 'DISC_M2123CNCV1' 'DIDC_M2123CNCV1' ; ... 
%                 'DISC_M2123CNCV2' 'DIDC_M2123CNCV2'; ...
%                 'DISC_M2123NC' 'DIDC_M2123NC';...
%                 'DISC_M2A' 'DIDC_M2A';...
%                 'DISC_M2A123' 'DIDC_M2A123';
%              };
                          
                    
% contrasts = {'SISC_EE' 'DISC_EE'; ... 
%              'DISC_EE' 'DIDC_EE'; ...
%              'SISC_EM2' 'DISC_EM2'; ... 
%              'DISC_EM2' 'DIDC_EM2'; ...
%              'DISC_M2M2' 'DIDC_M2M2';...
%              };
%  

contrasts = {'DISC_EM2' 'DIDC_EM2'};

cmaps2use = {[-.15 .15]};
cmaps2use = repmat(cmaps2use, 1, 50);
 
 
%perms2use = {'1-1' '1-1' '1-4' '1-4' '4-4'};
%perms2use = { '-4-0' '-4-0'};
%perms2use = {'4-4' '4-4' '4-4' '4-4' '4-4' '4-4' '4-4'};
perms2use = {'1-4'};
t2use = {'100_norm' '100_perm'};

 
cmapi = 1;
 
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym2use
    
    direct = fold(foldi);
    cd (direct.name)
 
    processFoldersWM;
    
    c = unique (contrasts);
    d = cellfun(@(x) [x '_id'], c, 'un', 0);
    for i = 1:length(c) 
        load([c{i} '.mat']);
        contrData{i,:} = eval(c{i});
        idData{i,:} = all_IDs;
    end

    region = 'pfc'; 
    noAv = 0;
    [out_c out_id] = averageSub_WM (c, d, contrData, idData, region, noAv);
    for i = 1:length(out_c) 
        eval([c{i} ' = out_c{i};']);
        eval([d{i} ' = out_id{i};']);
    end

    
    
 
 
    for permNoperm = 1:1 %1 = just plots (no perm) (2:2) just permutation
        for imi = 1:size(contrasts, 1)
 
            tic; clear all_cond1 all_cond2 all_cond1_A all_cond2_A;
 
            cond1 = contrasts{imi, 1};
            cond2 = contrasts{imi, 2};
 
 
            all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);
 
            %global parameters
            if strcmp(region, 'vvs')
                subj2exc        =       [18 22]; % pfc
            elseif strcmp(region, 'pfc')
                subj2exc        =       [1]; % pfc
            else
                subj2exc        =       []; % pfc
            end
            
            if permNoperm == 1
                runperm         =       0;  
            else
                runperm         =       1; 
            end 
            n_perm          =       1000;
            saveperm        =       1;
            cfg             =       [];
            cfg.clim        =       cmaps2use{cmapi}; cmapi = cmapi+1;  
            cfg.climT       =       [-3 3]; %color scale for t-map
            cfg.square      =       1;
            cfg.saveimg     =       1;
            cfg.enc_ret     =      'e';
            cfg.lim         =       'final'; %'no'  -   %'edge' - % 'final' -- 'jackk'
            cfg.res         =       t2use{permNoperm}; %'100_no_bas'
            cfg.cut2        =       perms2use{imi}; %4 3 2.5 2 
            cfg.cond1       =       cond1;
            cfg.cond2       =       cond2;
            cfg.runperm     =       runperm;
            test2use        =       'ttest'; %'wilcox'; %'ttest';
            cfg.remClust    =       1; %remove clusters from plot (only if run permutation is true)
            cfg.plot1clust  =       0; %to plot just one cluster selected in the following line
            cfg.clust2plot  =       11;
            cfg.subj2exc    =       subj2exc;
            cfg.lwd1        =       2; %baseline 
            cfg.lwd2        =       2; %significant outline
 
 
            %exclude subjects
            if subj2exc > 0
                all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
                all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
            end
            

            % % % % duplicate for EM2
            typec = strsplit(cond1, '_');
            if length(typec) > 1
                if strcmp(typec{2}, 'EM2')
                    dupSym = 0; % not duplicate in plot
                    for subji = 1:length(all_cond1_A)
                        d2c = all_cond1_A{subji};
                        for triali = 1:size(d2c, 1)
                            d2ct = squeeze(d2c(triali, :,:)); 
                            d2ct = triu(d2ct.',1) + tril(d2ct);
                            d2c(triali, :, :) = d2ct; 
                        end
                        all_cond1_A{subji} = d2c;

                        d2c = all_cond2_A{subji};
                        for triali = 1:size(d2c, 1)
                            d2ct = squeeze(d2c(triali, :,:)); 
                            d2ct = triu(d2ct.',1) + tril(d2ct);
                            d2c(triali, :, :) = d2ct; 
                        end
                        all_cond2_A{subji} = d2c;
                    end
                else
                    dupSym          =       1; % duplicate in plot
                end
            end

            
            
            cfg.all_cond1_A =       all_cond1_A;
            [cfg_plot]   =      set_reinst_plot_wm (cfg);
            % apply limits
            clear all_cond1 all_cond2;
            n_subj = length(all_cond1_A);
            for si = 1:n_subj   
                all_cond1{si} = all_cond1_A{si}(:,cfg_plot.mlimE,cfg_plot.mlimR);
                all_cond2{si} = all_cond2_A{si}(:,cfg_plot.mlimE,cfg_plot.mlimR);
            end
 
            %calculate real differences
            cfg_real.pT             =   0;          %pick samples from the condition with more trials
            cfg_real.pval           =   0.05; 
            cfg_real.test2use       =   test2use;
            cfg_real.all_cond1      =   all_cond1;
            cfg_real.all_cond2      =   all_cond2;
            cfg_real.connectivity   =   8;
 
 
            [out_real]         =   real_diff_reinst_wm(cfg_real);
 
 
 
            %%perm
            cfg_perm                    =       [];
            cfg_perm.savePerm           =       saveperm;
            cfg_perm.out_real           =       out_real;
            cfg_perm.all_cond1          =       all_cond1; 
            cfg_perm.all_cond2          =       all_cond2;
            cfg_perm.n_perm             =       n_perm; 
            cfg_perm.pT                 =       0; % select random samples from the condition with more trials
            cfg_perm.pval               =       0.05;
            cfg_perm.bins               =       size(all_cond1_A, 2);
            cfg_perm.cond1              =       cond1;
            cfg_perm.cond2              =       cond2;
 
            cfg_perm.ids_all_cond       =       [];
 
 
            if runperm 
                [out_perm] = myPerm(cfg_perm);
                cfg_plot.sigMH_thres = out_perm.sigMH_thres;
            else
                cfg_plot.sigMH_thres = out_real.sigMH_real;
            end
 
            cfg_plot.out_real       = out_real;
            cfg_plot.dupSym         = dupSym; 
            cfg_plot.lwd1        =  cfg.lwd1; %baseline 
            cfg_plot.lwd2        =  cfg.lwd2; %significant outline
            cfg_plot.plotD       = 0; 
 
            if strcmp(test2use, 'wilcox')
                cfg_plot.diff     = 'diff'; % or 'tmap';
            else
                cfg_plot.diff     = 'tmap'; % or 'tmap';
            end
 
            [myPlotR] = plot_reinst_map_wm(cfg_plot);
 
            if runperm & ~isempty (out_perm.max_clust_sum_real)
                max_clust_sum = out_perm.max_clust_sum;
                obs = max(out_real.all_clust_tsum_real(:,1));

                allAb = max_clust_sum(abs(max_clust_sum) > abs(obs));
                p =1 - (n_perm - (length (allAb)+1) )  /n_perm;
                disp (['p = ' num2str(p)]);
            end
 
            toc
 
 
 
        end
    end 
cd ..
end
 
cd ..
disp ('all plots done');
disp(string(datetime));
diary off
 
 

%% PLOT 2 CLEAN
%EEG = ALLEEG{2};
%eventChannel = 'Trigger';%'POL DC13'; % 'POL DC13'; %
%TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
chanids = [1:52];
eegplot(EEG.data(chanids, :), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 300, 'events', EEG.event, 'command', 'TMPREJ');%, 'command', 'TMPREJ'
 
 
 
%% Plot all trials
%%load cfg_contrasts
clear
tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load(sublist{subji});
    s_all{subji} = cfg_contrasts; 
    %size(cfg_contrasts.oneListPow)
end
s_all = s_all';
toc
 
%% plot by trials
bline               = [3 7];
acrossTrials       = 1;
zScType = 'allTrials'; 
 
for subji = 1:length(sublist)
   cfg_contrasts = s_all{subji};    
   cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
   oneListPow = cfg_contrasts.oneListPow;
   for triali = 1:size(oneListPow, 1)
       chanN = size(oneListPow, 2);
       pag2plot = ceil(chanN/100);
       cstri = 0;
       for pagi = 1:pag2plot
           figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
           csubi = 1;
           startChan = (cstri*100) +1
           endChan = (cstri*100) +100
           if endChan > chanN
               endChan = chanN
           end
           
           for chani = startChan:endChan
                subplot (10, 10, csubi)
                imagesc(squeeze(oneListPow(triali,chani,:,:))); 
                title(num2str(triali))
                set(gca,'YDir','normal')
                csubi = csubi + 1; 
           end
           
           cstri = cstri + 1;
           filename = [num2str(subji) '_tr_' num2str(triali) '_' num2str(pagi)];
           export_fig(pagi, filename,'-r200');    
           close all;  
       end
    end
   
end
 
 
 
 

%% plot by channels
bline               = [3 7];
acrossTrials       = 1;
zScType = 'allTrials'; 
 
for subji = 18:18 %1:length(sublist)
   cfg_contrasts = s_all{subji};    
   cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
   oneListPow = cfg_contrasts.oneListPow;
   for chani = 1:size(oneListPow, 2)
       trialN = size(oneListPow, 1);
       pag2plot = ceil(trialN/100);
       cstri = 0;
       for pagi = 1:pag2plot
           figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
           csubi = 1;
           startTri = (cstri*100) +1
           endTri = (cstri*100) +100
           if endTri > trialN
               endTri = trialN
           end
           
           for triali = startTri:endTri
                subplot (10, 10, csubi)
                imagesc(squeeze(oneListPow(triali,chani,:,:))); 
                title(num2str(triali))
                set(gca,'YDir','normal')
                csubi = csubi + 1; 
           end
           
           cstri = cstri + 1;
           filename = [num2str(subji) '_' num2str(chani) '_' num2str(pagi)];
           export_fig(pagi, filename,'-r200');    
           close all;  
       end
    end
   
end
 
 
%% plot by channels ONE
bline               = [3 7];
takeAllTrials       = 1;
 
cfg_contrasts = normalize_WM(cfg_contrasts, takeAllTrials, 'blo', bline);
oneListPow = cfg_contrasts.oneListPow;
for chani = 1:size(oneListPow, 2)
   trialN = size(oneListPow, 1);
   pag2plot = ceil(trialN/100);
   cstri = 0;
   for pagi = 1:pag2plot
       figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
       csubi = 1;
       startTri = (cstri*100) +1
       endTri = (cstri*100) +100
       if endTri > trialN
           endTri = trialN
       end

       for triali = startTri:endTri
            subplot (10, 10, csubi)
            imagesc(squeeze(oneListPow(triali,chani,:,:))); 
            title(num2str(triali))
            set(gca,'YDir','normal')
            csubi = csubi + 1; 
       end

       cstri = cstri + 1;
       filename = [num2str(subji) '_' num2str(chani) '_' num2str(pagi)];
       export_fig(pagi, filename,'-r200');    
       close all;  
   end
    
 
 
end



%%
 
imagesc(squeeze(oneListPow(10,20,:,:))); colorbar;
 
 
 
 
 
 
 
 
%% count electrodes
%%load cfg_contrasts
%clearvars
tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load(sublist{subji});
    %s_all{subji} = size(cfg_contrasts.oneListPow, 2); 
    s_all{subji} = cfg_contrasts.chanNames; 
    %size(cfg_contrasts.oneListPow)
end
s_all = s_all';
toc

%% count hippocapmal electrodes

clear s_ids
for subji = 1:length(sublist)

    chanNames = allChans{subji};
    c = string(chanNames(:,5));
    ids = strfind(c, '38');
    id = ~cellfun('isempty',ids);
    if ~isempty(find(id))
        s_ids(subji,:) = 1; 
    else
        s_ids(subji,:) = 0;
    end
    
end

sum(s_ids)
 
 
%%


 
 
 
 
%%
