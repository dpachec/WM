%% LOAD all conditions
%%
clearvars

region = 'vvs'; 
paths = load_paths_WM(region);

contrasts = {
              'DISC_M2M2' 'DIDC_M2M2';
              %'DISC_EE' 'DIDC_EE';
              %'DISC_M2A' 'DIDC_M2A';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
    idData{i,:} = all_IDs;
    %idData{i,:} = [];
end

noAv = 1;
[out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
for i = 1:length(out_c) 
    eval([c{i} ' = out_c{i};']);
    %eval([d{i} ' = out_id{i};']);
end



%% plot diagonal only 

cond1 = 'SISC_EE';
cond2 = 'DISC_EE';

subj2exc = 1; 

all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);
%exclude subjects
if subj2exc > 0
    all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end

cond1 = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond1_A, 'un', 0));
cond2 = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond2_A, 'un', 0));

%% 
mc1 = squeeze(mean(cond1)); 
mc2 = squeeze(mean(cond2)); 
figure()
plot(mc1, 'linewidth' ,2); hold on; 
plot(mc2, 'linewidth' ,2); 
plot(mc1-mc2, 'linewidth' ,2);


%% COND1 - COND2 (NEW LAYOUT 1x3) SHUFFLING AT THE TRIAL LEVEL
%clear all; 
tic; clear all_cond1 all_cond2 all_cond1_A all_cond2_A;

%define conditions 
cond1 = 'DISC_M2M2';
cond2 = 'DIDC_M2M2';
 


all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);
 
%global parameters
subj2exc        =       [1];% vvs;%[1] pfc
%subj2exc        =       [];% vvs;%[1] pfc
runperm         =       1;
plotClust       =       1; 
dupSym          =       1; 
n_perm          =       1000;
saveperm        =       1; 
cfg             =       [];
cfg.clim        =       [-.01 .0115];
cfg.climT       =       [-7 7]; %color scale for t-map
cfg.square      =       1;
cfg.saveimg     =       1;
cfg.enc_ret     =       'e';
cfg.lim         =       'final'; %'no'  -   %'edge' - % 'final' -- 'jackk'
cfg.res         =       '100_perm'; %'100_perm'; '100_norm'
cfg.cut2        =       '4-4'; %4 3 2.5 2 
cfg.cond1       =       cond1;
cfg.cond2       =       cond2;
cfg.runperm     =       runperm;
test2use        =       'ttest'; %'wilcox'; %'ttest';
cfg.remClust    =       1; %remove clusters from plot (only if run permutation is true)
cfg.plot1clust  =       0; %to plot just one cluster selected in the following line
cfg.clust2plot  =       1;
cfg.subj2exc    =       subj2exc;
cfg.lwd1        =       2; %baseline 
cfg.lwd2        =       2; %significant outline
 
 
%exclude subjects
if subj2exc > 0
    all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
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

cfg_plot.lwd1           = cfg.lwd1; 
cfg_plot.lwd2           = cfg.lwd2; 
cfg_plot.out_real       = out_real;
cfg_plot.dupSym         = dupSym; 
cfg_plot.plotD          = 0; 
cfg_plot.plotClust      = plotClust;

 
if strcmp(test2use, 'wilcox')
    cfg_plot.diff     = 'diff'; % or 'tmap';
else
    cfg_plot.diff     = 'tmap'; % or 'tmap';
end

[myPlotR] = plot_reinst_map_wm(cfg_plot);
 

    if runperm & ~isempty (out_perm.max_clust_sum_real)
        max_clust_sum = out_perm.max_clust_sum;

        %first find 
        obs = max(out_real.all_clust_tsum_real(:,1));
        %allAb = max_clust_sum(abs(max_clust_sum) > obs);
        allAb = max_clust_sum((max_clust_sum) > obs);
        p =1 - (n_perm - (length (allAb)-1) )  /n_perm;
        disp (['p = ' num2str(p)]);
        %x = find(max_clust_sum_ranked(:,2) == index)
    end

 
toc
 
%% average in period 
clear m1 m2
for subji = 1:length(all_cond1)
    %m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,6:15,6:45)));
    %m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,6:15,6:45)));
    m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,6:40,6:40)));
    m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,6:40,6:40)));

end

mm1 = squeeze(mean(m1,'omitnan'));
m1A = squeeze(mean(mean(m1,2,'omitnan'),3,'omitnan'));
%mm1 = triu(mm1.',1) + tril(mm1);
mmm1 = squeeze(mean(mm1,'omitnan'));
mm2 = squeeze(mean(m2,'omitnan'));
m2A = squeeze(mean(mean(m2,2,'omitnan'),3,'omitnan'));
%mm2 = triu(mm2.',1) + tril(mm2);
mmm2 = squeeze(mean(mm2,'omitnan'));

%%full maintenance period 
msD = m1A-m2A;
[h p ci ts] = ttest(msD);
t = ts.tstat
p



%% plot one bar
%data.data = [data_LC.data(:,1) data_LC.data(:,2) ];
%data.data = [diff_ERS_all diff_ERS_all_minus]; 
data.data = [msD]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.0025 .0075] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

%[h p ci t] = ttest (data.data(:,1), data.data(:,2));
[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r80');
close all;   

%% encoding - maintenance similarity (average full maintenance period)
ms1 = squeeze(mean(m1, 3, 'omitnan'));
ms2 = squeeze(mean(m2, 3, 'omitnan'));
msD = ms1-ms2; 

figure()
plot(msD'); hold on; 
plot(mean(msD), 'k', 'LineWidth',2)

[h p ci ts] = ttest(msD)
t = ts.tstat
%plot(h, 'LineWidth', 2)

%% 
clear m1 m2
for subji = 1:length(all_cond1)
    m1(subji, :,:) = squeeze(mean(all_cond1{subji}));
    m2(subji, :,:) = squeeze(mean(all_cond2{subji}));
end
ms1 = squeeze(mean(m1, 3, 'omitnan'));
ms2 = squeeze(mean(m2, 3, 'omitnan'));
msD = ms1-ms2; 

mD = mean(msD);
stdD = std(msD);
seD = stdD / sqrt(26);
figure()
shadedErrorBar(1:15, mD,seD, 'r', 1); hold on 
h(h==0) = nan; h(h==1) = 0;
plot(h, 'LineWidth', 2)

%% encoding - maintenance similarity (average full encoding period)
clearvars -except all_cond1 all_cond2 

for subji = 1:length(all_cond1)
    m1(subji, :,:) = squeeze(mean(all_cond1{subji}));
    m2(subji, :,:) = squeeze(mean(all_cond2{subji}));


end
mm1 = squeeze(mean(m1,'omitnan'));
m1A = squeeze(mean(mean(m1,2,'omitnan'),3,'omitnan'));
%mm1 = triu(mm1.',1) + tril(mm1);
mmm1 = squeeze(mean(mm1,'omitnan'));
mm2 = squeeze(mean(m2,'omitnan'));
m2A = squeeze(mean(mean(m2,2,'omitnan'),3,'omitnan'));
%mm2 = triu(mm2.',1) + tril(mm2);
mmm2 = squeeze(mean(mm2,'omitnan'));

ms1 = squeeze(mean(m1, 3, 'omitnan'));
ms2 = squeeze(mean(m2, 3, 'omitnan'));
msD = ms1-ms2; 
times = -.5:.1:0.9; 

figure(); set(gcf, 'Position', [100 100 500 200])
%plot(times, msD', 'color', [.5 .5 .5]); hold on; 
%plot(mean(msD), 'k', 'LineWidth',2)
mD = mean(msD);
stdD = std(msD);
seD = stdD / sqrt(26);
shadedErrorBar(times, mD,seD, 'k', 1); hold on 


[h p ci ts] = ttest(msD)
t = ts.tstat
clustinfo = bwconncomp(h);
h(h==0) = nan; h(h==1) = 0;
plot(times, h, 'r', 'LineWidth', 2)
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
set(gca, 'FontSize', 14)

exportgraphics(gcf, 'myPlot.png', 'Resolution', 300)

%% permutations 
nPerm = 1000; 

for permi = 1:nPerm
    
    for subji = 1:size(m1, 1)
        if rand<.5
            m1F(subji, :, :) = m1(subji, : ,:); 
            m2F(subji, :, :) = m2(subji, : ,:); 
        else
            m1F(subji, :, :) = m2(subji, : ,:); 
            m2F(subji, :, :) = m1(subji, : ,:); 
        end
        ms1 = squeeze(mean(m1F(:,6:15,:), 2, 'omitnan'));
        ms2 = squeeze(mean(m2F(:,6:15,:), 2, 'omitnan'));        
        msD = ms1-ms2; 
    end

    [hPerm p ci ts] = ttest(msD);
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


%% 
figure()
histogram(max_clust_sum)

%% permutations
n_perm = 1000;

%allAb = max_clust_sum(abs(max_clust_sum) > obsT);
allAb = max_clust_sum(max_clust_sum > obs);
p =1 - (n_perm - (length (allAb)+1) )  /n_perm;

disp (['p = ' num2str(p)]);
%x = find(max_clust_sum_ranked(:,2) == index)



%% 
figure()
imagesc(h)

%% 
mD = mean(msD);
stdD = std(msD);
seD = stdD / sqrt(26);
figure()
shadedErrorBar(1:45, mD,seD, 'r', 1); hold on 

%plot(h, 'LineWidth', 2)




%% LOAD all conditions IN LOOP

clearvars -except region
paths = load_paths_WM(region); 
currentDir = pwd; 

%cd (paths.results.bands)

disp(string(datetime));
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);

contrasts = {   
                  'SISC_EM2' 'DISC_EM2'; ...
%                 'DISC_EM2UV1' 'DIDC_EM2UV1'; ...
%                 'SISC_EM2UV2' 'DISC_EM2UV2'; ...
%                 'DISC_EM2UV2' 'DIDC_EM2UV2'; ...
%                 'DISC_M2123NC' 'DIDC_M2123NC'; ...
             };
% % use with 3 trials
cmaps2use = {[-.02 .02] [-.02 .02] [-.02 .02] [-.01 .015] [-.02 .02] ...
             [-.015 .015] [-.015 .015] [-.0125 .0125] [-.005 .0075] [-.025 .025] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.025 .025] [-.025 .025] [-.025 .025] [-.01 .015] [-.025 .025] ...
             [-.02 .02] [-.02 .02] [-.02 .02] [-.01 .015] [-.02 .02] ...
             [-.015 .015] [-.015 .015] [-.0125 .0125] [-.005 .0075] [-.025 .025] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.025 .025] [-.025 .025] [-.025 .025] [-.01 .015] [-.025 .025] ...
             [-.02 .02] [-.02 .02] [-.02 .02] [-.01 .015] [-.02 .02] ...
             [-.015 .015] [-.015 .015] [-.0125 .0125] [-.005 .0075] [-.025 .025] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.01 .015] [-.01 .015] [-.01 .01] [-.01 .01] [-.02 .02] ...
             [-.025 .025] [-.025 .025] [-.025 .025] [-.01 .015] [-.025 .025] ...
             };
   

perms2use = {   '1-4' ...
                '1-4' ...
                '1-4' ...
                '4-4' ...
                };
t2use = {'100_norm' '100_perm'};



cmapi = 1;
for foldi = 3:length(fold) %start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym
    
    direct = fold(foldi);
    cd (direct.name)
    %processFoldersWM;    
    
        
    c = unique (contrasts);
    d = cellfun(@(x) [x '_id'], c, 'un', 0);
    for i = 1:length(c) 
        load([c{i} '.mat']);
        contrData{i,:} = eval(c{i});
        idData{i,:} = all_IDs;
        %idData{i,:} = [];
    end
    
    noAv = 1;
    [out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv);
    for i = 1:length(out_c) 
        eval([c{i} ' = out_c{i};']);
        %eval([d{i} ' = out_id{i};']);
    end




    for permNoperm = 1:1 %1 = just plots (no perm) (2:2) just permutation
        for imi = 1:size(contrasts,1)

            clear all_cond1 all_cond2 all_cond1_A all_cond2_A;

            cond1 = contrasts{imi, 1};
            cond2 = contrasts{imi, 2};
            



            all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);

            %global parameters
            if strcmp(region, 'pfc')
                subj2exc        =       [1];
            elseif strcmp(region, 'vvs')
                 subj2exc        =       [18 22];
            else
                 subj2exc        =       [];
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
            cfg.climT       =       [-7 7]; %color scale for t-map
            plotClust       =       1; 
            dupSym          =       0; 
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
            cfg_plot.lwd1           = cfg.lwd1; 
            cfg_plot.lwd2           = cfg.lwd2; 
            cfg_plot.out_real       = out_real;
            cfg_plot.dupSym         = dupSym; 
            cfg_plot.plotD          = 0; 
            cfg_plot.plotClust      = plotClust; 


            if strcmp(test2use, 'wilcox')
                cfg_plot.diff     = 'diff'; % or 'tmap';
            else
                cfg_plot.diff     = 'tmap'; % or 'tmap';
            end

            [myPlotR] = plot_reinst_map_wm(cfg_plot);

            if runperm & ~isempty (out_perm.max_clust_sum_real)
                max_clust_sum = out_perm.max_clust_sum;
                
                %first find 
                obs = max(out_real.all_clust_tsum_real(:,1));
                allAb = max_clust_sum(abs(max_clust_sum) > obs);
                p =1 - (n_perm - (length (allAb)-1) )  /n_perm;
                disp (['p = ' num2str(p)]);
                %x = find(max_clust_sum_ranked(:,2) == index)
            end

            
        end
    end 
cd ..
end


disp ('all plots done');
disp(string(datetime));


%% Collect all diagonals in loop
clear 

cond1 = 'DISC_EE';
cond2 = 'DIDC_EE';

subj2exc = 1; 

disp(string(datetime));
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);

for foldi = 3:length(fold) %start at 3 cause 1 and 2 are . and ...
    
    direct = fold(foldi);
    cd (direct.name)
    load (cond1)
    load (cond2)

    all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);
    %exclude subjects
    if subj2exc > 0
        all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
        all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
    end

    c1(foldi-2,:,:) = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond1_A, 'un', 0));
    c2(foldi-2,:,:) = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond2_A, 'un', 0));
    
    cd ..
end

%% 
mc1 = squeeze(mean(c1, 2));
mc2 = squeeze(mean(c2, 2));
mcD = mc1-mc2;

times = -.25:.1:1.75;
figure()
plot(times, mcD', 'LineWidth', 2); hold on; 
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'),[0 0], 'k:', 'LineWidth', 2);
legend({'13-29Hz', '3-8Hz', '3-150Hz', '30-75Hz', '75-150Hz', '9-12Hz'})


%% permutations
n_perm = 1000;
max_clust_sum = out_perm.max_clust_sum;
%obs = max(abs(all_clust_tsum_real(:,1)));
obs = abs(out_perm.max_clust_sum_real); 

%allAb = max_clust_sum(abs(max_clust_sum) > obs);
allAb = max_clust_sum(max_clust_sum > obs);
p =1 - (n_perm - (length (allAb)+1) )  /n_perm;

 
disp (['p = ' num2str(p)]);
%x = find(max_clust_sum_ranked(:,2) == index)

 
 

%% plot histogram
t1 = [263.022340524431]
t2 = [154.218956436873]
%t1 =  16.7624423475231   % 17.2184940279255 % 16.7624423475231  % [118.918634126354];
[counts,centers] = hist(out_perm.max_clust_sum, 13);

figure(1); set(gcf, 'Position', [100 100 410 400])
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 25, 'ylim', [0 320] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 5,'Color', [.5 .0 .0] );
plot ([t2 t2], get(gca, 'ylim'), 'LineWidth', 5,'Color',  [.3 .3 .3] );
filename = ['histogram.png'];
export_fig(1, filename,'-transparent', '-r300');
%close all;



%% MEAN IN CLUSTER OF INTEREST
 
subj2exc    = [1]; 

rsa_cond1 = cell2mat(cellfun(@mean, DISC_EM2, 'un',0)); 
rsa_cond2 = cell2mat(cellfun(@mean, DIDC_EM2, 'un',0)); 

plotD           = 0;
takeclust       = 1;

roiX = 1:5; roiY = 13:25; %paper = 13-25
mlimE = 1:15; mlimR = 1:45;

% roiX = 75:135; roiY = 61:125; %used in 10ms analysis
% mlimE = 26:225; mlimR = 26:225;%no need to substract 250 to make it at th emiddle
% bins = 200; mfL = bins/4; 
 
pixel = 2;
%exclude subjects
if subj2exc > 0
    rsa_cond1(subj2exc, :, :) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    rsa_cond2(subj2exc, :, :) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end
 




%mlimR = 16:45;

rsa_cond1 = rsa_cond1(:,mlimE,mlimR);
rsa_cond2 = rsa_cond2(:,mlimE,mlimR);
 
nSubj = size(rsa_cond1, 1);
if takeclust
    ROI_cond1 = rsa_cond1(:,clustInfoReal.PixelIdxList{pixel}); 
    ROI_cond2 = rsa_cond2(:,clustInfoReal.PixelIdxList{pixel}); 
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
   
else
    ROI_cond1 = rsa_cond1(:,roiX, roiY); 
    ROI_cond2 = rsa_cond2(:,roiX, roiY); 
    ROI_cond1 = reshape (ROI_cond1, [nSubj size(roiX,2) * size(roiY, 2)]);
    ROI_cond2 = reshape (ROI_cond2, [nSubj size(roiX,2) * size(roiY, 2)]);
    m_roi_cond1 = mean(ROI_cond1, 2);
    m_roi_cond2 = mean(ROI_cond2, 2);
end
 
test = squeeze(mean(rsa_cond1, 1));
test(1:end) = 0;
test1 = test;
test(roiX,roiY,:) = .5;


if takeclust
    test(clustInfoReal.PixelIdxList{pixel}) = 1;
    test1(clustInfoReal.PixelIdxList{pixel}) = 1;
elseif takeoverlap
    test(clustInfoOverlapping.PixelIdxList{1}) = 1;
end
diag1 = diag(test);
if takeclust
    test(clustInfoReal.PixelIdxList{pixel}) = 1;
end
%test = flipud(test);
figure(1);imagesc(test);hold on; %axis square;hold on;
%contour(test, 1, 'lineWidth', 4, 'linecolor', 'k');axis square; %does not work with plotD

 
 
% if plotD
%  plot([mfL+0.5 mfL+0.5],get(gca,'ylim'),'w', 'LineWidth', 2); 
%  %plot([bins - (mfL-0.5) bins - (mfL-0.5)], get(gca,'xlim'), 'w', 'LineWidth', 2);
%  plot(get(gca,'xlim'), [bins - (mfL-0.5) bins - (mfL-0.5)],'w', 'LineWidth', 2);
%  plot(get(gca,'xlim'), [bins+.5 .5],'w', 'LineWidth', 2); % diagonal
% end 
% set (gca, 'clim', [0 1]);
% %remove labels
% axesHandles = findall(0, 'type', 'axes');
% for i=1:length(axesHandles)
%    set (axesHandles(i), 'visible', 'off'); % remove labels 
% end
 
 
%[h p ci t] = ttest (m_roi_cond1, m_roi_cond2);
%disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
 
[p,h,stats] = signrank(m_roi_cond1, m_roi_cond2);
disp (['W = ' num2str(stats.signedrank) '  ' ' p = ' num2str(p)]);

 
diff_ERS = m_roi_cond1-m_roi_cond2;
 
%set (gcf, 'Position', [200 200 400 450]);
%export_fig(2, '_cluster.png','-transparent', '-r80');
%close all;   
 


%% 2Bar 
%data.data = [data_LC.data(:,1) data_LC.data(:,2) ];
%data.data = [diff_ERS_all diff_ERS_all_minus]; 
data.data = [m_roi_cond1 m_roi_cond2]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(13); 
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.01 .015] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

%[h p ci t] = ttest (data.data(:,1), data.data(:,2));
[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

%export_fig(2, '_2.png','-transparent', '-r80');
%close all;   





























%% 




















%% RSA GLOBAL CREATE FOLDERS
% run in folder with gOBO_rsa files from all subj together
% first create the folders based on the conditions's names
clear, close all
sublist = dir('*_rsa.mat'); sublist = {sublist.name};
fname_tmp = 's0';
for filei=1:length(sublist)
   str = sublist{filei};
   fname = str(5:end-13);
   %fname = str(5:end-7);
   if ~strcmp(fname, fname_tmp)
      mkdir(fname)
      movefile(str, fname)
   else
      movefile(str, fname)
   end
   str_tmp = sublist{filei};
   fname_tmp = str_tmp (5:end-13);
   %fname_tmp = str_tmp (5:end-7);
end

%%GO THROUGH FOLDERS RECURSIVELY
% run in folder with subfolders SI_C, SI_F... SI_R
clear, close all
tic
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);
 
 
for foldi = 3:length(folders) %start at 3 cause 1 and 2 are . and ...
    direct = folders(foldi);
    cd (direct.name)
    sublist = dir('*_rsa.mat');
    sublist = {sublist.name};
    disp (['measurements -> ' num2str(length(sublist))]);
 
    for subji=1:length(sublist)
        load(sublist{subji});
        %chan2plot(subji)
        all{subji,1} = rsaZ; %when more than 1 electrode
        %all{subji,1} = rsaZ;
        if exist('allIDs')
            all_IDs{subji, :} = allIDs; 
        end
 
    end
 
    if exist ('timeBins')
        timeBins1 = timeBins (:, [1, end]);
    end
 
 
    cd .. % goes up one directory
    filename = [sublist{subji}(5:end-13)];
    %filename = [sublist{subji}(5:end-12)];
    eval([filename '= all;']);
    save (filename, filename, '-v7.3');
    %save (filename, filename);
    
    
end

toc

%averageSub_WM 
 

%% test average
clearvars

contrasts = {%'SISC_EE' 'DISC_EE'; ... 
             %'DISC_EE' 'DIDC_EE'; ...
             %'SISC_EM2' 'DISC_EM2'; ... 
             %'DISC_EM2' 'DIDC_EM2'; ...
             %'DISC_M2M2' 'DIDC_M2M2';...
             
             };


c = unique (contrasts);
for i = 1:length(c) 
    load([c{i} '.mat']);
    contrData{i,:} = eval(c{i});
end

out_contrasts = averageSub_WM (c, contrData, 1000)
for i = 1:length(out_contrasts) 
    eval([c{i} ' = out_contrasts{i};']);
end


clearvars -except DISC_EE DIDC_EE
all_cond1 = DISC_EE;
all_cond2 = DIDC_EE;

clear junts labels meanReal_cond1 meanReal_cond2;
for si = 1:length(all_cond1)
    junts = cat(1, all_cond1{si}, all_cond2{si});
    labels = [zeros(1, size(all_cond1{si}, 1))  ones(1, size(all_cond2{si}, 1))];
    %real labels
    meanReal_cond1(si, :, :) = squeeze(mean(junts(labels==0, :, :), 1, 'omitnan')); 
    meanReal_cond2(si, :, :) = squeeze(mean(junts(labels==1, :, :), 1, 'omitnan')); 
end

clear fin_cond1
for si = 1:length(all_cond1)
    si
    d2u = all_cond1{si};
    for tri = 1:size(d2u, 1)
        d2u(tri,:,:) = d2u(tri,:,:) - meanReal_cond2(si, :, :) ;
    end
    fin_cond1(si,:,:) = mean(d2u);
end



figure();
subplot(211)
d2p = squeeze(mean(fin_cond1)); 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('cond1')
lim1 = get(gca, 'clim');

subplot(212)
[h p ci t] = ttest(fin_cond1);
h = squeeze(h);t = squeeze(t.tstat);
d2p = h; 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('stats')


%% 
for subji = 1:74
    d2p = squeeze(mean(DISC_EE{subji}, 1));
    figure()
    imagesc(d2p); colorbar
end

%% 
figure()
%d2p = squeeze(rsaZ(1,:,:)); 
d2p = squeeze(mean(rsaZ,1, 'omitnan')); 
imagesc(d2p);colorbar


%%
figure();
subplot(421)
d2p = squeeze(mean(fin_cond1)); 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('cond1')
lim1 = get(gca, 'clim');

subplot(422)
d2p = squeeze(mean(meanReal_cond2)); 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('cond2')
set(gca, 'clim', lim1);

subplot(423)
diff3 = meanReal_cond1 - meanReal_cond2; 
d2p = squeeze(mean(diff3)); 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('diff')

subplot(424)
[h p ci t] = ttest(meanReal_cond1);
h = squeeze(h);t = squeeze(t.tstat);
d2p = h; 
imagesc(d2p);hold on; colorbar; axis square;
set (gca, 'xlim', [6 15], 'ylim', [6 15]);
title('stats')


