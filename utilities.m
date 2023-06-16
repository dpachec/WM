%% Plot all trials
%% load cfg_contrasts
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




































%% UTILITIES FROM pattern similarity script
% % % % % % % % 

%% Plot all trials
%%load cfg_contrasts
clear
tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:1%length(sublist)
    disp(['File > ' num2str(subji)]);
    load(sublist{subji});
    s_all{subji} = cfg_contrasts; 
    %size(cfg_contrasts.oneListPow)
end
s_all = s_all';
toc
 
%% plot by trials
bline               = [100 200];
acrossTrials       = 0;
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
           %export_fig(pagi, filename,'-r200');    
           close all;  
       end
    end
   
end
 

%% generate random vectors for plots

x = randn(1, 20)
figure(); set(gcf, 'Position', [50 100 575 25])
myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)

imagesc(x)
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'clim', [-1.6 1.6])
export_fig(gca, 'hola.png','-r300');  


%% 


x = randn(1, 75)
figure(); set(gcf, 'Position', [50 100 150 25])
myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)

imagesc(x)
set(gca, 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'clim', [-1.6 1.6])
export_fig(gca, 'hola.png','-r300');  


















 

%% plot by channels
bline               = [3 7];
acrossTrials       = 1;
zScType = 'allTrials'; 
 
for subji = 1:1 %1:length(sublist)
   cfg_contrasts = s_all{subji};    
   cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
   oneListPow = cfg_contrasts.oneListPow;
   for chani = 1:1%size(oneListPow, 2)
       trialN = size(oneListPow, 1);
       pag2plot = ceil(trialN/100);
       cstri = 0;
       for pagi = 1:1%pag2plot
           figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
           csubi = 1;
           startTri = (cstri*100) +1
           endTri = (cstri*100) +100
           if endTri > trialN
               endTri = trialN
           end
           
           for triali = startTri:endTri
                subplot (10, 10, csubi)
                myCmap = colormap(brewermap([],'*Spectral'));
                colormap(myCmap)
                d2p = squeeze(oneListPow(triali,chani,1:30,:)); 
                %imagesc(); 
                contourf(d2p, 40, 'linecolor', 'none')
                %title(num2str(triali))
                set(gca,'YDir','normal', 'clim', [-1.5 5], 'xlim', [200 400], 'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [])
                %colorbar
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
clearvars
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

elLength = cell2mat(cellfun(@length, s_all, 'un', 0))
sum(elLength)

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
 
 
%% ALL correlation analysis ENCODING
clear 
% loop through all folders with band specific results 
paths = load_paths_WM('pfc');
main_path_pfc = paths.results.bands; 
paths = load_paths_WM('vvs');
main_path_vvs = paths.results.bands;  

%pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
%vvs_fits = nnFit([7 9 13 18 19 20 21 23 27 28]); 

fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

for bandi = 1:6

    flist = dir([main_path_vvs fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_vvs fnames{bandi} '\' f2load])
    %ALL_EE = ALL_EE([7 9 13 18 19 20 21 23 27 28]); 
    %allPFC{bandi,:} = ALL_EE; % 6 bands and 10 subjects
    ALL_EM2 = ALL_EM2([7 9 13 18 19 20 21 23 27 28]); 
    allPFC{bandi,:} = ALL_EM2; % 6 bands and 10 subjects

    flist = dir([main_path_pfc fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_pfc fnames{bandi} '\' f2load])
    %ALL_EE = ALL_EE([2 3  5  9 10 11 12 14 15 16]);
    %allVVS{bandi,:} = ALL_EE;  % 6 bands and 10 subjects
    ALL_EM2 = ALL_EM2([2 3  5  9 10 11 12 14 15 16]);
    allVVS{bandi,:} = ALL_EM2;  % 6 bands and 10 subjects


end


 

%% 
clear c1VVS allPSPFC allPSVVS clear allRho
onlyDiag = 0; % is the pattern a 2D matrix or just matching time points
timeL = 6:15; 

for bandi = 1:6

    c1VVS_band = allVVS{bandi}; 
    c1PFC_band = allPFC{bandi}; 

    for subji = 1:10
    
        c1VVS_subj = c1VVS_band{subji}; 
        c1PFC_subj = c1PFC_band{subji}; 
        
        clear cVVSF cPFCF
        for triali = 1:size(c1VVS_subj, 1)
            
            cVVS_tr = squeeze(c1VVS_subj(triali, :, :));
            cPFC_tr = squeeze(c1PFC_subj(triali, :, :));

            % % % only diagonal (if matrix)
            if ~onlyDiag
                cVVS_tr_diag = diag(cVVS_tr(timeL, timeL)); % restricted to the period of significant dynamicity
                cPFC_tr_diag = diag(cPFC_tr(timeL, timeL)); % restricted to the period of significant dynamicity
            else %diagonal was extracted before
                cVVS_tr_diag = cVVS_tr(timeL); % restricted to the period of significant dynamicity
                cPFC_tr_diag = cPFC_tr(timeL); % restricted to the period of significant dynamicity
            end
           
            cVVSF(triali, :) = mean(cVVS_tr_diag); 
            cPFCF(triali, :) = mean(cPFC_tr_diag); 
    
        end
    
        allRho(bandi, subji,:) = corr(cVVSF, cPFCF, 'type', 's'); 
        

    end

    

end





%% 

[h p ci ts] = ttest(allRho')


%% ALL correlation analysis MAINTENANCE
clear 
% loop through all folders with band specific results 
paths = load_paths_WM('pfc');
main_path_pfc = paths.results.bands; 
paths = load_paths_WM('vvs');
main_path_vvs = paths.results.bands;  

%pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
%vvs_fits = nnFit([7 9 13 18 19 20 21 23 27 28]); 

fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

for bandi = 1:6

    flist = dir([main_path_vvs fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_vvs fnames{bandi} '\' f2load])
    ALL_EM2 = ALL_EM2([7 9 13 18 19 20 21 23 27 28]); 
    allPFC{bandi,:} = ALL_EM2; % 6 bands and 10 subjects

    flist = dir([main_path_pfc fnames{bandi}]); f2load = [flist.name]; f2load = f2load(4:end);
    load([main_path_pfc fnames{bandi} '\' f2load])
    ALL_EM2 = ALL_EM2([2 3  5  9 10 11 12 14 15 16]);
    allVVS{bandi,:} = ALL_EM2;  % 6 bands and 10 subjects


end


 

%% 
clear c1VVS allPSPFC allPSVVS clear allRho
onlyDiag = 0; % is the pattern a 2D matrix or just matching time points
timeL1 = 6:15; 
timeL2 = 6:15; 

for bandi = 1:6

    c1VVS_band = allVVS{bandi}; 
    c1PFC_band = allPFC{bandi}; 

    for subji = 1:10
    
        c1VVS_subj = c1VVS_band{subji}; 
        c1PFC_subj = c1PFC_band{subji}; 
        
        clear cVVSF cPFCF
        for triali = 1:size(c1VVS_subj, 1)
            
            cVVS_tr = squeeze(c1VVS_subj(triali, timeL1));
            cPFC_tr = squeeze(c1PFC_subj(triali, timeL1));

            cVVSF(triali, :) = mean(cVVS_tr, 'all'); 
            cPFCF(triali, :) = mean(cPFC_tr, 'all'); 
    
        end
    
        allRho(bandi, subji,:) = corr(cVVSF, cPFCF, 'type', 's'); 
        

    end

    

end





%%  utilities from GroupFX











%%
clearvars

region = 'vvs'; 
paths = load_paths_WM(region, []);

contrasts = {
              %'SISC_EEC' 'DISC_EEC';
              %'DISC_EE' 'DIDC_EE';
              'DISC_EM2' 'DIDC_EM2';
              %'DISC_M2M2' 'DIDC_M2M2';
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



%% plot simple

subj2exc = [18 22]; 
cond1 = 'DISC_EM2';
cond2 = 'DIDC_EM2';
lim1 = [3:17]; lim2 = [3:47];

all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);


%exclude subjects
if subj2exc > 0
    all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end
cond1 = cellfun(@(x) squeeze(mean(x)), all_cond1_A, 'un', 0);
c1 = cat(3, cond1{:}); c1 = permute(c1, [3 1 2]);c1 = c1 (: , lim1, lim2); 

cond2 = cellfun(@(x) squeeze(mean(x)), all_cond2_A, 'un', 0);
c2 = cat(3, cond2{:}); c2 = permute(c2, [3 1 2]);c2 = c2 (: , lim1, lim2); 

[h p ci ts] = ttest(c1, c2); 
h = squeeze(h); 
t = squeeze(ts.tstat); 
% 
% %remove non-sig clusters
% clustinfo = bwconncomp(h);
%  for pixi = 1:length(clustinfo.PixelIdxList)
%     h(clustinfo.PixelIdxList{pixi}) = 0;   
%  end
%  h(clustinfo.PixelIdxList{4}) = 1; % VVS
%  h(clustinfo.PixelIdxList{10}) = 1; % VVS
% % nothing in PFC

c1 = squeeze(mean(c1(:, 6:15, :), 2));
c2 = squeeze(mean(c2(:, 6:15, :), 2));
d = c1-c2; 
[h2 p2 ci2 ts2] = ttest(d); 

% %remove non-sig clusters
% clustinfo = bwconncomp(h2);
% for pixi = 1:length(clustinfo.PixelIdxList)
% h2(clustinfo.PixelIdxList{pixi}) = 0;   
% end
% h2(clustinfo.PixelIdxList{2}) = 1; % vvs
% h2(h2==0) = nan; h2(h2 == 1) = 0; 
% t2 = ts2.tstat ;
% h2Obs = sum(t2(clustinfo.PixelIdxList{2})); 

md = mean(d); 
stdd = std(d); 
sed = stdd/sqrt(size(c1, 1)); 


colormap(brewermap([],'YlOrRd')); 
%times = -.5:.1:3.99;
times = 1:45;
set(gcf, 'Position', [100 100 500 300])
tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'none'); 
nexttile
%shadedErrorBar(times, md,sed, 'k', 1); hold on 
errorbar(times, md, sed, 'k', 'LineWidth', 3); hold on 
%set(gca, 'xlim', [.5 45.5], 'ylim', [-.001 .003], 'xtick',[], 'ytick', []) % VVS
set(gca, 'xlim', [.5 45.5], 'ylim', [-.001 .003], 'xtick',[], 'ytick', []) % PFC
plot([4.95 4.95], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot([9.95 9.95], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(times, h2, 'Color', [.9 .7 .1], 'LineWidth', 10)
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 2);
nexttile
contourf(1:450, 1:150, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(1:450, 1:150, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
plot([45 45], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot([95 95], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'), [45 45], 'k:', 'LineWidth', 2);
set(gca,'xtick',[], 'ytick', [], 'clim', [-3 7]); %colorbar

exportgraphics(gcf, 'myP.png', 'Resolution', 150);


%% permutations 
nPerm = 1000; 

for permi = 1:nPerm
    
    for subji = 1:size(c1, 1)
        if rand<.5
            m1F(subji, :) = c1(subji,:); 
            m2F(subji, :) = c2(subji,:); 
        else
            m1F(subji, :) = c2(subji,:); 
            m2F(subji, :) = c1(subji,:); 
        end 
    end
    d = m1F-m2F; 
    [hPerm p ci ts] = ttest(d);
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


%% permutations
n_perm = 1000; 

%first find 
obs = sum(t(clustinfo.PixelIdxList{6}));

allAb = max_clust_sum(abs(max_clust_sum) > obs);
%allAb = max_clust_sum((max_clust_sum) > obs);
p = (1 - (n_perm - (length (allAb) ) )  /n_perm) + (1/ n_perm);
disp (['p = ' num2str(p)]);
%x = find(max_clust_sum_ranked(:,2) == index)

 

%% plot diagonal only


cond1 = 'DISC_M2M2';
cond2 = 'DIDC_M2M2';

subj2exc = [1]; 

all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);
%exclude subjects
if subj2exc > 0
    all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end

%% extract full map 
cond1 = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond1_A, 'un', 0));
cond2 = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), all_cond2_A, 'un', 0));

%% plot 

d2p1 = squeeze(mean(cond1))
d2p2 = squeeze(mean(cond2))

mDiff = d2p1-d2p2; 

figure; contourf(mDiff); hold on
h = ttest(cond1, cond2); 
%contour(squeeze(h))










%% average in period 
clear m1 m2
for subji = 1:length(all_cond1)
    %m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,6:15,6:45)));
    %m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,6:15,6:45)));

    lim2use = 30:40
    m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,6:40, lim2use)));
    m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,6:40, lim2use))); %for the 100_norm cut out

end

m1A = squeeze(mean(mean(m1,2,'omitnan'),3,'omitnan'));
m2A = squeeze(mean(mean(m2,2,'omitnan'),3,'omitnan'));


%%full maintenance period 
msD = m1A-m2A;
[h p ci ts] = ttest(msD);
t = ts.tstat
p


%% average in period (only one time dimension)
clear m1 m2
for subji = 1:length(all_cond1)
    m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,1:40,1:40)));
    m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,1:40,1:40)));


end

m1A = squeeze(mean(m1,2,'omitnan'));
m2A = squeeze(mean(m2,2,'omitnan'));


%%full maintenance period 
msD = m1A-m2A;
[h p ci ts] = ttest(msD);
t = ts.tstat
p

plot(mean(msD))


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
    'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.045 .045] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

%[h p ci t] = ttest (data.data(:,1), data.data(:,2));
[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r80');
close all;   

%% permutations
n_perm = 1000; 

max_clust_sum = out_perm.max_clust_sum;

%maintenance clusters : 192.369408876301 130.459520724292
%first find 
obs = max(all_clust_tsum_real(:,1));

allAb = max_clust_sum((max_clust_sum) > obs);
p = (1 - (n_perm - (length (allAb) ) )  /n_perm) + (1/ n_perm);
disp (['p = ' num2str(p)]);
%x = find(max_clust_sum_ranked(:,2) == index)

 

%% plot simple

subj2exc = [2]; 
cond1 = 'DISC_EM2';
cond2 = 'DIDC_EM2';
lim1 = [3:17]; lim2 = [3:47];

all_cond1_A = eval(cond1);all_cond2_A = eval(cond2);


%exclude subjects
if subj2exc > 0
    all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end
cond1 = cellfun(@(x) squeeze(mean(x)), all_cond1_A, 'un', 0);
c1 = cat(3, cond1{:}); c1 = permute(c1, [3 1 2]);c1 = c1 (: , lim1, lim2); 

cond2 = cellfun(@(x) squeeze(mean(x)), all_cond2_A, 'un', 0);
c2 = cat(3, cond2{:}); c2 = permute(c2, [3 1 2]);c2 = c2 (: , lim1, lim2); 

[h p ci ts] = ttest(c1, c2); 
h = squeeze(h); 
t = squeeze(ts.tstat); 

%remove non-sig clusters
clustinfo = bwconncomp(h);
for pixi = 1:length(clustinfo.PixelIdxList)
   h(clustinfo.PixelIdxList{pixi}) = 0;   
end
h(clustinfo.PixelIdxList{4}) = 1;
h(clustinfo.PixelIdxList{10}) = 1;

c1 = squeeze(mean(c1(:, 6:15, :), 2));
c2 = squeeze(mean(c2(:, 6:15, :), 2));
d = c1-c2; 
[h2 p2 ci2 ts2] = ttest(d); 
%remove non-sig clusters
clustinfo = bwconncomp(h2);
for pixi = 1:length(clustinfo.PixelIdxList)
   h2(clustinfo.PixelIdxList{pixi}) = 0;   
end
h2(clustinfo.PixelIdxList{2}) = 1;
h2(h2==0) = nan; h2(h2 == 1) = 0; 
t2 = ts2.tstat ;


md = mean(d); 
stdd = std(d); 
sed = stdd/sqrt(size(d, 1)); 


colormap(brewermap([],'YlOrRd')); 
%times = -.5:.1:3.99;
times = 1:45;
set(gcf, 'Position', [100 100 500 300])
tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'none'); 
nexttile
%shadedErrorBar(times, md,sed, 'k', 1); hold on 
errorbar(times, md, sed, 'k', 'LineWidth', 3); hold on 
plot([4.95 4.95], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(times, h2, 'LineWidth', 4)
set(gca, 'xlim', [.5 45.5], 'xtick',[], 'ytick', [])
nexttile
contourf(1:450, 1:150, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(1:450, 1:150, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
plot([45 45], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'), [50 50], 'k:', 'LineWidth', 2);
set(gca,'xtick',[], 'ytick', [], 'clim', [-3 7]); %colorbar

exportgraphics(gcf, 'myP.png', 'Resolution', 150)



%%
figure()
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'none'); 
for i=1:4    
    nexttile    
    plot(rand(1,10)); 
    set(gca,'xtick',[], 'ytick', [])
end





%% encoding - maintenance similarity (average full maintenance period)
ms1 = squeeze(mean(m1, 2, 'omitnan'));
ms2 = squeeze(mean(m2, 2, 'omitnan'));
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
   m1(subji, :,:) = squeeze(mean(all_cond1{subji}(:,6:15,1:45)));
   m2(subji, :,:) = squeeze(mean(all_cond2{subji}(:,6:15,1:45)));
end

mm1 = squeeze(mean(m1,'omitnan'));
m1A = squeeze(mean(mean(m1,2,'omitnan'),3,'omitnan'));
%mm1 = triu(mm1.',1) + tril(mm1);
mmm1 = squeeze(mean(mm1,'omitnan'));
mm2 = squeeze(mean(m2,'omitnan'));
m2A = squeeze(mean(mean(m2,2,'omitnan'),3,'omitnan'));
%mm2 = triu(mm2.',1) + tril(mm2);
mmm2 = squeeze(mean(mm2,'omitnan'));

ms1 = squeeze(mean(m1, 2, 'omitnan'));
ms2 = squeeze(mean(m2, 2, 'omitnan'));
msD = ms1-ms2; 
times = -.5:.1:3.9; 

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
set(gca, 'ylim', [-.0025 .005])
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
set(gca, 'FontSize', 18, 'xlim', [-.5 3.5])

exportgraphics(gcf, 'myPlot.png', 'Resolution', 300)

%% encoding - maintenance similarity (average full maintenance period)
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
                  'SISC_EE' 'DISC_EE'; ...
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




    for permNoperm = 2:2 %1 = just plots (no perm) (2:2) just permutation
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

obs = 27.4464968762709

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

data = squeeze(cfg_contrasts.oneListPow(1,1,:,:));

figure
contourf(data,  40, 'linecolor', 'none')

















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


%% plot all frequnecies maintenance

clear, clc
paths = load_paths_WM('vvs');
main_path = paths.results.bands; 

fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';

for bandi = 1:6

    flist = dir([main_path fnames{bandi}]); 
    load([main_path fnames{bandi} '\DISC_EM2.mat'])
    load([main_path fnames{bandi} '\DIDC_EM2.mat'])
    allVVS{bandi,1} = DISC_EM2; 
    allVVS{bandi,2} = DIDC_EM2; 

end




%%

subj2exc = [18 22]; 
lim1 = [3:17]; lim2 = [3:47];

clear allH2s allT2s allH2s allDs
for bandi = 1:6

    all_cond1_A = allVVS{bandi, 1};all_cond2_A = allVVS{bandi, 2};

    %exclude subjects
    if subj2exc > 0
        all_cond1_A(subj2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
        all_cond2_A(subj2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
    end
    cond1 = cellfun(@(x) squeeze(mean(x)), all_cond1_A, 'un', 0);
    c1 = cat(3, cond1{:}); c1 = permute(c1, [3 1 2]);c1 = c1 (: , lim1, lim2); 
    
    cond2 = cellfun(@(x) squeeze(mean(x)), all_cond2_A, 'un', 0);
    c2 = cat(3, cond2{:}); c2 = permute(c2, [3 1 2]);c2 = c2 (: , lim1, lim2); 
    
    
    c1 = squeeze(mean(c1(:, 6:15, :), 2));
    c2 = squeeze(mean(c2(:, 6:15, :), 2));
    d = c1-c2; 
    [h2 p2 ci2 ts2] = ttest(d); 
    ts2 = ts2.tstat; 

    md = mean(d); 
    stdd = std(d); 
    sed = stdd/sqrt(size(c1, 1)); 
    
    
    allT2s(bandi, :) = ts2; 
    allH2s(bandi, :) = h2; 
    allDs(bandi, :) = md;

% % % %     %colormap(brewermap([],'YlOrRd')); 
% % % %     times = 1:45;
% % % %     %set(gcf, 'Position', [100 100 500 300])
% % % %     %shadedErrorBar(times, md,sed, 'k', 1); hold on 
% % % %     plot(times, md, 'linewidth', 2); hold on; 
% % % %     set(gca, 'xlim', [.5 45.5], 'ylim', [-.002 .005], 'xtick',[], 'ytick', []) % PFC
% % % %     plot([4.95 4.95], get(gca, 'ylim'), 'k:', 'LineWidth', 1);
% % % %     plot([9.95 9.95], get(gca, 'ylim'), 'k:', 'LineWidth', 1);
% % % %     plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 1);
% % % %     %plot(times, h2, 'Color', [.9 .7 .1], 'LineWidth', 10)
    


    


end


%% 

allH2s(1,:) = []; 
allT2s(1,:) = []; 

colormap(brewermap([],'YlOrRd')); 

%imagesc(allT2s(2:end,:)); hold on; 
%imagesc(allH2s(2:end,:))


clustinfo = bwconncomp(allH2s);
allH2s(clustinfo.PixelIdxList{1}) = 0; 
allH2s(clustinfo.PixelIdxList{2}) = 0; 
allH2s(clustinfo.PixelIdxList{4}) = 0; 
allH2s(clustinfo.PixelIdxList{5}) = 0; 
allH2s(clustinfo.PixelIdxList{6}) = 0; 


%times = 1:450;
times = -.5:.01:3.99;
freqs = 1:50; 
contourf(times, freqs, myresizem(allT2s, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(allH2s, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set (gca, 'clim', [-5 5])

% d2p = allDs'; 
% 
% plot(d2p(:, 1:6))



%% utilities POWER ANALYSI
%% plot example trial 

size(allPow{1}.oneListPow)

d2p = squeeze(allPow{1}.oneListPow(1, 1,:,:));

figure
imagesc(d2p)










 