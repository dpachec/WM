%% generate maps for each frequency
%%
clear 
paths = load_paths_WM('pfc');
filelistSess = getFilesWM(paths.out_contrasts_path);


frequncies2test = [1:54]';
for i = 1:length(frequncies2test)
    fnames{i,:} = [num2str(i, '%02.f') 'Hz'];
end

avTime              = 0; 
win_width           = 5; 
mf                  = 1; 
meanInTime          = 1; 
meanInFreq          = 0; 
takeElec            = 0; 
takeFreq            = 0;    
TG                  = 0; %temporal generalization
contr2save          = { 'DISC_EE' 'DIDC_EE'}; %{};
bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 100;
n2s                 = 1000000;
loadSurr            = 0; 
zScType             = 'allTrials'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks
avMeth              = 'pow';  


tic

 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.out_contrasts_path filelistSess{sessi}]);   
    
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
        fname = fnames{freqi};
        mkdir ([paths.results.tf_res fname]);
        cd ([paths.results.tf_res fname]);

        if iscell(frequncies2test)
            f           = frequncies2test{freqi}; 
        else
            f           = frequncies2test(freqi); 
        end
 
        rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, takeElec, takeFreq, idxCH, idxF, sessi, TG, avTime)
 
        cd ..
    end
 
end


clear
paths = load_paths_WM(region); 
currentDir = pwd; 

cd (paths.results_path)
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

%%
clear 
paths = load_paths_WM; 
currentFolder = pwd
cd (paths.results.frequency_resolved); 
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
fold([1 2]) = [];

for foldi = 1:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    direct = fold(foldi);
    cd (direct.name)
    
    sublist = dir('*.mat'); sublist = {sublist.name};

    load(sublist{1});
    load(sublist{2});
    
    eval(['C2 = ' sublist{1}(1:end-4)])
    eval(['C1 = ' sublist{2}(1:end-4)])

    C1 = averSub2(C1, 'pfc');
    C2 = averSub2(C2, 'pfc');

    C1A{foldi,:} = C1; 
    C2A{foldi,:} = C2; 
    cd ..
end

cd (currentFolder)

toc

%% 
c1AM = cell2mat(cellfun(@(x) mean(x, 'omitnan'), C1A, 'un',false))
figure()
plot(c1AM)



%% 
clear U2 U21 DISC_H DIDC_H DISC_M2M2_H DIDC_M2M2_H
nFreqs = length(frequncies2test);
for freqi = 1:nFreqs
     
    U2 = DISC_M2M2_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DISC_M2M2_H(freqi, :)  = U21; % average_xGM(U21, 'vvs');
    
    U2 = DIDC_M2M2_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DIDC_M2M2_H(freqi, :)  = U21;  %average_xGMF(U21, 'vvs');
    
    
end

%% 

sub2exc = [1]; 

diff = DISC_M2M2_H - DIDC_M2M2_H;
diff(:, sub2exc) = [] ;
diff = diff';
[h p ci ts] = ttest(diff); 
h = squeeze(h); t = squeeze(ts.tstat); 

clustinfo = bwconncomp(h);
clear allSTs max_clust_sum_obs
for pxi = 1:length(clustinfo.PixelIdxList)
    allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs = allSTs(id); 

hL = h; hL(hL==0) = nan; hL(hL==1) = 0; 
d2p = mean(diff); 
stdu = std(diff); 
seu = stdu / sqrt(size(diff, 1));
figure()
shadedErrorBar([1:29 30:5:150], d2p, seu, 'r', 1); hold on;
plot (hL, 'linewidth', 4); 
set(gca, 'Fontsize', 18, 'xlim', [.5 155]);
set(gca, 'clim', [-3.5 3.5], 'yLim', [-.065 .075] )

exportgraphics(gcf, 'figure.png', 'Resolution', 300);

%% permutations 

sub2exc = [18 22]; 
nPerm = 1000; 

junts = cat(1, DISC_M2M2_H', DIDC_M2M2_H');
labs  = logical([ones(1,28) zeros(1,28)]);

clear diff allSTs max_clust_sum_perm
for permi = 1:nPerm
    labsP = labs(randperm(length(labs)));
    c1 = junts(labsP,:); 
    c2 = junts(~labsP,:);
    diff = c1-c2; 
    diff(sub2exc,:) = [] ;
    [h p ci ts] = ttest(diff); 
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi, :)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 

    %d2p = mean(diff); 
    %stdu = std(diff); 
    %seu = stdu / sqrt(size(diff, 1));
    %figure()
    %shadedErrorBar(1:nFreqs, d2p, seu, 'r', 1); hold on;
    %hL = h; hL(hL==0) = nan; hL(hL==1) = 0; 
    %plot (hL, 'linewidth', 4); 
    %set(gca, 'Fontsize', 18, 'xlim', [.5 nFreqs + .5]);
    %set(gca, 'clim', [-3.5 3.5], 'yLim', [-.045 .0475] )
    
end


%%
t1 = max_clust_sum_obs; 
[counts,centers] = hist(max_clust_sum_perm, 13);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 700] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
%close all;


%% get rank (for stats)

clear p

mcsR = max_clust_sum_obs; 
mcsP = max_clust_sum_perm; 
allAb = mcsP(abs(mcsP) > abs(mcsR));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% Only Broadband

sub2exc = [22]; 

diff = DISC_M2M2_H - DIDC_M2M2_H;
diff(:, sub2exc) = [] ;
diff = diff';
[h p ci ts] = ttest(diff); 
h = squeeze(h); t = squeeze(ts.tstat); 
d2u = [diff(:,1)]
h = [h(1)]; 
p = [p(1)] 


%% One Bar (all freqs)
data.data = [d2u]; 


figure(2); set(gcf,'Position', [0 0 550 500]); 
mean_S = mean(data.data, 1);
hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', 'o', 'MarkerSize',15);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 ],'XTickLabel',{'', ''}, ...
    'FontSize', 22, 'linew',2, 'xlim', [0 2], 'ylim', [-.1 .13] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   



%% 5 Bar 
data.data = [d2u]; 


figure(2); set(gcf,'Position', [0 0 550 500]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2 3 4 5 ], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',25);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4 5 ],'XTickLabel',{'', ''}, ...
    'FontSize', 22, 'linew',2, 'xlim', [0 6], 'ylim', [-.03 .07] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   


%% 54 Bar 
data.data = [diff]; 
data.data(18,:) = [];


figure(2); set(gcf,'Position', [0 0 1200 1200]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(13); 
hb = plot ([1:54], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',5);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 1);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:54],'XTickLabel',{'', ''}, ...
    'FontSize', 20, 'linew',2, 'xlim', [0 55], 'ylim', [-.02 .05] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   


%%

tic
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);

clear SISC_EE_ALL DISC_EE_ALL DIDC_EE_ALL DISC_M2M2_ALL DIDC_M2M2_ALL DISC_M2M2_ALL DIDC_M2M2_ALL
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    direct = fold(foldi);
    cd (direct.name)
    
    load SISC_EE
    load DISC_EE
    load DIDC_EE
    load DISC_M2M2
    load DIDC_M2M2
    load DISC_M2M2
    load DIDC_M2M2
    
    SISC_EE_ALL{foldi,:} = SISC_EE; 
    DISC_EE_ALL{foldi,:} = DISC_EE; 
    DIDC_EE_ALL{foldi,:} = DIDC_EE; 
    DISC_M2M2_ALL{foldi,:} = DISC_M2M2; 
    DIDC_M2M2_ALL{foldi,:} = DIDC_M2M2; 
    DISC_M2M2_ALL{foldi,:} = DISC_M2M2; 
    DIDC_M2M2_ALL{foldi,:} = DIDC_M2M2; 
    cd ..
end


SISC_EE_ALL([1 2]) = []; 
DISC_EE_ALL([1 2]) = []; 
DIDC_EE_ALL([1 2]) = []; 
DISC_M2M2_ALL([1 2]) = []; 
DIDC_M2M2_ALL([1 2]) = []; 
DISC_M2M2_ALL([1 2]) = []; 
DIDC_M2M2_ALL([1 2]) = []; 

toc



%% 
clearvars -except  SISC_EE_ALL DISC_EE_ALL DIDC_EE_ALL DISC_M2M2_ALL DIDC_M2M2_ALL DISC_M2M2_ALL DIDC_M2M2_ALL
for freqi = 1:7
    U2 = SISC_EE_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    SISC_EE_H(freqi, :,:)  = average_xGM(U21, 'vvs');
    
    U2 = DISC_EE_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DISC_EE_H(freqi, :,:)  = average_xGM(U21, 'vvs');
    
    U2 = DIDC_EE_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DIDC_EE_H(freqi, :,:)  = average_xGM(U21, 'vvs');
    
     U2 = DISC_M2M2_ALL{freqi};
     U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
     DISC_M2M2_H(freqi, :,:)  = average_xGM(U21, 'vvs');
     
     U2 = DIDC_M2M2_ALL{freqi};
     U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
     DIDC_M2M2_H(freqi, :,:)  = average_xGM(U21, 'vvs');
     
    U2 = DISC_M2M2_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DISC_M2M2_H(freqi, :,:)  = average_xGM(U21, 'vvs');
    
    U2 = DIDC_M2M2_ALL{freqi};
    U21 = cell2mat(cellfun(@(x) mean(x, 'omitnan'), U2, 'un',0)); 
    DIDC_M2M2_H(freqi, :,:)  = average_xGM(U21, 'vvs');
    
    
end

%% 

diff = DISC_M2M2_H - DIDC_M2M2_H;
diff = diff';
[h p ci ts] = ttest(diff); 
h = squeeze(h); t = squeeze(ts.tstat); 

d2p = mean(diff) ; 
d2p = t;
figure()
plot(d2p); hold on; 
plot (h)

set(gca, 'clim', [-3.5 3.5])




%% 6 Bar 
data.data = [diff]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(13); 
hb = plot ([1 2 3 4 5 6 7], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',25);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3 4 5 6 7],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 8], 'ylim', [-.01 .03] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   



%% 54 Bar 
data.data = [diff]; 
data.data(18,:) = [];


figure(2); set(gcf,'Position', [0 0 1200 1200]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(13); 
hb = plot ([1:54], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',5);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 1);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:54],'XTickLabel',{'', ''}, ...
    'FontSize', 20, 'linew',2, 'xlim', [0 55], 'ylim', [-.02 .05] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   
















%%






































%%