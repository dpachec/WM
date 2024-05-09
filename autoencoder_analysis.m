%% AUTOENCODER 
%% 
clear, clc
f2sav =  'AE0000N_hipp_E123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
cfg = getParams(f2sav);
sessi = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

subj_ch_fr = 1; 
[ACT_FR layNames] = load_AE_activ(cfg, sessi, subj_ch_fr, paths, cfg.net2load);
subj_ch_fr = 20; 
ACT_CH = load_AE_activ(cfg, sessi, subj_ch_fr, paths, cfg.net2load);

    
%% Representational consistency all layers / time points 

ACT = ACT_FR; 
nLays = size(ACT, 1); 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS_FR = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_FR(nLays+1,:) = nan; allMS_FR(:,nLays+1) = nan; %need this for the pColor below
%allMS_FR = tril(allMS_FR, -1); allMS_FR(allMS_FR==0) = nan; 

ACT = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS_CH = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_CH(nLays+1,:) = nan; allMS_CH(:,nLays+1) = nan; %need this for the pColor below
%allMS_CH = tril(allMS_CH, -1); allMS_CH(allMS_CH==0) = nan; 

allMS2 = cat(3, allMS_FR, allMS_CH);
allMS = squeeze(mean(allMS2, 3))

% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 10);
set(gca, 'xtick',  (1:nLays) +.5, 'ytick', (1:nLays) +.5, 'xticklabels', layNames, 'yticklabels', layNames, 'TickLabelInterpreter','none')%'clim', [0.4 1],
colormap(myCmap)


exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%


ACT = ACT_FR; 
nLays = size(ACT, 1); 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS_FR = corr(act3', 'type', 's'); %allMS = allM.^2;

ACT = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS_CH = corr(act3', 'type', 's'); %allMS = allM.^2;

allMS2 = cat(3, allMS_FR, allMS_CH);
allMS = squeeze(mean(allMS2, 3))

figure(); set(gcf, 'Position', [100 100 800 800])
d2p = 1- allMS;
cols = zeros(nLays, 3); cols(1:10, 1) = 1; cols(11:19, 3) = 1; 
[rdmMDS] = cmdscale(d2p);
o2u = .03; 
scatter(rdmMDS(:,1),rdmMDS(:,2), 1550, cols, '.'); axis square
text(rdmMDS(:,1)+o2u,rdmMDS(:,2),layNames, 'Interpreter','none');
%set(gca, 'xlim', [-.5 1], 'ylim', [-.2 .6])
%set(gca, 'xlim', [-.4 .8], 'ylim', [-.3 .4])
%set(gca, 'xlim', [-.4 .8], 'ylim', [-.3 .5])
set(gca, 'xlim', [-.4 .8], 'ylim', [-.1 .6])

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% RNN all RDMS

ACT = ACT_FR; 
nLays = size(ACT, 1); 
figure()
for layi = 1:nLays
   subplot (5, 5, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [-.5 .5])
   %colorbar
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   title(layNames(layi), 'Interpreter','none')
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = nLays+1;
for rowi = 1:5
    for coli = 1:5
        %set(ha(count),'position',[0+coli/6 0+rowi/6 n n ])
        count = count+-1;
        
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS

ACT = ACT_FR; 
nLays = size(ACT, 1); 
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:nLays
   subplot (5, 5, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   title(layNames(layi), 'Interpreter','none')
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = nLays+1;
for rowi = 1:5
    for coli = 1:5
        %set(ha(count),'position',[0+coli/6 0+rowi/6 n n ])
        count = count-1;
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer 

clc
clearvars -except ACT_FR ACT_CH layNames

%create cateogy model
M = kron(eye(6), ones(10));

% freiburg
ACT = ACT_FR;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCW_FR(layi) = mWithin;
    CCB_FR(layi) = mAcross;
    CCI_FR(layi) = mWithin-mAcross ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_CH(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_CH(:, layi) = tmp;

end

%china
ACT = ACT_CH;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    CCW_CH(layi) = mWithin;
    CCB_CH(layi) = mAcross;
    CCI_CH(layi) = mWithin-mAcross ;
    %CCI_CH(layi) = (mAcross - mWithin) ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_FR(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_FR(:, layi) = tmp;

end

CCI2 = [CCI_FR; CCI_CH]; CCW2 = [CCW_FR; CCW_CH]; CCB2 = [CCB_FR; CCB_CH]; 
CCI = mean(CCI2); CCW = mean(CCW2); CCB = mean(CCB2);

aW = squeeze(mean(cat(3, aW_CH, aW_FR), 3)); 
aB = squeeze(mean(cat(3, aB_CH, aB_FR), 3)); 


figure(); set(gcf, 'Position', [100 100 800 800])
plot ([1:10], CCW([1:10]),'b', 'Linewidth', 3); hold on; %axis square
plot ([12:20], CCW([11:19]), 'b', 'Linewidth', 3); 
plot ([1:10], CCB([1:10]), 'r', 'Linewidth', 3); 
plot ([12:20], CCB([11:19]), 'r', 'Linewidth', 3); 
plot ([1:10], CCI([1:10]), 'k', 'Linewidth', 3); 
plot ([12:20], CCI([11:19]), 'k', 'Linewidth', 3); 
plot([11 11], get(gca, 'ylim'), 'k:', LineWidth=2)
set(gca, 'xlim', [0 21])

%legend({'Within category', 'Between category', 'Within - between'},   'Location','northwest', 'FontSize',16); 
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% compute CCI for each layer INCLUDING CLASSIFICATION HEAD

clc
clearvars -except ACT_FR ACT_CH layNames

%create cateogy model
M = kron(eye(6), ones(10));

% freiburg
ACT = ACT_FR;
for layi = 1:size(ACT, 1)
    rdmMDS = squeeze(ACT(layi, :,:)); 
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCW_FR(layi) = mWithin;
    CCB_FR(layi) = mAcross;
    CCI_FR(layi) = mWithin-mAcross ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_CH(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_CH(:, layi) = tmp;

end

%china
ACT = ACT_CH;
for layi = 1:size(ACT, 1)
    rdmMDS = squeeze(ACT(layi, :,:)); 
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    CCW_CH(layi) = mWithin;
    CCB_CH(layi) = mAcross;
    CCI_CH(layi) = mWithin-mAcross ;
    %CCI_CH(layi) = (mAcross - mWithin) ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_FR(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_FR(:, layi) = tmp;

end

CCI2 = [CCI_FR; CCI_CH]; CCW2 = [CCW_FR; CCW_CH]; CCB2 = [CCB_FR; CCB_CH]; 
CCI = mean(CCI2); CCW = mean(CCW2); CCB = mean(CCB2);

aW = squeeze(mean(cat(3, aW_CH, aW_FR), 3)); 
aB = squeeze(mean(cat(3, aB_CH, aB_FR), 3)); 


figure(); set(gcf, 'Position', [100 100 800 800])
plot ([1:10], CCW([1:10]),'b', 'Linewidth', 3); hold on; %axis square
plot ([12:20], CCW([11:19]), 'b', 'Linewidth', 3); 
plot ([1:10], CCB([1:10]), 'r', 'Linewidth', 3); 
plot ([12:20], CCB([11:19]), 'r', 'Linewidth', 3); 
plot ([1:10], CCI([1:10]), 'k', 'Linewidth', 3); 
plot ([12:20], CCI([11:19]), 'k', 'Linewidth', 3); 


plot ([22:24], CCW([20:22]), 'b', 'Linewidth', 3); 
plot ([22:24], CCB([20:22]), 'r', 'Linewidth', 3); 
plot ([22:24], CCI([20:22]), 'k', 'Linewidth', 3); 

set(gca, 'FontSize', 28, 'xlim', [0 25], 'ylim', [-.2 1], 'xtick', [5 16 23],'xticklabels', [{'Enc', 'Dec', 'Class'}] )
plot([11 11], get(gca, 'ylim'), 'k:', LineWidth=2)
plot([21 21], get(gca, 'ylim'), 'k:', LineWidth=2)
%legend({'Within category', 'Between category', 'Within - between'},   'Location','northwest', 'FontSize',16); 
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);



%%

x = repelem(1:9, 270);
aWH = aW(:, 1:9)
y = aWH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 9.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 270);
aWH = aW(:, 10:18)
y = aWH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 9.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 1500);
aBH = aB(:, 1:9); 
y = aBH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 1500);
aBH = aB(:, 10:18); 
y = aBH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);







%% ENCODER

close all
clear bCMR wCMR


x0 = 1:9;
for i = 1:1500
    
    y0 = aB(i,1:9); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    bCMR(i,:) = b; 

    %y = b(1) + x0*b(2)
    %figure()
    %plot(x0,y0,'o')
    %hold on
    %plot(x0,y,'--r')
    %h1 = lsline
    
end

for i = 1:270

    y0 = aW(i,1:9); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(wCMR(:, 2));
t = ts.tstat;
disp(['t = ' num2str(t) ' // p = ' num2str(p)])


%% DECODER 
close all
clear bCMR wCMR


x0 = 1:9;
for i = 1:1500
    
    y0 = aB(i,10:18); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    bCMR(i,:) = b; 

    %y = b(1) + x0*b(2)
    %figure()
    %plot(x0,y0,'o')
    %hold on
    %plot(x0,y,'--r')
    %h1 = lsline
    
end

for i = 1:270

    y0 = aW(i,10:18); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(bCMR(:, 2));
t = ts.tstat;
disp(['t = ' num2str(t) ' // p = ' num2str(p)])




%% AUTOENCODER load network trained on ecoset
% % % 
clear, clc
%f2sav =  'AE-0000EN_hipp_E123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
f2sav =  'AE-1000EN_hipp_E123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
cfg = getParams(f2sav);
sessi = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
net2load = strsplit(cfg.net2load, '-'); 
net2load = net2load{2};

[ACT layNames] = load_AE_ECO_activ(cfg, sessi, paths, net2load);

%% Representational consistency all layers / time points 

nLays = size(ACT, 1); 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS_FR = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_FR(nLays+1,:) = nan; allMS_FR(:,nLays+1) = nan; %need this for the pColor below
%allMS_FR = tril(allMS_FR, -1); allMS_FR(allMS_FR==0) = nan; 

allMS = squeeze(mean(allMS_FR, 3))

% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 10);
set(gca, 'xtick',  (1:nLays) +.5, 'ytick', (1:nLays) +.5, 'xticklabels', layNames, 'yticklabels', layNames, 'TickLabelInterpreter','none', 'clim', [0 1])%
colormap(myCmap)


exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

nLays = size(ACT, 1); 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, nLays, []); act3 = act2(:,all(~isnan(act2)));  
allMS = corr(act3', 'type', 's'); %allMS = allM.^2;


figure(); set(gcf, 'Position', [100 100 800 800])
d2p = 1- allMS;
cols = zeros(nLays, 3); cols(1:10, 1) = 1; cols(11:19, 3) = 1; 
[rdmMDS] = cmdscale(d2p);
o2u = .03; 
scatter(rdmMDS(:,1),rdmMDS(:,2), 1550, cols, '.'); axis square
text(rdmMDS(:,1)+o2u,rdmMDS(:,2),layNames, 'Interpreter','none');
set(gca, 'xlim', [-.6 .6], 'ylim', [-.6 .6])

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% RNN all RDMS

nLays = size(ACT, 1); 
figure()
for layi = 1:nLays
   subplot (5, 5, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [-.5 .5])
   %set(gca,'cLim', [0 1])
   %colorbar
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   title(layNames(layi), 'Interpreter','none')
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = nLays+1;
for rowi = 1:5
    for coli = 1:5
        %set(ha(count),'position',[0+coli/6 0+rowi/6 n n ])
        count = count+-1;
        
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS

nLays = size(ACT, 1); 
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:nLays
   subplot (5, 5, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   title(layNames(layi), 'Interpreter','none')
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = nLays+1;
for rowi = 1:5
    for coli = 1:5
        %set(ha(count),'position',[0+coli/6 0+rowi/6 n n ])
        count = count-1;
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer 

clc
clearvars -except ACT

%create cateogy model
M = kron(eye(6), ones(10));
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCW(layi) = mWithin;
    CCB(layi) = mAcross;
    CCI(layi) = mWithin-mAcross ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_CH(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_CH(:, layi) = tmp;

end


figure(); set(gcf, 'Position', [100 100 800 800])
plot ([1:10], CCW([1:10]),'b', 'Linewidth', 3); hold on; %axis square
plot ([12:20], CCW([11:19]), 'b', 'Linewidth', 3); 
plot ([1:10], CCB([1:10]), 'r', 'Linewidth', 3); 
plot ([12:20], CCB([11:19]), 'r', 'Linewidth', 3); 
plot ([1:10], CCI([1:10]), 'k', 'Linewidth', 3); 
plot ([12:20], CCI([11:19]), 'k', 'Linewidth', 3); 
plot([11 11], get(gca, 'ylim'), 'k:', LineWidth=2)
set(gca, 'xlim', [0 21])

%legend({'Within category', 'Between category', 'Within - between'},   'Location','northwest', 'FontSize',16); 
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% compute CCI for each layer INCLUDING CLASSIFICATION HEAD

clc
clearvars -except ACT

%create cateogy model
M = kron(eye(6), ones(10));

for layi = 1:size(ACT, 1)
    rdmMDS = squeeze(ACT(layi, :,:)); 
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCW(layi) = mWithin;
    CCB(layi) = mAcross;
    CCI(layi) = mWithin-mAcross ;

    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_CH(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_CH(:, layi) = tmp;

end

figure(); set(gcf, 'Position', [100 100 800 800])
plot ([1:10], CCW([1:10]),'b', 'Linewidth', 3); hold on; %axis square
plot ([12:20], CCW([11:19]), 'b', 'Linewidth', 3); 
plot ([1:10], CCB([1:10]), 'r', 'Linewidth', 3); 
plot ([12:20], CCB([11:19]), 'r', 'Linewidth', 3); 
plot ([1:10], CCI([1:10]), 'k', 'Linewidth', 3); 
plot ([12:20], CCI([11:19]), 'k', 'Linewidth', 3); 


plot ([22:24], CCW([20:22]), 'b', 'Linewidth', 3); 
plot ([22:24], CCB([20:22]), 'r', 'Linewidth', 3); 
plot ([22:24], CCI([20:22]), 'k', 'Linewidth', 3); 

set(gca, 'FontSize', 28, 'xlim', [0 25], 'ylim', [-.2 1], 'xtick', [5 16 23],'xticklabels', [{'Enc', 'Dec', 'Class'}] )
plot([11 11], get(gca, 'ylim'), 'k:', LineWidth=2)
plot([21 21], get(gca, 'ylim'), 'k:', LineWidth=2)
%legend({'Within category', 'Between category', 'Within - between'},   'Location','northwest', 'FontSize',16); 
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);



%%

x = repelem(1:9, 270);
aWH = aW(:, 1:9)
y = aWH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 9.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 270);
aWH = aW(:, 10:18)
y = aWH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 9.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 1500);
aBH = aB(:, 1:9); 
y = aBH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:9, 1500);
aBH = aB(:, 10:18); 
y = aBH(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);







%% ENCODER

close all
clear bCMR wCMR


x0 = 1:9;
for i = 1:1500
    
    y0 = aB(i,1:9); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    bCMR(i,:) = b; 

    %y = b(1) + x0*b(2)
    %figure()
    %plot(x0,y0,'o')
    %hold on
    %plot(x0,y,'--r')
    %h1 = lsline
    
end

for i = 1:270

    y0 = aW(i,1:9); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(wCMR(:, 2));
t = ts.tstat;
disp(['t = ' num2str(t) ' // p = ' num2str(p)])


%% DECODER 
close all
clear bCMR wCMR


x0 = 1:9;
for i = 1:1500
    
    y0 = aB(i,10:18); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    bCMR(i,:) = b; 

    %y = b(1) + x0*b(2)
    %figure()
    %plot(x0,y0,'o')
    %hold on
    %plot(x0,y,'--r')
    %h1 = lsline
    
end

for i = 1:270

    y0 = aW(i,10:18); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(bCMR(:, 2));
t = ts.tstat;
disp(['t = ' num2str(t) ' // p = ' num2str(p)])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%FIT TO NEURAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc

%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
   
listF2sav = {

'AE0000N_hipp_E123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE1000N_hipp_E123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE0000N_vvs_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE1000N_vvs_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE0000N_pfc_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE1000N_pfc_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE0000N_hipp_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
'AE1000N_hipp_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 


};   

t1 = datetime; 
for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFilesWM(paths.powerFromRT);
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.powerFromRT filelistSess{sessi}]);   
        
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
        if strcmp(cfg.meth, 'LATERAL')
            cfg_contrasts = take_only_lateral_electrodes(cfg_contrasts);
        end
        if length(cfg_contrasts.oneListIds) > 1 & size(cfg_contrasts.chanNames, 1) > 1
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);

            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
                
            nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
            nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
        end
    end
    
    save([paths.results.DNNs f2sav '.mat'], 'nnFit');

end
t2 = datetime; 
etime(datevec(t2), datevec(t1))



%%  Load and plot results 
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav = 'AE1000N_pfc_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 



cLim = [-5 5]; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{3}, 1) % nnFit{1} is empty in PFC
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
         if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
           end
         end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    % if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
    %     h(:, 1:13) = 0; % only sum p-values in clusters after the baseline
    % end
    
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
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', cLim, 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', cLim);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    end

    if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', cLim, 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        plot([5+8 5+8],get(gca,'ylim'), 'k:','lineWidth', 2);
        % set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        % set(gca, 'xlim', [6 25], 'clim', [-5 5], 'FontSize', 10);
        
        
    end
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



















































%%

































%%
