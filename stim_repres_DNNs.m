%% STIMULUS REPRESENTATION IN DNNs
% % ALEXNET 
% % % 
clear, clc
f2sav = 'Alex_pfc_E123_[1-8]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

subj_ch_fr = 1; 
ACT_FR = load_alex_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
subj_ch_fr = 20; 
ACT_CH = load_alex_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


%% Plot all layers Alexnet one line horizontal

ACT = ACT_CH; 

figure(); set(gcf, 'Position', [100 100 1500 500]);
myCmap = brewermap([], '*Spectral') 
colormap(myCmap)

for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot MDS one Line Vertical

ACT = ACT_FR; 

cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);
%plot(1,1,'.','color',cols(11,:), 'Markersize', 2000) % check the color
% also nice palette here: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5a6dcde-4a80-11e4-9553-005056977bd0/a64ed616-e9ce-4b1d-8d18-5118cc03f8d2/images/screenshot.png

figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
 %legend
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Representational consistency all layers / time points Alexnet

ACT = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  
allMS_FR = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_FR(9,:) = nan; allMS_FR(:,9) = nan; %need this for the pColor below
allMS_FR = tril(allMS_FR, -1); allMS_FR(allMS_FR==0) = nan; 

ACT = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  
allMS_CH = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_CH(9,:) = nan; allMS_CH(:,9) = nan; %need this for the pColor below
allMS_CH = tril(allMS_CH, -1); allMS_CH(allMS_CH==0) = nan; 

allMS2 = cat(3, allMS_FR, allMS_CH);
allMS = squeeze(mean(allMS2, 3))

% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0.5 1], 'xtick',  (1:8) +.5, 'ytick', (1:8) +.5, ...
                'xticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'}, ...
                'yticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer Alexnet

clc
clearvars -except ACT_FR ACT_CH 

%create cateogy model
M = kron(eye(6), ones(10)); %equivalent to
%M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

% freiburg
ACT = ACT_FR;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_FR(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_FR(:, layi) = tmp;
    CCW_FR(layi) = mWithin;
    CCB_FR(layi) = mAcross;
    CCI_FR(layi) = mWithin-mAcross ;
end

%china
ACT = ACT_CH;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    tmp = rdmMDS(M == 1);
    tmp(any(isnan(tmp), 2), :) = [];
    aW_CH(:, layi) = tmp
    tmp = rdmMDS(M == 0);
    tmp(any(isnan(tmp), 2), :) = [];
    aB_CH(:, layi) = tmp;
    CCW_CH(layi) = mWithin;
    CCB_CH(layi) = mAcross;
    CCI_CH(layi) = mWithin-mAcross ;
    %CCI_CH(layi) = (mAcross - mWithin) ;
end


aW = squeeze(mean(cat(3, aW_CH, aW_FR), 3)); 
aB = squeeze(mean(cat(3, aB_CH, aB_FR), 3)); 


CCI2 = [CCI_FR; CCI_CH]; CCW2 = [CCW_FR; CCW_CH]; CCB2 = [CCB_FR; CCB_CH]; 
CCI = mean(CCI2); CCW = mean(CCW2); CCB = mean(CCB2);

figure(); set(gcf, 'Position', [100 100 500 400])
plot (CCW, 'Linewidth', 3); hold on; %axis square
plot (CCB, 'Linewidth', 3); %axis square
plot (CCI, 'Linewidth', 3); hold on; %axis square
set(gca, 'FontSize', 28, 'xlim', [0 9], 'ylim', [0 1.1])
%legend({'Within category', 'Between category', 'Within - between'})
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%
close all
clear bCMR wCMR


x0 = 1:8;
for i = 1:1500
    
    y0 = aB(i,:); 
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

    y0 = aW(i,:); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(bCMR(:, 2))
%t = ts.tstat 
%disp([num2str(p) ' // ' num2str(t)])

%% 
[h p ci ts] = ttest(wCMR(:, 2))
%t = ts.tstat; 
%disp([num2str(p) ' // ' num2str(t)])



%%

%%

x = repelem(1:8, 270);
y = aW(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 8.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:8, 1500);
y = aB(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 8.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%

x = repelem(1:8, 270);
y = aW(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 8.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);



%%  permutations
clearvars -except ACT CCI

nPerm = 1000; 

for permi = 1:nPerm
    %create cateogy model
    idSh = randperm(60*60);
    M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 
    M = M(idSh); M = reshape(M, [60 60]);
    
    for layi = 1:size(ACT, 1)
        d2p = squeeze(ACT(layi, :,:)); 
        d2p = 1- d2p;
        rdmMDS = d2p; 
        rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
        mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
        mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
        CCIP(permi, layi) = (mAcross - mWithin);
    
    end
        
end

%% 
lay = 8; 
obsT = CCI(lay); 
permT = CCIP(:, lay);


figure
histogram(permT); hold on; 
scatter(obsT,0, 'filled','r');

%% compute mean at each layer 


for layi = 1:8 

    mPT(layi) = mean(CCIP(:, layi));




end


%% % % % % % %% % % % % %% % % % % %% % % % % %% % % % % %% % % % % %% % % % % %% % RNN 
% % % 
clear, clc
f2sav = 'BLNETi_pfc_E123_[1-56]_3-54_0_0_0_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

subj_ch_fr = 1; 
[ACT_FR] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet

subj_ch_fr = 20; 
[ACT_CH] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet




%% RNN all RDMS

ACT = ACT_FR;

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   title(num2str(layi))
end

%exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% RNN all RDMS

ACT = ACT_FR; 

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   %title(num2str(layi)) % check order
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+coli/9 0+rowi/9 n n ])
        count = count+-1;
        
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS

ACT = ACT_FR; 

cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+rowi/9 0+coli/9 n n ])
        count = count-1;
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer (final timepoint only)

clc
clearvars -except ACT_FR ACT_CH 

%create cateogy model
M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

% freiburg
ACT = ACT_FR(8:8:56, :, :);
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
ACT = ACT_CH(8:8:56, :, :);;
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


figure(); set(gcf, 'Position', [100 100 500 400])
plot (CCW, 'Linewidth', 3); hold on; %axis square
plot (CCB, 'Linewidth', 3); %axis square
plot (CCI, 'Linewidth', 3); hold on; %axis square
set(gca, 'FontSize', 28, 'xlim', [0 8], 'ylim', [0 1])
legend({'Within category', 'Between category', 'Within - between'})
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%

x = repelem(1:7, 270);
y = aW(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

x = repelem(1:7, 1500);
y = aB(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);





%%
close all
clear bCMR wCMR


x0 = 1:7;
for i = 1:1500
    
    y0 = aB(i,:); 
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

    y0 = aW(i,:); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(bCMR(:, 2))

%% 
[h p ci ts] = ttest(wCMR(:, 2))



%% compute CCI for each layer ALL TIME POINTS

clc
clearvars -except ACT_FR ACT_CH 

%create cateogy model
M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

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


%%

x = repelem(1:7, 270);
%y = aW(:,1:8); y = (:);

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%%

%x = repelem(1:7, 1500);
%y = aB(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 7.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);




%% compute CCI for each layer / timepoint


% First create cateogy model
M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

%Freiburg data
ACT = ACT_FR;         
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    CCI_FR(layi) = (mAcross - mWithin) ;
end

%China data
ACT = ACT_CH;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    CCI_CH(layi) = (mAcross - mWithin) ;
end

CCI2 = [CCI_FR; CCI_CH]; 
CCI = mean(CCI2); 


clear c1 c2 c3 cols
c1 = (1:7)/7'; c2 = repmat(zeros(1), 7, 1)'; c3 = repmat(ones(1), 7, 1)'; cols = [c1 ; c2 ; c3]';

figure(); set(gcf, 'Position', [100 100 500 400])
lw = 4;
plot (CCI(1:8), 'Linewidth', lw, 'Color', cols(1, :)); hold on; 
plot (CCI(9:16), 'Linewidth', lw, 'Color', cols(2, :)); hold on; 
plot (CCI(17:24), 'Linewidth', lw, 'Color', cols(3, :)); hold on; 
plot (CCI(25:32), 'Linewidth', lw, 'Color', cols(4, :)); hold on; 
plot (CCI(33:40), 'Linewidth', lw, 'Color', cols(5, :)); hold on; 
plot (CCI(41:48), 'Linewidth', lw, 'Color', cols(6, :)); hold on; 
plot (CCI(49:56), 'Linewidth', lw, 'Color', cols(7, :)); hold on; 
set(gca, 'FontSize', 22, 'xlim', [0.5 8.5], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% compute CCI for each layer / timepoint FOR ECO AND IMAGE

clc
clearvars -except ACT_ECO ACT_IMA

%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(ACT_ECO, 1)
    d2p_ECO = squeeze(ACT_ECO(layi, :,:)); 
    d2p_ECO = 1- d2p_ECO;
    rdmMDS_ECO = d2p_ECO; 
    rdmMDS_ECO(find(eye(size(rdmMDS_ECO)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS_ECO(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS_ECO(M == 0), 'all', 'omitnan');
    CCI_ECO(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    %CCI_ECO(layi) = (mAcross - mWithin) ;
    
    d2p_IMA = squeeze(ACT_IMA(layi, :,:)); 
    d2p_IMA = 1- d2p_IMA;
    rdmMDS_IMA = d2p_IMA; 
    rdmMDS_IMA(find(eye(size(rdmMDS_IMA)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS_IMA(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS_IMA(M == 0), 'all', 'omitnan');
    CCI_IMA(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    %CCI_IMA(layi) = (mAcross - mWithin) ;
    
    
end

%% 

figure()
%plot (CCI_ECO, 'Linewidth', 2); hold on; 
%plot (CCI_IMA, 'Linewidth', 2); hold on; 
plot(CCI_ECO-CCI_IMA, 'Linewidth', 2); hold on; 
set(gca, 'FontSize', 20)




%%  permutations
clearvars -except act_CH act_FR

nPerm = 1000; 

for permi = 1:nPerm
%create cateogy model
idSh = randperm(60*60);
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 
M = M(idSh); M = reshape(M, [60 60]);

    for layi = 1:size(act_CH, 1)
        d2p = squeeze(act_FR(layi, :,:)); 
        d2p = 1- d2p;
        rdmMDS = d2p; 
        rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
        mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
        mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
        %CCIP(permi, layi) = (mAcross - mWithin) / (mAcross + mWithin);
        CCIP(permi, layi) = (mAcross - mWithin);

    end
    
end





%% plot only last time point one line

ACT = ACT_FR;
figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *6  0.1 n n ])
 set(ha(2),'position',[.09 * 5 0.1 n n ])
 set(ha(3),'position',[.09 * 4 0.1 n n])
 set(ha(4),'position',[.09 * 3 0.1 n n ])
 set(ha(5),'position',[.09 * 2 0.1 n n ])
 set(ha(6),'position',[.09 0.1 n n])
 set(ha(7),'position',[.0 0.1 n n ])

myCmap = brewermap([], '*Spectral') 
colormap (myCmap)

exportgraphics(gcf, 'allM.png', 'Resolution', 300);





%% plot MDS only last time point one line
cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 56];
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *6  0.1 n n ])
 set(ha(2),'position',[.09 * 5 0.1 n n ])
 set(ha(3),'position',[.09 * 4 0.1 n n])
 set(ha(4),'position',[.09 * 3 0.1 n n ])
 set(ha(5),'position',[.09 * 2 0.1 n n ])
 set(ha(6),'position',[.09 0.1 n n])
 set(ha(7),'position',[.0 0.1 n n ])
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% MDS all layers (black and yellow circles plot)

act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
c1 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c2 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c3 = repmat(zeros(1), 56, 1);

cols = [c1 c2 c3];



d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);
figure()
plot(rdmMDS(1:8,1),rdmMDS(1:8,2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS(9:16,1),rdmMDS(9:16,2),'k', 'linewidth', 2)
plot(rdmMDS(17:24,1),rdmMDS(17:24,2),'k', 'linewidth', 2)
plot(rdmMDS(25:32,1),rdmMDS(25:32,2),'k', 'linewidth', 2)
plot(rdmMDS(33:40,1),rdmMDS(33:40,2),'k', 'linewidth', 2)
plot(rdmMDS(41:48,1),rdmMDS(41:48,2),'k', 'linewidth', 2)
plot(rdmMDS(49:56,1),rdmMDS(49:56,2),'k', 'linewidth', 2)
scatter(rdmMDS(:,1),rdmMDS(:,2),3500,cols, '.'); 
set(gca,'FontSize', 26);



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

%exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% MDS all layers (pink and blue triangle-circles plot) only for first and last time points

act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  
allM_FR = corr(act3', 'type', 'p');%allMS = allM.^2;

act = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  
allM_CH = corr(act3', 'type', 'p');%allMS = allM.^2;

allM2 = cat(3, allM_FR, allM_CH);

allM = squeeze(mean(allM2, 3));
%allM = allM_CH;

clear c1 c2 c3 cols
c1 = (1:7)/7'; c2 = repmat(zeros(1), 7, 1)';c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);


allM = allM([1 8 9 16 17 24 25 32 33 40 41 48 49 56], [1 8 9 16 17 24 25 32 33 40 41 48 49 56]); 

d2p = 1- allM;

[rdmMDS] = cmdscale(d2p);
figure()
plot(rdmMDS([1 2],1),rdmMDS([1 2],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([3 4],1),rdmMDS([3 4],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([5 6],1),rdmMDS([5 6],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([7 8],1),rdmMDS([7 8],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([9 10],1),rdmMDS([9 10],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([11 12],1),rdmMDS([11 12],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([13 14],1),rdmMDS([13 14],2),'k', 'linewidth', 2);hold on; 



fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 7);


for i = 1:14
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gcf, 'Position', [100 100 500 400])
    set(gca,'FontSize', 22, 'xlim', [-.4 .4], 'ylim', [-.2 .4]); 
end



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% Representational consistency all layers / time points

act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 15)


%% Representational consistency only first and last time point

%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 54]; 


act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  
allM_FR = corr(act3', 'type', 'p');%allMS = allM.^2;

act = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  
allM_CH = corr(act3', 'type', 'p');%allMS = allM.^2;

allM2 = cat(3, allM_FR, allM_CH);

allM = squeeze(mean(allM2, 3));


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);
allMS(8,:) = nan; allMS(:,8) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 

% % % matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse');
set(gca, 'clim', [0.5 1], 'xtick',  (1:7) +.5, 'ytick', (1:7) +.5, 'xticklabels', {'1', '2', '3', '4', '5', '6', '7'}, ...
                'yticklabels', {'1', '2', '3', '4', '5', '6', '7'})

set(gca, 'FontSize', 25)
colormap(myCmap)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% MDS Representational consistency only first and last time point

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 


act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);

[rdmMDS] = cmdscale(allMS);

clear c1 c2 c3 cols
c1 = (1:7)/7';
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';

figure()
for i = 1:7
   scatter(rdmMDS(i,1),rdmMDS(i,2),5350, cols(i, :), '.'); hold on; axis square
   
end
set(gca,'FontSize', 40, 'xlim', [-.5 .5], 'ylim', [-.5 .5]); 

exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Changes in representational consistency 

act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  

c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';

clear allM2_FR
for tlyi = 1:55
    allM2_FR(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');
end

act = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  

c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';

clear allM2_CH
for tlyi = 1:55
    allM2_CH(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');
end

allM2 = [allM2_FR  allM2_CH];
allM2 = squeeze(mean(allM2, 2))


figure(); set(gcf, 'Position', [100 100 500 400])
lw = 3;

plot(abs(allM2(1:7)), 'Linewidth', lw, 'Color', [cols(1,:)]); hold on; 
plot(abs(allM2(9:15)), 'Linewidth', lw, 'Color', [cols(2,:)])
plot(abs(allM2(17:23)), 'Linewidth', lw, 'Color', [cols(3,:)])
plot(abs(allM2(25:31)), 'Linewidth', lw, 'Color', [cols(4,:)])
plot(abs(allM2(33:39)), 'Linewidth', lw, 'Color', [cols(5,:)])
plot(abs(allM2(41:47)), 'Linewidth', lw, 'Color', [cols(6,:)])
plot(abs(allM2(49:55)), 'Linewidth', lw, 'Color', [cols(7,:)])
%set(gca, 'ylim', [0.8 1], 'xlim', [0.5 7.5], 'xtick', [1:7], 'xticklabels', {'1-2' '2-3' '3-4' '4-5' '5-6' '6-7' '7-8'}, 'FontSize', 22)
set(gca, 'ylim', [0.8 1], 'xlim', [0.5 7.5], 'xtick', [1:7], 'xticklabels', {'1' '2' '3' '4' '5' '6' '7'}, 'FontSize', 22)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%%
for i = 1:54
   
    allps(i,:) = compare_correlation_coefficients(allM2(i), allM2(i+1), 1770, 1770)
    allp2s(i,:) = corr_rtest(allM2(i), allM2(i+1), 1770, 1770)

end


%%
d2p = allM2(1:7); 
d2p2 = stdfilt(d2p)
%diff(d2p)

plot(d2p); hold on; 
plot(d2p2)



%% 
% % CORNET 
% % % 
clear, clc
%cfg = getParams(f2sav);
cfg.lays2load = 1:20;
sessi = 1; 
paths = load_paths_WM('pfc', 'CORrt');

subj_ch_fr = 1; 
[ACT_FR] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
%ACT = ACT([1:18 20:22 24:26 29:30 33:34], :, :);
subj_ch_fr = 20; 
[ACT_CH] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet




%% corNET all RDMS

ACT = ACT_FR;

figure()
for layi = 1:20
   subplot (7, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   set(gca,'cLim', [0 1])
   title(num2str(layi))
end

%exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);



%% Plot all layers CORNET one line horizontal

figure(); set(gcf, 'Position', [100 100 1500 500]);
myCmap = brewermap([], '*Spectral') 
colormap(myCmap)

%l2l = [5 10 14 18 22 26 30 34];
l2l = [5 10 15 20];

for layi = 1:4
   subplot (1, 4, layi)
   d2p = squeeze(ACT(l2l(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% plot MDS one Line Vertical
cols = brewermap(6, 'Accent') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);
%plot(1,1,'.','color',cols(11,:), 'Markersize', 2000) % check the color
% also nice palette here: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5a6dcde-4a80-11e4-9553-005056977bd0/a64ed616-e9ce-4b1d-8d18-5118cc03f8d2/images/screenshot.png

ACT = ACT_FR;
l2l = [5 10 15 20];
figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:4
   subplot (1, 4, layi)
   d2p = squeeze(ACT(l2l(layi), :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
  
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Representational consistency all layers / time points CORNET 
clearvars -except ACT_CH ACT_FR

%l2l = [5 14 22 30]; % first time-point
l2l = [5 10 15 20]; % last time-point


ACT = ACT_CH; 
ACT = ACT(l2l,:,:);
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 4, []); act3 = act2(:,all(~isnan(act2)));  
allMS_CH = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_CH(5,:) = nan; allMS_CH(:,5) = nan; %need this for the pColor below
allMS_CH = tril(allMS_CH, -1); allMS_CH(allMS_CH==0) = nan; 



ACT = ACT_FR; 
ACT = ACT(l2l,:,:);
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 4, []); act3 = act2(:,all(~isnan(act2)));  
allMS_FR = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS_FR(5,:) = nan; allMS_FR(:,5) = nan; %need this for the pColor below
allMS_FR = tril(allMS_FR, -1); allMS_FR(allMS_FR==0) = nan; 



allMS2 = cat(3, allMS_FR, allMS_CH);
allMS = squeeze(mean(allMS2, 3))

% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0.5 1], 'xtick',  (1:4) +.5, 'ytick', (1:4) +.5, 'xticklabels', {'V1','V2', 'V4','IT'}, 'yticklabels', {'V1','V2', 'V4','IT'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);








%% MDS all layers (pink and blue triangle-circles plot) only for INPUT AND OUTPUT CORNET

clearvars -except ACT_FR ACT_CH


act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 20, []); act3 = act2(:,any(~isnan(act2)));  
allM_FR = corr(act3', 'type', 's');

act = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 20, []); act3 = act2(:,any(~isnan(act2)));  
allM_CH = corr(act3', 'type', 's');

allM2 = cat(3, allM_FR, allM_CH);

allM = squeeze(mean(allM2, 3));

%allM(isnan(allM(:, 1)), :) = [];
%allM = allM(:,all(~isnan(allM)));   % for nan - columns


l2l = [1 5 7 10 13 15 19 20 ]; 
allM = allM(l2l, l2l);


d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);

rdmMDS = rdmMDS*2; % for visualization

figure()
plot(rdmMDS([1 2],1),rdmMDS([1 2],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([3 4],1),rdmMDS([3 4],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([5 6],1),rdmMDS([5 6],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([7 8],1),rdmMDS([7 8],2),'k', 'linewidth', 2);hold on; 

fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 8);

clear c1 c2 c3 cols
c1 = (1:4)/4'; c2 = repmat(zeros(1), 4, 1)';c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';
cols = repelem(cols,2, 1)

for i = 1:8
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gcf, 'Position', [100 100 500 400])
    set(gca,'FontSize', 28, 'xlim', [-1 1], 'ylim', [-1 1]); 
end


exportgraphics(gcf, 'allM.png', 'Resolution', 300);








%% compute CCI for each layer and timepoint Cornet 

clearvars -except ACT_FR ACT_CH 

%create cateogy model
M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

% freiburg
ACT = ACT_FR;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI_FR(layi) = (mAcross - mWithin) ;
end

%china
ACT = ACT_CH;
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    CCI_CH(layi) = (mAcross - mWithin) ;
end

CCI2 = [CCI_FR; CCI_CH]; 
CCI = mean(CCI2); 

% % % mds
clear c1 c2 c3 cols
c1 = (1:4)/4'; % sorted by time point
c2 = repmat(zeros(1), 4, 1)';
c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';


figure(); set(gcf, 'Position', [100 100 500 400])
plot (1:5, CCI(1:5), 'Linewidth', 3, 'Color', [cols(1,:)]); hold on; %axis square
plot (2:5,CCI(7:10), 'Linewidth', 3, 'Color', [cols(2,:)]); %axis square
plot (3:5,CCI(13:15), 'Linewidth', 3, 'Color', [cols(3,:)]); %axis square
plot (4:5,CCI(19:20), 'Linewidth', 3, 'Color', [cols(4,:)]); %axis square
set(gca, 'FontSize', 22, 'xlim', [0 6], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);







%% compute CCI for each layer FINAL TIME POINT

clearvars -except ACT_FR ACT_CH 

l2l = [5 10 15 20];

%create cateogy model
M = zeros (60); M(1:10, 1:10) = 1; M(11:20, 11:20) = 1; M(21:30, 21:30) = 1; M(31:40, 31:40) = 1; M(41:50, 41:50) = 1; M(51:60, 51:60) = 1; 

% freiburg
ACT = ACT_FR;
for layi = 1:size(l2l, 2)
    d2p = squeeze(ACT(l2l(layi), :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    rdmMDS = tril(rdmMDS, -1); rdmMDS(rdmMDS==0) = nan; 
    
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
for layi = 1:size(l2l, 2)
    d2p = squeeze(ACT(l2l(layi), :,:)); 
    %d2p = 1- d2p;
    rdmMDS = d2p; 
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
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

figure(); set(gcf, 'Position', [100 100 500 400])
plot (CCW, 'Linewidth', 3); hold on; %axis square
plot (CCB, 'Linewidth', 3); %axis square
plot (CCI, 'Linewidth', 3); hold on; %axis square
set(gca, 'FontSize', 28, 'xlim', [0 5], 'ylim', [0 1.1])
%legend({'Within category', 'Between category', 'Within - between'})
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);














%%  permutations
clearvars -except ACT CCI

nPerm = 1000; 

for permi = 1:nPerm
    %create cateogy model
    idSh = randperm(60*60);
    M = zeros (60); 
    M(1:10, 1:10) = 1; 
    M(11:20, 11:20) = 1; 
    M(21:30, 21:30) = 1; 
    M(31:40, 31:40) = 1; 
    M(41:50, 41:50) = 1; 
    M(51:60, 51:60) = 1; 
    M = M(idSh); M = reshape(M, [60 60]);
    
    for layi = 1:size(ACT, 1)
        d2p = squeeze(ACT(layi, :,:)); 
        d2p = 1- d2p;
        rdmMDS = d2p; 
        rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
        mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
        mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
        %CCIP(permi, layi) = (mAcross - mWithin) / (mAcross + mWithin);
        CCIP(permi, layi) = (mAcross - mWithin);
    
    end
        
end

%% 
lay = 1; 
obsT = CCI(lay); 
permT = CCIP(:, lay);


figure
histogram(permT); hold on; 
scatter(obsT,0, 'filled','r');

%% Changes in representational consistency CORNET

clearvars -except ACT_FR ACT_CH

act = ACT_FR; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 20, []); act3 = act2(:,any(~isnan(act2)));  
m = squeeze(act2(20,:,:));
act2_1 = reshape(m, 1, []); 

clear allM2_FR
for tlyi = 1:19
    allM2_FR(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');
end

act = ACT_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 20, []); act3 = act2(:,any(~isnan(act2)));  
clear allM2_CH
for tlyi = 1:19
    allM2_CH(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');
end

allM2 = [allM2_FR  allM2_CH];
allM2 = squeeze(mean(allM2, 2))

c1 = (1:4)/4'; c2 = repmat(zeros(1), 4, 1)';c3 = repmat(ones(1), 4, 1)';cols = [c1 ; c2 ; c3]';

figure(); set(gcf, 'Position', [100 100 500 400])
lw = 3;
plot(1:4, abs(allM2(1:4)), 'Linewidth', lw, 'Color', [cols(1,:)]); hold on; 
plot(2:4, abs(allM2(7:9)), 'Linewidth', lw, 'Color', [cols(2,:)])
plot(3:4, abs(allM2(13:14)), 'Linewidth', lw, 'Color', [cols(3,:)])
scatter (4, allM2(19), 1000, [cols(4,:)], '.')


set(gca, 'ylim', [0.75 1], 'xlim', [0.5 4.5], 'xtick', [1:4], 'xticklabels', ...
            {'1' '2' '3' '4'}, 'FontSize', 28)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%%
close all
clear bCMR wCMR


x0 = 1:4;
for i = 1:1500
    
    y0 = aB(i,:); 
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

    y0 = aW(i,:); 
    X1 = [ones(length(x0),1)  x0'];
    b = X1\y0';
    wCMR(i,:) = b; 

end

%%
[h p ci ts] = ttest(bCMR(:, 2))

%% 
[h p ci ts] = ttest(wCMR(:, 2))




%%

x = repelem(1:4, 270);
y = aW(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 4.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%%

x = repelem(1:4, 1500);
y = aB(:); 

figure()
scatter(x, y, '+')
l1 = lsline
l1.LineWidth = 3; l1.Color = 'black';
set(gca, 'xlim', [0.5 4.5], 'FontSize', 22)
[rho, pval] = corr(x', y, 'type', 's')
ylabel('Rho'); xlabel('Layer'); title(['Rho = ' num2str(rho) '  p = ' num2str(pval)])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);






%% correlate imagenet and ecoset activations 

for layi = 1:56

    rdmECO = squeeze(ACT_ECO(layi, :, :)); 
    rdmECO = vectorizeRDM(rdmECO);
    rdmIMA = squeeze(ACT_IMA(layi, :, :)); 
    rdmIMA = vectorizeRDM(rdmIMA);
    allRs(layi, :) = corr(rdmECO', rdmIMA', 'type', 's');

end 

%% 
figure()
plot(allRs, 'LineWidth', 3)
xlabel('Layer/Time')
ylabel('Rho')
set(gca, 'FontSize', 14)

%%

figure(); 
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';
lw = 3;
plot (allRs(1:8), 'Linewidth', lw, 'Color', cols(1, :)); hold on; 
plot (allRs(9:16), 'Linewidth', lw, 'Color', cols(2, :)); hold on; 
plot (allRs(17:24), 'Linewidth', lw, 'Color', cols(3, :)); hold on; 
plot (allRs(25:32), 'Linewidth', lw, 'Color', cols(4, :)); hold on; 
plot (allRs(33:40), 'Linewidth', lw, 'Color', cols(5, :)); hold on; 
plot (allRs(41:48), 'Linewidth', lw, 'Color', cols(6, :)); hold on; 
plot (allRs(49:56), 'Linewidth', lw, 'Color', cols(7, :)); hold on; 
set(gca, 'FontSize', 16, 'xlim', [1 8])
xlabel('Time')
ylabel('Rho')


%% load RDMS 2

clear 
load rdms

for layi = 1:8

    for layj = 1:56

        x = vectorizeRDM(squeeze(ACT_CH_Alex(layi, :, :))); 
        y = vectorizeRDM(squeeze(ACT_CH_BLNET(layj, :, :))); 
        allR(layi, layj) = corr(x, y, 'type', 's');
        
    end

end

%% 
clear 

lay2load = 8:8:56;
load rdms

for layi = 1:8

    for layj = 1:length(lay2load)

        x = vectorizeRDM(squeeze(ACT_CH_Alex(layi, :, :))); 
        y = vectorizeRDM(squeeze(ACT_CH_BLNET(lay2load(layj), :, :))); 
        [allR(layi, layj) allP(layi, layj)]= corr(x, y, 'type', 's');
        
    end

end



%% 
figure
imagesc(allR); axis square; colorbar



%% LOAD WITHIN AND BETWEEN CORRELATIONS IN ALL LAYERS 3 MODELS
clear 
load all_wb_corr


%% COMPARE MEAN SLOPES 
x = mean(alexnet_wCMR(:, 2));fprintf(1, '\n\n\n');
disp (['Alexnet within-category:   ' num2str(x)]);
x = mean(blnet_wCMR(:, 2));
disp (['BLNET within-category:   ' num2str(x)]);
x = mean(cornet_wCMR(:, 2));
disp (['cornet within-category:   ' num2str(x)]);
x = mean(alexnet_bCMR(:, 2));
disp (['Alexnet between-category:   ' num2str(x)]);
x = mean(blnet_bCMR(:, 2));
disp (['BLNET between-category:   ' num2str(x)]);
x = mean(cornet_bCMR(:, 2));
disp (['cornet between-category:   ' num2str(x)]);

%% COMPARE MEAN SLOPES 
x = mean(alexnet_wCMR(:, 2));fprintf(1, '\n\n\n');
disp (['Alexnet between-category:   ' num2str(x)]);
y = mean(blnet_bCMR(:, 2));
disp (['BLNET between-category:   ' num2str(y)]);
x-y


%% COMPARE SLOPES WITHIN EACH MODE
[h p ci ts] = ttest(alexnet_wCMR(:, 2)); fprintf(1, '\n\n\n');
disp (['Alexnet within-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(blnet_wCMR(:, 2)); disp(''); 
disp (['BLNET within-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(cornet_wCMR(:, 2)); disp(''); 
disp (['corNET-RT within-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

fprintf(1, '\n\n\n');
[h p ci ts] = ttest(alexnet_bCMR(:, 2)); fprintf(1, '\n\n\n');
disp (['Alexnet between-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(blnet_bCMR(:, 2)); disp(''); 
disp (['BLNET between-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(cornet_bCMR(:, 2)); disp(''); 
disp (['corNET-RT between-category:   t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% COMPARE CORRELATIONS ACROSS THREE MODELS BETWEEN
clear 
load all_wb_corr

fprintf(1, '\n\n\n');
[h p ci ts] = ttest(alexnet_bCMR(:, 2), blnet_bCMR(:, 2));
disp (['BETWEEN (ALEXNET vs BLNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);
[h p ci ts] = ttest(alexnet_bCMR(:, 2), cornet_bCMR(:, 2));
disp (['BETWEEN (ALEXNET vs CORNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);
[h p ci ts] = ttest(blnet_bCMR(:, 2), cornet_bCMR(:, 2));
disp (['BETWEEN (BLNET vs CORNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% COMPARE CORRELATIONS ACROSS THREE MODELS WITHIN
clear , clc
load all_wb_corr

fprintf(1, '\n\n\n');
[h p ci ts] = ttest(alexnet_wCMR(:, 2), blnet_wCMR(:, 2));fprintf(1, '\n');
disp (['WITHIN (ALEXNET vs BLNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);
[h p ci ts] = ttest(alexnet_wCMR(:, 2), cornet_wCMR(:, 2));fprintf(1, '\n');
disp (['WITHIN (ALEXNET vs CORNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);
[h p ci ts] = ttest(blnet_wCMR(:, 2), cornet_wCMR(:, 2));
disp (['WITHIN (BLNET vs CORNET):  t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% plot 2 bars
clear data
%data.data = [alexnet_wCMR(:, 2)  blnet_wCMR(:, 2)]; 
%data.data = [alexnet_wCMR(:, 2)  cornet_wCMR(:, 2)]; 
%data.data = [blnet_wCMR(:, 2)  cornet_wCMR(:, 2)]; 
%data.data = [alexnet_bCMR(:, 2)  blnet_bCMR(:, 2)]; 
%data.data = [alexnet_bCMR(:, 2)  cornet_bCMR(:, 2)]; 
data.data = [blnet_bCMR(:, 2)  cornet_bCMR(:, 2)]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2], data.data, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);





