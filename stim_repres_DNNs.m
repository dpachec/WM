%% STIMULUS REPRESENTATION IN DNNs
% % ALEXNET 
% % % 
clear, clc
f2sav = 'Alex_pfc_E123_[1-8]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
[ACT] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


%% Plot all layers Alexnet one line horizontal

figure(); set(gcf, 'Position', [100 100 1500 500]);
myCmap = brewermap([], 'Spectral') 
colormap(myCmap)

for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(ACT(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [-.2 .6])
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


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  


allMS = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS(9,:) = nan; allMS(:,9) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 


% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0 .9], 'xtick',  (1:8) +.5, 'ytick', (1:8) +.5, ...
                'xticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'}, ...
                'yticklabels', {'Conv1', 'Conv2', 'Conv3', 'Conv4', 'Conv5', 'Fc6', 'Fc7', 'Fc8'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for each layer Alexnet

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(act_CH, 1)
    d2p = squeeze(act_CH(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


figure()
plot (CCI, 'Linewidth', 3); %axis square
set(gca, 'FontSize', 25, 'xlim', [0 9], 'ylim', [0 .75])
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
lay = 8; 
obsT = CCI(lay); 
permT = CCIP(:, lay);


figure
histogram(permT); hold on; 
scatter(obsT,0, 'filled','r');



%% % RNN 
% % % 
clear, clc
f2sav = 'RNN_pfc_E123_[1-56]_3-8_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav);
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM(cfg.brainROI);
[ACT] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


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

%% RNN all RDMS

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

%% compute CCI for each layer / timepoint

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(act_CH, 1)
    d2p = squeeze(act_CH(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


figure()
plot (CCI, 'Linewidth', 2)
set(gca, 'FontSize', 20)

figure
lw = 4;
plot (CCI(1:8), 'Linewidth', lw, 'Color', cols(1, :)); hold on; 
plot (CCI(9:16), 'Linewidth', lw, 'Color', cols(2, :)); hold on; 
plot (CCI(17:24), 'Linewidth', lw, 'Color', cols(3, :)); hold on; 
plot (CCI(25:32), 'Linewidth', lw, 'Color', cols(4, :)); hold on; 
plot (CCI(33:40), 'Linewidth', lw, 'Color', cols(5, :)); hold on; 
plot (CCI(41:48), 'Linewidth', lw, 'Color', cols(6, :)); hold on; 
plot (CCI(49:56), 'Linewidth', lw, 'Color', cols(7, :)); hold on; 
set(gca, 'FontSize', 25, 'xlim', [0 9], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% compute CCI for the last time point in each layer 

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(act_CH, 1)
    d2p = squeeze(act_CH(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


%figure()
%plot (CCI, 'Linewidth', 2)
%set(gca, 'FontSize', 20)

figure
lw = 4;
plot (CCI([8:8:56]), 'Linewidth', lw); hold on; 
set(gca, 'FontSize', 25, 'xlim', [.5 7.5], 'ylim', [0 .5])
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

figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 .8])
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




%% plot only last time point one line vertical

figure(); set(gcf, 'Position', [100 100 700 700]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(ACT(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .1; 
 set(ha(1),'position',[0 6/9 n n ])
 set(ha(2),'position',[0 5/9 n n ])
 set(ha(3),'position',[0 4/9 n n])
 set(ha(4),'position',[0 3/9 n n ])
 set(ha(5),'position',[0 2/9 n n ])
 set(ha(6),'position',[0 1/9 n n])
 set(ha(7),'position',[0 .0 n n ])
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

act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');%allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
%c2 = repmat(zeros(1), 7, 1)';

c2 = repmat(zeros(1), 7, 1)';

c3 = repmat(ones(1), 7, 1)';

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
    set(gca,'FontSize', 30, 'xlim', [-.4 .4], 'ylim', [-.2 .4]); 
end



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% Representational consistency all layers / time points

act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 15)


%% Representational consistency only first and last time point

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 


act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 

% % % matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse');
set(gca, 'FontSize', 15, 'clim', [0 1], 'xtick',  (1:7) +.5, 'ytick', (1:7) +.5, 'xticklabels', {'1', '2', '3', '4', '5', '6', '7', '8'}, ...
                'yticklabels', {'1', '2', '3', '4', '5', '6', '7', '8'})

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

act = ACT; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  

c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


clear allM2
for tlyi = 1:55

    allM2(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');allMS = allM;
    
end

figure()
lw = 3;

plot(abs(allM2(1:7)), 'Linewidth', lw, 'Color', [cols(1,:)]); hold on; 
plot(abs(allM2(9:15)), 'Linewidth', lw, 'Color', [cols(2,:)])
plot(abs(allM2(17:23)), 'Linewidth', lw, 'Color', [cols(3,:)])
plot(abs(allM2(25:31)), 'Linewidth', lw, 'Color', [cols(4,:)])
plot(abs(allM2(33:39)), 'Linewidth', lw, 'Color', [cols(5,:)])
plot(abs(allM2(41:47)), 'Linewidth', lw, 'Color', [cols(6,:)])
plot(abs(allM2(49:55)), 'Linewidth', lw, 'Color', [cols(7,:)])
%legend({'Layer 1' 'Layer 2' 'Layer 3' 'Layer 4' 'Layer 5' 'Layer 6' 'Layer 7' })
set(gca, 'ylim', [0.8 1], 'xlim', [0 8], 'xtick', [1:7], 'xticklabels', {'1-2' '2-3' '3-4' '4-5' '5-6' '6-7' '7-8'}, 'FontSize', 25)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% 
% % CORNET 
% % % 
clear, clc
%cfg = getParams(f2sav);
cfg.lays2load = 1:34;
sessi = 1; 
subj_ch_fr = 1; 
paths = load_paths_WM('pfc', 'EXAMPLE');
%[ACT] = load_CORr_TV(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
%ACT = ACT([1:18 20:22 24:26 29:30 33:34], :, :);

%% RNN all RDMS


figure()
for layi = 1:34
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
l2l = [10 18 26 34];

for layi = 1:4
   subplot (1, 4, layi)
   d2p = squeeze(ACT(l2l(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 .8])
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
clear, clc
cfg.lays2load = 1:34;
sessi = 1; subj_ch_fr = 1; 
paths = load_paths_WM('pfc', 'EXAMPLE');
[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet

l2l = [5 14 22 30]; % first time-point
%l2l = [10 18 26 34]; % last time-point


ACT = ACT(l2l,:,:);
act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 4, []); act3 = act2(:,all(~isnan(act2)));  


allMS = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS(5,:) = nan; allMS(:,5) = nan; %need this for the pColor below
allMS = tril(allMS, -1); allMS(allMS==0) = nan; 


% % % plot matrix
figure()
myCmap = brewermap([], '*Spectral') 
%imagesc(allMS); axis square; colorbar
s = pcolor(allMS); axis square; colorbar
set(s, 'edgecolor', 'none'); 
set(gca, 'ydir', 'reverse', 'FontSize', 20);
set(gca, 'clim', [0 .9], 'xtick',  (1:4) +.5, 'ytick', (1:4) +.5, ...
                'xticklabels', {'V1','V2', 'V4','IT'}, ...
                'yticklabels', {'V1','V2', 'V4','IT'})

colormap(myCmap)

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);






%% compute CCI for each layer Cornet

clc
clearvars -except act_CH act_FR ACT
act_CH = ACT;
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(act_CH, 1)
    d2p = squeeze(act_CH(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


figure()
plot (CCI, 'Linewidth', 3); %axis square
set(gca, 'FontSize', 25, 'xlim', [0 5], 'ylim', [0 .5], 'xtick', [1 2 3 4], 'xticklabels', {'V1' 'V2' 'V4' 'IT'})
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



%% MDS all layers (pink and blue triangle-circles plot) only for INPUT AND OUTPUT CORNET

clear, clc
cfg.lays2load = 1:34;
sessi = 1; subj_ch_fr = 1; 
paths = load_paths_WM('pfc', 'EXAMPLE');
[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
l2l = [5 10 14 18 22 26 30 34];

act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  
allM = corr(act3', 'type', 's'); %allMS = allM.^2;

clear c1 c2 c3 cols
c1 = (1:4)/4'; % sorted by time point
c2 = repmat(zeros(1), 4, 1)';
c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);

d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);

figure()
plot(rdmMDS([1 2],1),rdmMDS([1 2],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([3 4],1),rdmMDS([3 4],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([5 6],1),rdmMDS([5 6],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([7 8],1),rdmMDS([7 8],2),'k', 'linewidth', 2);hold on; 

fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 8);

for i = 1:8
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gca,'FontSize', 24, 'xlim', [-1 1], 'ylim', [-1 1]); 
end


exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% MDS all layers (pink and blue triangle-circles plot) only for OUTPUT
clear, clc
cfg.lays2load = 1:34;
sessi = 1; subj_ch_fr = 1; 
paths = load_paths_WM('pfc', 'EXAMPLE');
[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
l2l = [5 10 14 18 22 26 30 34];


act1 = arrayfun(@(i)tril(squeeze(ACT(i,:,:)), -1), 1:size(ACT,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  
allM = corr(act3', 'type', 's'); %allMS = allM.^2;
allMS = tril(allM, -1); allM(allM==0) = nan; 

clear c1 c2 c3 cols
c1 = (1:4)/4'; % sorted by time point
c2 = repmat(zeros(1), 4, 1)';
c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);

d2p = 1- allMS;
[rdmMDS] = cmdscale(allM);

figure()
fmt = {'o'; 'o'};
fmt = repelem(fmt, 1, 8);

for i = [2 4 6 8]
    scatter(rdmMDS(i,1),rdmMDS(i,2),300,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); hold on; 
    set(gca,'FontSize', 24, 'xlim', [-.35 .55], 'ylim', [-.25 .5]); 
end
box on; 


%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);





%% compute CCI for each layer Cornet WITH DIFFERENT TIMES

clear, clc
cfg.lays2load = 1:34;
sessi = 1; subj_ch_fr = 1; 
paths = load_paths_WM('pfc', 'EXAMPLE');
[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet
l2l = [5 10 14 18 22 26 30 34];
%create cateogy model
M = zeros (60); 
M(1:10, 1:10) = 1; 
M(11:20, 11:20) = 1; 
M(21:30, 21:30) = 1; 
M(31:40, 31:40) = 1; 
M(41:50, 41:50) = 1; 
M(51:60, 51:60) = 1; 

        
for layi = 1:size(ACT, 1)
    d2p = squeeze(ACT(layi, :,:)); 
    d2p = 1- d2p;
    rdmMDS = d2p; 
    %[rdmMDS] = cmdscale(d2p);
    %nDims = size(rdmMDS,2);
    rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
    mWithin = mean(rdmMDS(M == 1), 'all', 'omitnan');
    mAcross = mean(rdmMDS(M == 0), 'all', 'omitnan');
    %CCI(layi) = (mAcross - mWithin) / (mAcross + mWithin);
    CCI(layi) = (mAcross - mWithin) ;
    
end


% % % mds
clear c1 c2 c3 cols
c1 = (1:4)/4'; % sorted by time point
c2 = repmat(zeros(1), 4, 1)';
c3 = repmat(ones(1), 4, 1)';
cols = [c1 ; c2 ; c3]';


figure()
plot (1:5, CCI(6:10), 'Linewidth', 3, 'Color', [cols(1,:)]); hold on; %axis square
plot (1:4,CCI(15:18), 'Linewidth', 3, 'Color', [cols(2,:)]); %axis square
plot (1:3,CCI(24:26), 'Linewidth', 3, 'Color', [cols(3,:)]); %axis square
plot (1:2,CCI(33:34), 'Linewidth', 3, 'Color', [cols(4,:)]); %axis square
set(gca, 'FontSize', 25, 'xlim', [0 6], 'ylim', [0 .5])
exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);














