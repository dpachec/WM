%%
%%

clear, clc

region = 'pfc';
paths = load_paths_WM(region);

contrasts = {
              'DISC_EE' 
              'DIDC_EE';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEE c{i} '.mat']);
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



%% plot one example trial 

%d2p = squeeze(DISC_EE{1}(1, :, :)); 
d2p = squeeze(mean(DISC_EE{1}, 1, 'omitnan')); 
d2p =  triu(d2p.',1) + tril(d2p); 
imagesc(flipud(d2p)); axis square

%10ms
%set(gca, 'xtick', [0 51 150] , 'xticklabels', {'-.5' '0' '1.5'})
%set(gca, 'ytick', [1 102 151] , 'yticklabels', {'1.5' '0' '-.5'})

set(gca, 'xtick', [0 5 15] , 'xticklabels', {'-.5' '0' '1.5'})
set(gca, 'ytick', [6 16 22] , 'yticklabels', {'1.5' '0' '-.5'})
set(gca, 'FontSize', 22)




%% Compute dynamicity scores at the group level 
clearvars -except DISC_EE DIDC_EE out_c region



timeL =  [1:21]; 


meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
%exclude subjects
if strcmp(region, 'vvs')
    sub2exc = [18 22];
elseif strcmp(region, 'pfc')
    sub2exc = [1];
elseif strcmp(region, 'hipp')
    sub2exc = [2]
end
if sub2exc > 0
    meanC1(sub2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    meanC2(sub2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end
mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
nSubj = length(meanC1);
lM = size(meanC1{1},1 ); 
mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [lM lM nSubj]); mC1FD = permute(mC1FD, [3 1 2]);
mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [lM lM nSubj]); mC2FD = permute(mC2FD, [3 1 2]);

mDiff = mC1FD - mC2FD; 
mDiff = mDiff(:,timeL, timeL);




for timei = 1:size(mDiff, 2)

    for timej = 1:size(mDiff, 2)

        if timei~=timej
            offDiag = mDiff(:, timei, timej); 
            onDiag1 = mDiff(:, timei, timei); 
            onDiag2 = mDiff(:, timej, timej); 
            [h1 p1 ci1 ts1] = ttest(onDiag1, offDiag); 
            [h2 p2 ci2 ts2] = ttest(onDiag2, offDiag); 

             % % % % just 2 check the direction of the effect
%             if ts1.tstat > 0
%                 break
%             end


            if h1 & h2 & ts1.tstat > 0 & ts2.tstat > 0
                dynMatrx(timei, timej) = 1; 
            else
                dynMatrx(timei, timej) = 0; 
            end
        end
    end
end


clustinfo = bwconncomp(dynMatrx);
[clusL id] = max(cellfun(@length, clustinfo.PixelIdxList))







%% 
figure; 
imagesc(flipud(dynMatrx)); axis square
%set(gca, 'xtick', [0 51 150] , 'xticklabels', {'-.5' '0' '1.5'})
%set(gca, 'ytick', [1 102 151] , 'yticklabels', {'1.5' '0' '-.5'})

%set(gca, 'xtick', [0 5 15] , 'xticklabels', {'-.5' '0' '1.5'})
%set(gca, 'ytick', [])

set(gca, 'FontSize', 22)


%% plot each subject Diff matrix
figure; 
d2p = squeeze(mDiff(9, :, :));
imagesc(flipud(d2p)); axis square
%set(gca, 'xtick', [0 51 150] , 'xticklabels', {'-.5' '0' '1.5'})
%set(gca, 'ytick', [1 102 151] , 'yticklabels', {'1.5' '0' '-.5'})

%set(gca, 'xtick', [0 5 15] , 'xticklabels', {'-.5' '0' '1.5'})
%set(gca, 'ytick', [])

set(gca, 'FontSize', 22)






%% 
d2p = mean(dynMatrx, 1); 

figure
plot(d2p)





%% permutations 

nPerm = 100;

clear dynMatrxPerm mDiffPerm

meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 

%exclude subjects
if strcmp(region, 'vvs')
    sub2exc = [18 22];
elseif strcmp(region, 'pfc')
    sub2exc = [1];
elseif strcmp(region, 'hipp')
    sub2exc = [2];
end
if sub2exc > 0
    meanC1(sub2exc) = []; %all_cond1 = all_cond1(~cellfun('isempty',all_cond1));
    meanC2(sub2exc) = []; %all_cond2 = all_cond2(~cellfun('isempty',all_cond2));
end

mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
nSubj = length(meanC1);
lM = size(meanC1{1},1 ); 
mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [lM lM nSubj]); mC1FD = permute(mC1FD, [3 1 2]);
mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [lM lM nSubj]); mC2FD = permute(mC2FD, [3 1 2]);

mDiff = mC1FD - mC2FD; 
mDiff = mDiff(:, timeL,timeL );%take limits from previous analysis

for permi = 1:nPerm
    
    % % % shuffle matrix for all subjects 
    for subji = 1:nSubj
        mDiffSubj = squeeze(mDiff(subji, :, :)); 
        onDiag = diag(mDiffSubj); 
        for timei = 1:size(mDiff, 2)
            idR = randi(size(mDiff, 2)); 
            tmp = mDiffSubj(timei, idR); 
            onDiagPerm(timei) = tmp; 
            mDiffSubj(timei, idR) = onDiag(timei); 
            mDiffSubj(timei, timei) = onDiagPerm(timei); 
        end
        mDiffPerm(subji, :, :) = mDiffSubj ;
    end
    
    for timei = 1:size(mDiff, 2)
    
        for timej = 1:size(mDiff, 2)
    
            if timei~=timej
    
                
                offDiag = mDiffPerm(:, timei, timej); 
                onDiag1 = mDiffPerm(:, timei, timei); 
                onDiag2 = mDiffPerm(:, timej, timej); 
                [h1 p1 ci1 ts1] = ttest(onDiag1, offDiag); 
                [h2 p2 ci2 ts2] = ttest(onDiag2, offDiag); 
               
                
                if h1 & h2 & ts1.tstat > 0 & ts2.tstat > 0
                    %disp('one')
                    dynMatrxPerm(permi, timei, timej) = 1; 
                else
                    %disp('zero')
                    dynMatrxPerm(permi, timei, timej) = 0; 
                end
            end
        end
    end

    clustinfo = bwconncomp(squeeze(dynMatrxPerm(permi, :, :)));
    if ~isempty(clustinfo.PixelIdxList)
        clusLPerm(permi,:) = max(cellfun(@length, clustinfo.PixelIdxList));
    end

end



%% 
d2p = squeeze(dynMatrxPerm(1,:,:));
figure
contourf(d2p); axis square
clustinfo = bwconncomp(d2p)


%% 

figure
histogram(clusLPerm); hold on
scatter(clusL, 0, 2000, 'r.')


max(clusLPerm)







%% Compute dynamicity scores at the trial level 

% % % % % load data

clear, clc

region = 'vvs';
paths = load_paths_WM(region);

contrasts = {
              'DISC_EE' 
              'DIDC_EE';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSessEE c{i} '.mat']);
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


clearvars -except DISC_EE DIDC_EE out_c region




%% compute dynamicity scores

timeL = [1:151]; 

C1 = DISC_EE; 
C2 = DIDC_EE; 

% exclude subjects
if strcmp(region, 'vvs')
    sub2exc = [18 22];
elseif strcmp(region, 'pfc')
    sub2exc = [1];
elseif strcmp(region, 'hipp')
    sub2exc = [2]
end
if sub2exc > 0
    C1(sub2exc) = []; 
    C2(sub2exc) = []; 
end


for subji = 1:length(C1)

    clear clustL dyNB dynMTrial
    C1Sub = C1{subji}; 
    C2Sub = C2{subji}; 



    % % % compute mean C2 and substract same mean matrix to all trials ? ? ? 
     meanC2 = mean(C2Sub); 
     nTrials = size(C1Sub, 1); 
     meanC2 = repmat(meanC2, nTrials, 1, 1);
     %C1Sub = C1Sub - meanC2; 

    for triali = 1:size(C1Sub, 1)

        C1ST = squeeze(C1Sub(triali, :, :));
        C1ST =  triu(C1ST.',1) + tril(C1ST ); 
        C1ST = C1ST(timeL, timeL);
        dynMatrx = zeros(length(timeL));
           
        for timei = 1:size(C1ST, 2)
        
            for timej = 1:size(C1ST, 2)
        
                if timei~=timej
                    offDiag = C1ST(timei, timej); 
                    onDiag1 = C1ST(timei, timei); 
                    onDiag2 = C1ST(timej, timej); 
                    h1 = onDiag1 > offDiag;
                    h2 = onDiag2 > offDiag; 
                    
                    if h1 & h2 
                        dynMatrx(timei, timej) = 1; 
                    else
                        dynMatrx(timei, timej) = 0; 
                    end
                end
            end
        end
        
        
        clustinfo = bwconncomp(dynMatrx);
        biggestClust = max(cellfun(@length, clustinfo.PixelIdxList)); 
        if ~isempty(biggestClust)
            clustL(triali,:) = biggestClust; 
        else
            clustL(triali,:) = 0; 
        end
    
        dynMTrial(triali, : ,:) = dynMatrx; 
        dyNB(triali,:) = sum(dynMatrx, 'all');
    end

    allDynMatrix{subji,:} = dynMTrial; 

    allSClust{subji,:} = clustL;
    allSClustB{subji,:} = dyNB; 

end

%% plot example dyn matrix

d2p = squeeze(allDynMatrix{15}(5, :, :));

figure; imagesc(flipud(d2p)); axis square


%% 

meanDYM = cellfun(@(x) squeeze(mean(x)), allDynMatrix, 'un', 0); 
nSubj = size(meanDYM, 1);
lM = size(meanDYM{1},1 ); 
mDYMS = cat(3, meanDYM{:}); mDYMS = reshape(mDYMS, [lM lM nSubj]); mDYMS = permute(mDYMS, [3 1 2]);

d2p = squeeze(mean(mDYMS))

imagesc(flipud(d2p)); axis square

%d2p(triu(d2p) == 0) = nan
%contourf(d2p, 100, 'linecolor', 'none'); axis square; colorbar
set(gca, 'clim', [.2 .45]); colorbar;

%% 

[h p ci ts] = ttest(mDYMS); 
h = squeeze(h); t = squeeze(ts.tstat);
figure(); 
imagesc(flipud(h)); axis square; colorbar


%% 
%contourf(meanDYM{6}); axis square; colorbar

imagesc(flipud(meanDYM{22})); axis square; colorbar



%% compute PS (on diagonal, off diagonal , full matrix) 

clearvars -except DISC_EE DIDC_EE out_c allSClust allSClustB region

timeL = [6:15]; 

C1 = DISC_EE; 

% exclude subjects
if strcmp(region, 'vvs')
    sub2exc = [18 22];
elseif strcmp(region, 'pfc')
    sub2exc = [1];
elseif strcmp(region, 'hipp')
    sub2exc = [2]
end
if sub2exc > 0
    C1(sub2exc) = []; 
end

for subji = 1:size(C1, 1)
    
    clear psT

    C1Sub = C1{subji}; 
    
    for triali = 1:size(C1Sub, 1)

        C1ST = squeeze(C1Sub(triali, :, :));

        % % % entire matrix 
%          C1ST =  triu(C1ST.',1) + tril(C1ST ); 
%          C1ST = C1ST(timeL, timeL);
%          psT(triali, :) = mean(C1ST, 'all'); 

        % % % % only off-diagonal
%          C1ST = squeeze(C1Sub(triali, timeL, timeL));
%          C1ST(eye(size(C1ST))>0) = nan; 
%          psT(triali, :) = mean(C1ST, 'all', 'omitnan'); 

        
        % % % only diagonal
       C1STD = diag(C1ST(timeL, timeL)); % restricted to the period of significant dynamicity
       psT(triali, :) = mean(C1STD); 

    end

    allPST{subji,:} = psT; 


end



%% correlate dynamicity with PS

clear allR allSlopes

tiledlayout(6, 5); set(gcf, 'Position', [10 10 1000 1000])
for subji = 1:size(allPST, 1)
    clear dy 

    set(gca, 'xlim', [-.15 .15], 'ylim', [-10 220])
    ax= nexttile

    ps = allPST{subji}; 
    dy = allSClustB{subji}; 


    allR(subji,:) = corr(ps, dy, 'type', 's');


    scatter(ps, dy,'.')
    h2 = lsline(ax);h2.LineWidth = 2;h2.Color = [.5 .5 .5 ]
    B = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
    allSlopes(subji, :) = B(2);
    allIntercepts(subji, :) = B(1);




end


[h p ci ts] = ttest(allR); 
obsT = ts.tstat



%% plot one bar
clear data
ylim = [-1 1];


data.data = allR;
figure(); set (gcf, 'Position', [300 300 520 650]);
mean_S = mean(data.data, 1);std_S = std(data.data, [], 1); h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on; 
set(hb, 'Marker', '.', 'MarkerSize',60);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);set(hb,'linestyle','none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{''}, 'FontSize', 30, 'linew',2, 'ylim', ylim, 'xlim', [0 2]);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);
ylabel('Rho')
%exportgraphics(gcf, 'myP.png', 'Resolution', 300)





%% compute after shuffling the trial labels 

nPerm = 100

for permi = 1:nPerm

    for subji = 1:size(allPST, 1)
        ps = allPST{subji}; 
        dy = allSClust{subji}; 
    
        psPerm = ps(randperm(length(ps))); 
        allRPerm(subji) = corr(psPerm, dy, 'type', 's'); 
        
    end
    [h p ci ts] = ttest(allRPerm); 
    allTPerm(permi,:) = ts.tstat; 
end






















































%%