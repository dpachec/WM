%%
%%

clear, clc

region = 'vvs';
paths = load_paths_WM(region);

contrasts = {
              'DISC_EE' 'DIDC_EE';
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


%% Compute dynamicity scores at the group level 
clearvars -except DISC_EE DIDC_EE out_c



timeL =  [51:150]; 


meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
nSubj = size(out_c{1}, 1);
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

            
            if h1 & h2 
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
%imagesc(dynMatrx)
contourf(dynMatrx); axis square

%% 
d2p = mean(dynMatrx, 1); 

figure
plot(d2p)





%% permutations 

nPerm = 100;

clear dynMatrxPerm mDiffPerm

for permi = 1:nPerm
    
    permi 
    meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
    meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
    mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
    mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
    nSubj = size(out_c{1}, 1);
    lM = size(meanC1{1},1 ); 
    mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [lM lM nSubj]); mC1FD = permute(mC1FD, [3 1 2]);
    mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [lM lM nSubj]); mC2FD = permute(mC2FD, [3 1 2]);

    mDiff = mC1FD - mC2FD; 
    mDiff = mDiff(:, 1:15, 1:15);

    % % % shuffle matrix for all subjects 
    for subji = 1:nSubj
        mDiffSubj = squeeze(mDiff(subji, :, :)); 
        onDiag = diag(mDiffSubj); 
        for timei = 1:size(mDiff, 2)
            idR = randi(size(mDiff, 2)); 
            tmp = mDiffSubj(timei, idR); 
            onDiagPerm(timei) = tmp; 
            mDiffSubj(timei, idR) = diag(timei); 
            %mDiffSubj(timei, :) = diag(timei); 
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
               
                
                if h1 & h2 
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
    clusLPerm(permi,:) = max(cellfun(@length, clustinfo.PixelIdxList));

end



%% 
d2p = squeeze(dynMatrxPerm(22,:,:));
figure
contourf(d2p)
clustinfo = bwconncomp(d2p)












%% Compute dynamicity scores at the trial level 

clearvars -except DISC_EE DIDC_EE out_c

timeL = [1:150]; 

C1 = DISC_EE; 
C2 = DIDC_EE; 


for subji = 1:length(C1)

    clear clustL
    C1Sub = C1{subji}; 
    C2Sub = C2{subji}; 

    % % % compute mean C2 and substract same mean matrix to all trials ? ? ? 
    meanC2 = mean(C2Sub); 
    nTrials = size(C1Sub, 1); 
    meanC2 = repmat(meanC2, nTrials, 1, 1);
    C1Sub = C1Sub - meanC2; 

    for triali = 1:size(C1Sub, 1)

        C1ST = squeeze(C1Sub(triali, :, :));
        C1ST =  triu(C1ST.',1) + tril(C1ST ); 
        C1ST = C1ST(timeL, timeL);
           
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
    end

    allDynMatrix{subji,:} = dynMTrial; 

    allSClust{subji,:} = clustL;

end

%% 

meanDYM = cellfun(@(x) squeeze(mean(x)), allDynMatrix, 'un', 0); 
nSubj = size(meanDYM, 1);
lM = size(meanDYM{1},1 ); 
mDYMS = cat(3, meanDYM{:}); mDYMS = reshape(mDYMS, [lM lM nSubj]); mDYMS = permute(mDYMS, [3 1 2]);

d2p = squeeze(mean(mDYMS))

%imagesc(d2p); axis square

d2p(triu(d2p) == 0) = nan
contourf(d2p, 100, 'linecolor', 'none'); axis square; colorbar
set(gca, 'clim', [.2 .45]); 

%% 

[h p ci ts] = ttest(mDYMS); 
h = squeeze(h); t = squeeze(ts.tstat);
imagesc(h); axis square; colorbar


%% 
contourf(meanDYM{6}); axis square; colorbar



%% compute PS (on diagonal, off diagonal , full matrix) 

clearvars -except DISC_EE DIDC_EE out_c allSClust

timeL = [6:15]; 

C1 = DISC_EE; 

for subji = 1:size(C1, 1)
    
    clear psT

    C1Sub = C1{subji}; 
    
    for triali = 1:size(C1Sub, 1)

        C1ST = squeeze(C1Sub(triali, :, :));

        % % % entire matrix 
%         C1ST =  triu(C1ST.',1) + tril(C1ST ); 
%         C1ST = C1ST(timeL, timeL);
%         psT(triali, :) = mean(C1ST, 'all'); 

        % % % % only off-diagonal
        C1ST = squeeze(C1Sub(triali, timeL, timeL));
        C1ST(eye(size(C1ST))>0) = nan; 
        psT(triali, :) = mean(C1ST, 'all', 'omitnan'); 

        
        % % % only diagonal
        %C1STD = diag(C1ST(timeL, timeL)); % restricted to the period of significant dynamicity
        %psT(triali, :) = mean(C1STD); 

    end

    allPST{subji,:} = psT; 


end



%% correlate dynamicity with PS


for subji = 1:size(allPST, 1)

    ps = allPST{subji}; 
    dy = allSClust{subji}; 


    allR(subji) = corr(ps, dy, 'type', 's');


end


[h p ci ts] = ttest(allR); 
obsT = ts.tstat



%% plot one bar
clear data
ylim = [-1 1];


data.data = allR';
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