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


%% ENCODING 


meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
nSubj = size(out_c{1}, 1);
lM = size(meanC1{1},1 ); 
mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [lM lM nSubj]); mC1FD = permute(mC1FD, [3 1 2]);
mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [lM lM nSubj]); mC2FD = permute(mC2FD, [3 1 2]);

mDiff = mC1FD - mC2FD; 


for timei = 1:lM

    for timej = 1:lM

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
contourf(dynMatrx); 

%% 
d2p = mean(dynMatrx, 1); 

figure
plot(d2p)





%% permutations 

nPerm = 100;

clear clusL dynMatrxPerm mDiffPerm

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
    

    % % % shuffle matrix for all subjects 
    for subji = 1:nSubj
        mDiffSubj = squeeze(mDiff(subji, :, :)); 
        onDiag = diag(mDiffSubj); 
        for timei = 1:lM
            idR = randi(lM); 
            tmp = mDiffSubj(timei, idR); 
            onDiagPerm(timei) = tmp; 
            mDiffSubj(timei, idR) = diag(timei); 
            %mDiffSubj(timei, :) = diag(timei); 
            mDiffSubj(timei, timei) = onDiagPerm(timei); 
        end
        mDiffPerm(subji, :, :) = mDiffSubj ;
    end
    
    for timei = 1:lM
    
        for timej = 1:lM
    
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
d2p = squeeze(dynMatrxPerm(5,:,:));
figure
contourf(d2p)
clustinfo = bwconncomp(d2p)

%%
d2p = squeeze(dynMatrxPerm);
figure
imagesc(d2p)











%%


























%%