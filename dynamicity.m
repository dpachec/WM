%%
clearvars

region = 'pfc'; 
paths = load_paths_WM(region);

contrasts = {
              'DISC_EE' 'DIDC_EE';
             };

c = unique (contrasts);
d = cellfun(@(x) [x '_id'], c, 'un', 0);
for i = 1:length(c) 
    load([paths.currSess c{i} '.mat']);
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


%%
clear mC1
meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [21 21 28]); mC1FD = permute(mC1FD, [3 1 2]);
mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [21 21 28]); mC2FD = permute(mC2FD, [3 1 2]);

mDiff = mC1FD - mC2FD; 


for timei = 6:15

    for timej = 6:15

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
clusL = max(cellfun(@length, clustinfo.PixelIdxList))





%% test it 
figure
imagesc(squeeze(mean(DISC_EE{1})))


figure
imagesc(squeeze(mC1(1,:,:)))


%% 
figure; 
%imagesc(dynMatrx)
contourf(dynMatrx)


%% permutations 

nPerm = 10

clear clusL dynMatrxPerm

for permi = 1:nPerm
    clear dynMatrxPerm

    meanC1 = cellfun(@(x) squeeze(mean(x)), DISC_EE, 'un', 0); 
    meanC2 = cellfun(@(x) squeeze(mean(x)), DIDC_EE, 'un', 0); 
    mC1F  = cellfun(@(x) triu(x.',1) + tril(x), meanC1, 'un', 0); 
    mC2F  = cellfun(@(x) triu(x.',1) + tril(x), meanC2, 'un', 0); 
    mC1FD = cat(3, mC1F{:}); mC1FD = reshape(mC1FD, [21 21 28]); mC1FD = permute(mC1FD, [3 1 2]);
    mC2FD = cat(3, mC2F{:}); mC2FD = reshape(mC2FD, [21 21 28]); mC2FD = permute(mC2FD, [3 1 2]);
    
    mDiff = mC1FD - mC2FD; 
    

    
    for timei = 6:15
    
        for timej = 6:15
    
            if timei~=timej
    
                    offDiag = mDiff(:, timei, timej); 
                    onDiag1 = mDiff(:, timei, timei); 
                    onDiag2 = mDiff(:, timej, timej); 

                    all3 = [onDiag1 onDiag2 offDiag]; 
                    all3 = all3(:, randperm(3));

                    [h1 p1 ci1 ts1] = ttest(all3(:, 1), all3(:,3)); 
                    [h2 p2 ci2 ts2] = ttest(all3(:, 2), all3(:,3)); 
               
                
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
    clusL(permi,:) = max(cellfun(@length, clustinfo.PixelIdxList));

end







%% 
d2p = squeeze(dynMatrxPerm(10,:,:));
figure
imagesc(d2p)
clustinfo = bwconncomp(d2p)

%%
d2p = squeeze(dynMatrxPerm);
figure
imagesc(d2p)











%%


























%%