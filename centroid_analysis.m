%% (load only part of big file)
tic
clear
it = 5
f = 6
example = matfile('rAll_spear_noAV_10ms')
tmp = example.rAll(f,it);
rAll{f,it} = tmp{1};
toc

%% create RDM 
clearvars 
f2sav       = 'RNN_vvs_MaintAll_noAv_B_Real_nT_';
f           = 3:54; %in case B
load_parameters_WMFRA;



all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    f2t = strsplit(f2sav, '_'); 
    if strcmp(f2t{5}, '54') | strcmp(f2t{5}, '150') 
        for freqi = 1:length(freqs2test)
            f  = freqs2test(freqi);
            rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
        end
    elseif strcmp(f2t{5}, 'B') 
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
        if subji == 1 % change name only once
            f2sav = [f2sav num2str(f(1)) '-' num2str(f(end)) 'Hz'];
        end
    end
    
    rAll{subji} = rdm_prev;
    
end

cd (currentF)

save(f2sav, 'rAll');


%% CCI (no centroid)

clearvars -except rAll
nTimes = size(rAll{1}.rdm, 3);
nSubj = length(rAll);
md2CTR = zeros(nSubj, nTimes, 6); %mean distance to centroid for each category
md2ARE = zeros(nSubj, nTimes, 6); %cluster area for each cat
md2_TR = zeros(nSubj, nTimes, 4, 6); %split into different conditions
areaCtr = zeros(nSubj, nTimes); %are formed by the centroids


tic
for subji = 1:nSubj
    
    rdm = rAll{subji}.rdm;
    ids = rAll{subji}.ids;
    
    nTrials = length(ids);
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
    %ids0 = double(string(ids0));
    %ind = floor(ids0/100);
    

    id19 = cellfun(@(x) x(19), ids, 'UniformOutput', false);
    id19 = double(string(id19));
    id20 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
    id20 = double(string(id20));
    id_I = id19 == 0 & id20 == 0;
    id_CI = id19 == 1;
    id_CC = id20 == 1;   
    id_CCNI = id19 == 0 & id20 == 1;
    disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])
    nTrls(subji, 1) = sum(id_CI);nTrls(subji, 2) = sum(id_CC); nTrls(subji, 3) = sum(id_I);nTrls(subji, 4) = sum(id_CCNI);
    
    M = zeros(length(ids0));
    ids0 = string(ids0);
    for i = 1:length(M)
        for j = 1:length(M)
            if strcmp(ids0{i}(1), ids0{j}(1)) 
                M(i, j) = 1;
            end
        end 
    end    
    
    for timei = 1:nTimes
         rdmMDS = 1- squeeze(rdm(:,:,timei));
         nDims = size(rdmMDS,2);
         rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
         mWithinCorrect(subji,timei)   = mean(rdmMDS(M == 1 & id_CC == 1), 'all', 'omitnan');
         mAcrossCorrect(subji,timei)   = mean(rdmMDS(M == 0 & id_CC == 1), 'all', 'omitnan');
         mWithinInCorrect(subji,timei) = mean(rdmMDS(M == 1 & id_I == 1), 'all', 'omitnan');
         mAcrossInCorrect(subji,timei) = mean(rdmMDS(M == 0& id_I == 1), 'all', 'omitnan');
         mWithin(subji,timei) = mean(rdmMDS(M == 1), 'all', 'omitnan');
         mAcross(subji,timei) = mean(rdmMDS(M == 0), 'all', 'omitnan');

    end
    
    
    
end

toc

%% CCI (no centroid) ALL TRIALS

clearvars -except rAll
nTimes = size(rAll{1}.rdm, 3);
nSubj = length(rAll);
md2CTR = zeros(nSubj, nTimes, 6); %mean distance to centroid for each category
md2ARE = zeros(nSubj, nTimes, 6); %cluster area for each cat
md2_TR = zeros(nSubj, nTimes, 4, 6); %split into different conditions
areaCtr = zeros(nSubj, nTimes); %are formed by the centroids


tic
for subji = 1:nSubj
    
    rdm = rAll{subji}.rdm;
    ids = rAll{subji}.ids;
    
    nTrials = length(ids);
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
    %ids0 = double(string(ids0));
    %ind = floor(ids0/100);
    

    id19 = cellfun(@(x) x(19), ids, 'UniformOutput', false);
    id19 = double(string(id19));
    id20 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
    id20 = double(string(id20));
    id_I = id19 == 0 & id20 == 0;
    id_CI = id19 == 1;
    id_CC = id20 == 1;   
    id_CCNI = id19 == 0 & id20 == 1;
    disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])
    nTrls(subji, 1) = sum(id_CI);nTrls(subji, 2) = sum(id_CC); nTrls(subji, 3) = sum(id_I);nTrls(subji, 4) = sum(id_CCNI);
    
    M = zeros(length(ids0));
    ids0 = string(ids0);
    for i = 1:length(M)
        i3i = strsplit(ids0{i}, '_'); 
        for j = 1:length(M)
            i3j = strsplit(ids0{j},'_'); 
            junts = [floor(double(string(i3i)) / 100) floor(double(string(i3j)) / 100)];
            %[c1 c2] = unique(junts);
%            xi = setdiff(1:6, c2); 
%             cath = floor (double(string(junts(xi))) / 100);
            %length(unique(junts)) 
            if length(unique(junts)) == 4
                M(i, j) = 1;
            end
        end 
    end    
    
    
    for timei = 1:nTimes
         rdmMDS = 1- squeeze(rdm(:,:,timei));
         nDims = size(rdmMDS,2);
         rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
         mWithinCorrect(subji,timei)   = mean(rdmMDS(M == 1 & id_CC == 1), 'all', 'omitnan');
         mAcrossCorrect(subji,timei)   = mean(rdmMDS(M == 0 & id_CC == 1), 'all', 'omitnan');
         mWithinInCorrect(subji,timei) = mean(rdmMDS(M == 1 & id_I == 1), 'all', 'omitnan');
         mAcrossInCorrect(subji,timei) = mean(rdmMDS(M == 0& id_I == 1), 'all', 'omitnan');
         mWithin(subji,timei) = mean(rdmMDS(M == 1), 'all', 'omitnan');
         mAcross(subji,timei) = mean(rdmMDS(M == 0), 'all', 'omitnan');

    end
    
    
    
end

toc

%% category cluster index for all trials 

sub2exc = [18 22];
catClst = (mAcross - mWithin) ./ (mWithin + mAcross); 

times = 1:size(catClst, 2);

catClst(sub2exc, :) = [];

figure(1); set(gcf, 'Position', [100 100 600 450]);
m1 = mean(catClst, 'omitnan');
std1 = std(catClst,[], 1, 'omitnan');
se1 = std1/sqrt(size(catClst,1));
shadedErrorBar(times, m1,se1, 'k', 1); hold on; 
plot([5 5], get(gca, 'ylim'), 'k:', 'LineWidth', 4);
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 4);

set(gca, 'FontSize', 20);%'xlim', [1 15]
set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
exportgraphics(gcf, 'test1.png', 'Resolution', 300)

%%
mbl = mean(catClst(:,1:5),'all');  
[h p ci tstats] = ttest(catClst);
t = tstats.tstat; 

clear allSTs
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs= allSTs(id)

%%
nTrls_A = nTrls;
sub2exc1 = nTrls_A(:, 3) < 12;
sub2exc2 = zeros(28,1); sub2exc2([18 22]) = 1;sub2exc2 = logical(sub2exc2);
sub2exc = sub2exc1|sub2exc2

catClstCorr = (mAcrossCorrect - mWithinCorrect) ./ (mWithinCorrect + mAcrossCorrect); 
catClstInc = (mAcrossInCorrect - mWithinInCorrect) ./ (mWithinInCorrect + mAcrossInCorrect); 


catClstCorr(sub2exc, :) = [];
catClstInc(sub2exc, :) = [];
figure(1); set(gcf, 'Position', [100 100 600 450]);
m1 = mean(catClstCorr, 'omitnan');
std1 = std(catClstCorr,[], 1, 'omitnan');
se1 = std1/sqrt(size(catClstCorr,1));
shadedErrorBar(times, m1,se1, 'k', 1); hold on; 
m1 = mean(catClstInc, 'omitnan');
std1 = std(catClstInc,[], 1, 'omitnan');
se1 = std1/sqrt(size(catClstInc,1));
shadedErrorBar(times, m1,se1, 'r', 1); hold on; 

diff = catClstCorr - catClstInc;
[hObs pObs ciOBs t] = ttest(diff);
tObs = t.tstat;
clear allObsClust
clustinfo = bwconncomp(hObs);
for pxi = 1:length(clustinfo.PixelIdxList)
   allObsClust(pxi) = sum(tObs(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

bL = get(gca, 'ylim'); bL = bL(1);
h = hObs; 
h(h==0) = nan; h(h==1) = bL;
plot (times, h, 'LineWidth', 5);hold on; 

plot([5 5], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 2);
set(gca, 'FontSize', 20);%'xlim', [1 15]
%set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
exportgraphics(gcf, 'test1.png', 'Resolution', 300)




%% Frequency_resolved CCI analysis

clearvars 
f2sav       = 'rAll_vvs_M&C_noAv_B_Real_nT_'
load_parameters_WMFRA;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    parfor freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1);

        rdm = rdm_prev.rdm;
        ids = rdm_prev.ids;

        nTrials = length(ids);

        ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
        ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
        ids0 = double(string(ids0));
        ind = floor(ids0/100);

        M = zeros(length(ids0));
        ids0 = string(ids0);
        for i = 1:length(M)
            for j = 1:length(M)
                if strcmp(ids0{i}(1), ids0{j}(1)) 
                    M(i, j) = 1;
                end
            end 
        end    

        for timei = 1:nTimes
             rdmMDS = 1- squeeze(rdm(:,:,timei));
             nDims = size(rdmMDS,2);
             rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
             mWithin(subji,freqi, timei) = mean(rdmMDS(M == 1), 'all', 'omitnan');
             mAcross(subji,freqi, timei) = mean(rdmMDS(M == 0), 'all', 'omitnan');
        end

    end
    
end

toc

cd (currentF)

save(f2sav, 'mWithin', 'mAcross');

%% 

sub2exc = [1];
catClst = (mAcross - mWithin) ./ (mWithin + mAcross); 
times = 1:size(catClst, 2);

catClst(sub2exc, :,:) = [];
figure(1)
m1 = squeeze(mean(catClst, 'omitnan'));
imagesc(m1)



%%
%mbl = mean(catClst(:,1:5),'all');  
[h p ci tstats] = ttest(catClst);
h = squeeze(h); t = squeeze(tstats.tstat); 

clear allSTs
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs= allSTs(id)

times = 1:size(catClst, 3)*10;
freqs = 1:540;
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);




%% Trial based centroid analysis 

clearvars -except rAll
nTimes = size(rAll{1}.rdm, 3);
nSubj = length(rAll);
md2CTR = zeros(nSubj, nTimes, 6); %mean distance to centroid for each category
md2ARE = zeros(nSubj, nTimes, 6); %cluster area for each cat
md2_TR = zeros(nSubj, nTimes, 4, 6); %split into different conditions
areaCtr = zeros(nSubj, nTimes); %are formed by the centroids


tic
for subji = 1:nSubj
    
    rdm = rAll{subji}.rdm;
    ids = rAll{subji}.ids;
    
    nTrials = length(ids);
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
    ids0 = double(string(ids0));
    ind = floor(ids0/100);

    id19 = cellfun(@(x) x(19), ids, 'UniformOutput', false);
    id19 = double(string(id19));
    id20 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
    id20 = double(string(id20));
    id_I = id19 == 0 & id20 == 0;
    id_CI = id19 == 1;
    id_CC = id20 == 1;   
    id_CCNI = id19 == 0 & id20 == 1;
    disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])
    nTrls(subji, 1) = sum(id_CI);nTrls(subji, 2) = sum(id_CC); nTrls(subji, 3) = sum(id_I);nTrls(subji, 4) = sum(id_CCNI);
    

    % % match trial numbers for correct and incorrect by randomly removing%from correct trials
%     nnIDCC = find(id_CC);
%     x = randsample(length(nnIDCC),sum(id_I));
%     id_CC_new = zeros(length(id_CC), 1);
%     id_CC_new(nnIDCC(x)) = 1; 
%     id_CC = logical(id_CC_new);



    for timei = 1:nTimes
        
%         % % % MDS
          %d2p = 1- squeeze(rdm(:,:,timei));
          %[rdmMDS] = cmdscale(d2p);
          %nDims = 3; %size(rdmMDS,2);
%         
%             % % % pca
%             d2p = squeeze(rdm(:,:,timei));
%             [rdmMDS,score,latent,tsquared,explained,mu] = pca(d2p);
%             B = cumsum(explained);
%             nDims =  min(find (B > 99.9));
%             dim2check(subji, timei) = nDims;

        % % % % no reduction
         rdmMDS = 1- squeeze(rdm(:,:,timei));
         nDims = size(rdmMDS,2);
         rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal


        %centroid coordinates for each category
        X = rdmMDS(:,1:nDims);
        clear cX
        for i = 1:nDims
            cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
        end
        
% %         % % % % this is just to check
% %         clear i2u % repeat centroid coordinates for each item
% %         for i = 1:nDims
% %             i2u(:,i) = cX(ind, i); 
% %         end
% %         dis2c =  sqrt (sum ( (X - i2u) .^2 , 2));
% %         
        
        %coordinates of the correspondent category centroid are stored in
        %the first dimension of i2O
        clear allCat iO
        allCat = repmat(1:6, nTrials, 1);
        for triali = 1:nTrials
           h = allCat(triali,:) == ind(triali) ;
           iO(triali, :) = allCat(triali,~h) ;
        end
        iO = [ind iO];
        iO = iO(:);

        clear i2O i2oR
        for i = 1:nDims
           i2O(:,i) = cX(iO, i); 
        end
        i2O = reshape (i2O, size(X, 1),[],  nDims);
        
        
        clear dis2c dis2
        for i = 1:6
            i2Otmp = squeeze(i2O(:,i,:));
            dis2c(:,i,:) =  sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
        end
        
        
        % % % %normalize matrix z-score
        %dis2c = ( dis2c - mean(dis2c(:)) ) ./ std(dis2c(:));
        
        % % % % normalize matrix range
        %dis2c = dis2c - min(dis2c(:));
        %dis2c = dis2c ./ max(dis2c(:)); % normalize distances
        
        % %% mean distance to centroid for corresponding category
        d2u = dis2c(:,1);
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TR(subji, timei, 1, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 2, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 3, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 4, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 5, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
        
        % %% mean distance to centroid for non-corresponding categories
        d2u = squeeze(mean(dis2c(:,2:6),2));
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TR(subji, timei, 6, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 7, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 8, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 9, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 10, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
         
         
        
        %centroid area
%        areaCtr(subji, timei, :) = polyarea(cX(:, 1), cX(:, 2)); 
        
%         clear categs categsH areaCat hullCat
%         for i = 1:6 categs{i} = rdmMDS(ind==i,1:2); end 
%         for i = 1:6 if length(categs{i}) > 2 hullCat{i} = convhull(categs{i}); end ; end
%         if length(hullCat) < 6 hullCat{6} = NaN;end
%         for i = 1:6 categsH{i} = categs{i}(hullCat{i}(1:end-1),:,:); end 
%         for i = 1:6 areaCat(i) = polyarea(categsH{i}(:, 1), categsH{i}(:, 2)); end % only 2D area for now
%         md2ARE(subji, timei,:) = areaCat; 
              
        
  
    end
    
    
    
end

toc

%% Trial based centroid analysis ALL TRIALS

clearvars -except rAll
nTimes = size(rAll{1}.rdm, 3);
nSubj = length(rAll);
md2CTR = zeros(nSubj, nTimes, 6); %mean distance to centroid for each category
md2ARE = zeros(nSubj, nTimes, 6); %cluster area for each cat
md2_TR = zeros(nSubj, nTimes, 4, 6); %split into different conditions
areaCtr = zeros(nSubj, nTimes); %are formed by the centroids


tic
for subji = 1:nSubj
    
    rdm = rAll{subji}.rdm;
    ids = rAll{subji}.ids;
    
    nTrials = length(ids);
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids0 = cell2mat(cellfun(@(x) double(strsplit(x(3),'_')), ids, 'UniformOutput', false));
    
    ind = floor(ids0/100);

    id19 = cellfun(@(x) x(19), ids, 'UniformOutput', false);
    id19 = double(string(id19));
    id20 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
    id20 = double(string(id20));
    id_I = id19 == 0 & id20 == 0;
    id_CI = id19 == 1;
    id_CC = id20 == 1;   
    id_CCNI = id19 == 0 & id20 == 1;
    disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])
    nTrls(subji, 1) = sum(id_CI);nTrls(subji, 2) = sum(id_CC); nTrls(subji, 3) = sum(id_I);nTrls(subji, 4) = sum(id_CCNI);
    

    % % match trial numbers for correct and incorrect by randomly removing%from correct trials
%     nnIDCC = find(id_CC);
%     x = randsample(length(nnIDCC),sum(id_I));
%     id_CC_new = zeros(length(id_CC), 1);
%     id_CC_new(nnIDCC(x)) = 1; 
%     id_CC = logical(id_CC_new);



    for timei = 1:nTimes
        
%         % % % MDS
          %d2p = 1- squeeze(rdm(:,:,timei));
          %[rdmMDS] = cmdscale(d2p);
          %nDims = 3; %size(rdmMDS,2);
%         
%             % % % pca
%             d2p = squeeze(rdm(:,:,timei));
%             [rdmMDS,score,latent,tsquared,explained,mu] = pca(d2p);
%             B = cumsum(explained);
%             nDims =  min(find (B > 99.9));
%             dim2check(subji, timei) = nDims;

        % % % % no reduction
         rdmMDS = 1- squeeze(rdm(:,:,timei));
         nDims = size(rdmMDS,2);
         %rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal


        % % % % centroid coordinates for each category
        clear c2t cX
        for cati = 1:6
            c2t = logical(sum(floor(ids0/100) == cati,2)); 
            for i = 1:nDims
                cX(cati, :) =  mean(rdmMDS(c2t,:), 'all', 'omitnan');
            end
        end

        % % for each trial calculate distance to 3 same centroids vs 3 different centroids
        
        tc = floor(ids0/100); 
        for triali = 1:length(rdmMDS)
            catsh2u = tc(triali,:);
            nonCatSh2u = setdiff(1:6, catsh2u);
            sameCatCentrMean = mean(cX(catsh2u,:), 'all');
            diffCatCentrMean = mean(cX(nonCatSh2u,:),'all');
            trialPos = rdmMDS(triali,:);
            mWithin(subji, triali, timei) =  sqrt (sum ( (sameCatCentrMean - trialPos) .^2 , 2, 'omitnan'));
            mAcross(subji, triali, timei) =  sqrt (sum ( (diffCatCentrMean - trialPos) .^2 , 2, 'omitnan'));
        end
        
    end
    
end

toc


%% category cluster index for all trials 
mWithin(mWithin==0) = nan; mAcross(mAcross==0) = nan;
mWithin2 = mean(mWithin,2,'omitnan');
mAcross2 = mean(mAcross,2,'omitnan');

sub2exc = [18 22];
catClst = squeeze((mAcross2 - mWithin2) ./ (mWithin2 + mAcross2)); 

times = 1:size(catClst, 2);

catClst(sub2exc, :) = [];

figure(1); set(gcf, 'Position', [100 100 600 450]);
m1 = mean(catClst, 'omitnan');
std1 = std(catClst,[], 1, 'omitnan');
se1 = std1/sqrt(size(catClst,1));

shadedErrorBar(times, m1,se1, 'k', 1); hold on; 
plot([5 5], get(gca, 'ylim'), 'k:', 'LineWidth', 4);
plot(get(gca, 'xlim'), [0 0], 'k:', 'LineWidth', 4);

set(gca, 'FontSize', 20);%'xlim', [1 15]
%set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
exportgraphics(gcf, 'test1.png', 'Resolution', 300)



%% category cluster index for all trials 
md2_TR_A =md2_TR; 
nTrls_A = nTrls;
sub2exc = [1];
catClst = (squeeze(mean(md2_TR_A(:,:,10,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,5,:), 4, 'omitnan'))) ./ ...
 (squeeze(mean(md2_TR_A(:,:,10,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,5,:), 4, 'omitnan')));
times = 1:size(catClst, 2);

catClst(sub2exc, :) = [];
figure(1)
m1 = mean(catClst, 'omitnan');
std1 = std(catClst,[], 1, 'omitnan');
se1 = std1/sqrt(size(catClst,1));
shadedErrorBar(times, m1,se1, 'm', 1); hold on; 

%%
mbl = mean(catClst(:,1:5),'all');  
[h p ci tstats] = ttest(catClst, mbl);
t = tstats.tstat; 

clear allSTs
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs= allSTs(id)




%% group data analysis (no normalization) 

md2_TR_A =md2_TR; 
nTrls_A = nTrls;
%md2_TR_A([22 28],:,:,:) = []; 
%nTrls_A([22 28],:) = []; 
sub2exc = nTrls_A(:, 3) < 12;
%sub2exc(22) = 1; %22 has only 1 channel


%dCC = squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan'));
%dI = squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan')) ;

%dCC = squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan'));
%dI = squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')) ;

%dCC = squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan'));
%dI = squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan'));
% 
dCC = (squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan'))) ./ ...
 (squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')));
dI = (squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan'))) ./ ... 
  ( squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')));

dCC(sub2exc,:) = []; dI(sub2exc,:) = [];

%times = [-.75:0.01:4.25];
times =1:size(dCC,2);

figure(1)
m1 = mean(dCC, 'omitnan');
std1 = std(dCC,[], 1, 'omitnan');
se1 = std1/sqrt(size(dCC,1));
shadedErrorBar(times, m1,se1, 'm', 1); hold on; 

m2 = mean(dI, 'omitnan');
std2 = std(dI,[], 1, 'omitnan');
se2 = std2/sqrt(size(dI,1));
shadedErrorBar(times, m2,se2, 'b', 1); hold on; 
% 
diff = dCC - dI;
% m3 =  mean(diff, 'omitnan'); 
% std3 = std(diff,[], 1, 'omitnan');
% se3 = std3/sqrt(size(diff,1));
% shadedErrorBar (times,m3, se3, 'y', 1)

[hObs pObs ciOBs t] = ttest(diff);
tObs = t.tstat;
clear allObsClust
clustinfo = bwconncomp(hObs);
for pxi = 1:length(clustinfo.PixelIdxList)
   allObsClust(pxi) = sum(tObs(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
end

bL = get(gca, 'ylim'); bL = bL(1);
h = hObs; 
h(h==0) = nan; h(h==1) = bL;
plot (times, h, 'LineWidth', 3);hold on; 
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)
%plot([0.8 0.8], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)
plot([3.5 3.5], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)
set(gca, 'FontSize', 12)
%set(gca, 'xlim', [-.5 4])
%set(gca, 'xlim', [-.5 1.2])

%exportgraphics(gcf, 'myPlot.png', 'Resolution', 150)





%% PERMUTATIONS
nPerm = 1000;
%clearvars -except rAll
f = 6;
it = 4; 
nTimes = size(rAll{f, it}{1}.rdm, 3);
nSubj = length(rAll{f, it});
md2CTRP = zeros(nPerm, nSubj, nTimes, 6); %mean distance to centroid for each category
md2AREP = zeros(nPerm, nSubj, nTimes, 6); %cluster area for each cat
md2_TRP = zeros(nPerm, nSubj, nTimes, 4, 6); %split into different conditions
areaCtrP = zeros(nPerm, nSubj, nTimes); %are formed by the centroids

tic
clear md2_TR  md2CTR ind_CI ind_I areaCtr id_CCNI

for permi = 1:nPerm
    permi 
    
    for subji = 1:nSubj
        rdm = rAll{f, it}{subji}.rdm;
        ids = rAll{f, it}{subji}.ids;
        
        nTrials = length(ids);
        
        ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
        ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
        ids0 = double(string(ids0));
        ind = floor(ids0/100);

        id19 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
        id19 = double(string(id19));
        id20 = cellfun(@(x) x(21), ids, 'UniformOutput', false);
        id20 = double(string(id20));
        id_I = id19 == 0 & id20 == 0;
        id_CI = id19 == 1;
        id_CC = id20 == 1;   
        id_CCNI = id19 == 0 & id20 == 1;
        %disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])

        
        [idR t] = sort(ind);
        rdm = rdm(t, t, :);
        
        id4p = randperm(length(ind));
        ind = ind(id4p);
        
        
        for timei = 1:nTimes

%             d2p = 1- squeeze(rdm(:,:,timei));
%             [rdmMDS] = cmdscale(d2p);
%             nDims = size(rdmMDS,2);
%             
%             d2p = squeeze(rdm(:,:,timei));
%             [rdmMDS,score,latent,tsquared,explained,mu] = pca(d2p);
%             B = cumsum(explained);
%             nDims =  min(find (B > 99));
        
        % % % % no reduction
          rdmMDS = 1- squeeze(rdm(:,:,timei));
          nDims = size(rdmMDS,2);
          rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal
         



        %centroid coordinates for each category
        X = rdmMDS(:,1:nDims);
        clear cX
        for i = 1:nDims
            cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
        end
        
% %         % % % % this is just to check
% %         clear i2u % repeat centroid coordinates for each item
% %         for i = 1:nDims
% %             i2u(:,i) = cX(ind, i); 
% %         end
% %         dis2c =  sqrt (sum ( (X - i2u) .^2 , 2));
% %         
        
        %coordinates of the correspondent category centroid are stored in
        %the first dimension of i2O
        clear allCat iO
        allCat = repmat(1:6, nTrials, 1);
        for triali = 1:nTrials
           h = allCat(triali,:) == ind(triali) ;
           iO(triali, :) = allCat(triali,~h) ;
        end
        iO = [ind iO];
        iO = iO(:);

        clear i2O i2oR
        for i = 1:nDims
           i2O(:,i) = cX(iO, i); 
        end
        i2O = reshape (i2O, size(X, 1),[],  nDims);
        
        
        clear dis2c dis2
        for i = 1:6
            i2Otmp = squeeze(i2O(:,i,:));
            dis2c(:,i,:) = sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
        end
        
        
        % %% mean distance to centroid for corresponding category
        d2u = dis2c(:,1);
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TRP(permi,subji, timei, 1, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 2, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 3, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 4, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 5, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
        
        % %% mean distance to centroid for non-corresponding categories
        d2u = squeeze(mean(dis2c(:,2:6),2));
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TRP(permi,subji, timei, 6, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 7, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 8, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi,subji, timei, 9, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TRP(permi, subji, timei, 10, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
         




%             areaCtrP(permi, subji, timei, :) = polyarea(cX(:, 1), cX(:, 2)); 
% 
%             clear categs categsH areaCat hullCat
%             for i = 1:6 categs{i} = rdmMDS(ind==i,1:2); end 
%             for i = 1:6 if length(categs{i}) > 2 hullCat{i} = convhull(categs{i}); end ; end
%             if length(hullCat) < 6 hullCat{6} = NaN;end
%             for i = 1:6 categsH{i} = categs{i}(hullCat{i}(1:end-1),:,:); end 
%             for i = 1:6 areaCat(i) = polyarea(categsH{i}(:, 1), categsH{i}(:, 2)); end % only 2D area for now
%             md2AREP(permi, subji, timei,:) = areaCat; 

% 
% 
% 
%             md2CTRP(permi, subji, timei,:) = accumarray(ind, dis2c, [6 1], @mean); %mean distance to mean in all trials from same category




        end

    end
    
end

filename = ['myPermFile_maint_' num2str(nPerm)];
save(filename, 'md2_TRP');

toc


%% check one perm (no normalization)
p = 4;
md2_TR_A = squeeze(md2_TRP(p,:,:,:,:)); 
sub2exc = nTrls(:, 3) < 10; 


bline = 25:175; 


dCC = (squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan'))) ./ ...
      (squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')));
  
dI = (squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan'))) ./ ... 
       ( squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')));

dCC(sub2exc,:) = []; dI(sub2exc,:) = [];


times = [-.75:0.01:4.25];
figure()
m1 = mean(dCC, 'omitnan');
std1 = std(dCC,[], 1, 'omitnan');
se1 = std1/sqrt(size(dCC,1));
shadedErrorBar(times, m1,se1, 'm', 1); hold on; 

m2 = mean(dI, 'omitnan');
std2 = std(dI,[], 1, 'omitnan');
se2 = std2/sqrt(size(dI,1));
shadedErrorBar(times, m2,se2, 'b', 1); hold on; 

diff = dCC - dI;
h = ttest(diff);

bL = get(gca, 'ylim'); bL = bL(1);
h(h==0) = nan; h(h==1) = bL;
plot (times, h, 'LineWidth', 3);hold on; 
plot([0 0], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)
%plot([0.8 0.8], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)
plot([3.5 3.5], get(gca, 'ylim'), 'k:', 'LineWidth' ,1)

set(gca, 'xlim', [-.5 4])
%set(gca, 'xlim', [-.5 1.2])




%% count t values

restrictedTime = 100:400;  

clear h t allSTs max_clust_sum
for permi = 1:nPerm;
    
md2_TR_A = squeeze(md2_TRP(permi,:,:,:,:)); 

sub2exc = nTrls(:, 3) < 10; 

% % % % no normalization
 dCC = (squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan'))) ./ ...
       (squeeze(mean(md2_TR_A(:,:,1,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,6,:), 4, 'omitnan')));
   
 dI = (squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')) - squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan'))) ./ ... 
        ( squeeze(mean(md2_TR_A(:,:,3,:), 4, 'omitnan')) + squeeze(mean(md2_TR_A(:,:,8,:), 4, 'omitnan')));
 
 dCC(sub2exc,:) = []; dI(sub2exc,:) = [];

diff = dCC - dI;

[h p ci tstats] = ttest(diff);
t = tstats.tstat; 
h = h(restrictedTime);
t = t(restrictedTime);

clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
        end
else
    allSTs(permi, :) = 0;
end

[maxh id] = max(abs(allSTs(permi, :)));
max_clust_sum(permi,:) = allSTs(permi,id);
    
    
end


%%
t1 = [[118.283]];

close all;
figure(1)
%[counts] = hist(max_clust_sum, 12); hold on;
[counts,centers] = hist(max_clust_sum, 14);
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 500] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
export_fig(1, filename,'-transparent', '-r300');
%close all;



%% get rank (for stats)
max_clust_sum_real = t1;
allAb = max_clust_sum(abs(max_clust_sum) > abs(max_clust_sum_real));
p =1 - (nPerm - (length (allAb)))  /nPerm




%% Trial based analysis (ALL ITEMS)

clearvars -except rAll
f = 6;
it = 7; 
nTimes = size(rAll{f, it}{1}.rdm, 3);
nSubj = length(rAll{f, it});
md2CTR = zeros(nSubj, nTimes, 6); %mean distance to centroid for each category
md2ARE = zeros(nSubj, nTimes, 6); %cluster area for each cat
md2_TR = zeros(nSubj, nTimes, 4, 6); %split into different conditions
areaCtr = zeros(nSubj, nTimes); %are formed by the centroids


tic
for subji = 1:nSubj
    
    rdm = rAll{f, it}{subji}.rdm;
    ids = rAll{f, it}{subji}.ids;
    
    nTrials = length(ids);
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);
    ids0 = double(string(ids0));
    ind = floor(ids0/100);

    id19 = cellfun(@(x) x(19), ids, 'UniformOutput', false);
    id19 = double(string(id19));
    id20 = cellfun(@(x) x(20), ids, 'UniformOutput', false);
    id20 = double(string(id20));
    id_I = id19 == 0 & id20 == 0;
    id_CI = id19 == 1;
    id_CC = id20 == 1;   
    id_CCNI = id19 == 0 & id20 == 1;
    disp(['Subj : ' num2str(subji) ' > CI = ' num2str(sum(id_CI))  ' > CC = ' num2str(sum(id_CC))  ' > I = ' num2str(sum(id_I)) ' > CCNI = ' num2str(sum(id_CCNI)) ])
    nTrls(subji, 1) = sum(id_CI);nTrls(subji, 2) = sum(id_CC); nTrls(subji, 3) = sum(id_I);nTrls(subji, 4) = sum(id_CCNI);
    

    % % match trial numbers for correct and incorrect by randomly removing%from correct trials
%     nnIDCC = find(id_CC);
%     x = randsample(length(nnIDCC),sum(id_I));
%     id_CC_new = zeros(length(id_CC), 1);
%     id_CC_new(nnIDCC(x)) = 1; 
%     id_CC = logical(id_CC_new);



    for timei = 1:nTimes
        
%         % % % MDS
          %d2p = 1- squeeze(rdm(:,:,timei));
          %[rdmMDS] = cmdscale(d2p);
          %nDims = 3; %size(rdmMDS,2);
%         
%             % % % pca
%             d2p = squeeze(rdm(:,:,timei));
%             [rdmMDS,score,latent,tsquared,explained,mu] = pca(d2p);
%             B = cumsum(explained);
%             nDims =  min(find (B > 99.9));
%             dim2check(subji, timei) = nDims;

        % % % % no reduction
         rdmMDS = 1- squeeze(rdm(:,:,timei));
         nDims = size(rdmMDS,2);
         rdmMDS(find(eye(size(rdmMDS)))) = nan; % % % remove zero coordinates on the diagonal


        %centroid coordinates for each category
        X = rdmMDS(:,1:nDims);
        clear cX
        for i = 1:nDims
            cX(:, i) =  accumarray(ind, rdmMDS(:,i), [6 1], @mean);    % coordinates of the centroid for each category
        end
        
% %         % % % % this is just to check
% %         clear i2u % repeat centroid coordinates for each item
% %         for i = 1:nDims
% %             i2u(:,i) = cX(ind, i); 
% %         end
% %         dis2c =  sqrt (sum ( (X - i2u) .^2 , 2));
% %         
        
        %coordinates of the correspondent category centroid are stored in
        %the first dimension of i2O
        clear allCat iO
        allCat = repmat(1:6, nTrials, 1);
        for triali = 1:nTrials
           h = allCat(triali,:) == ind(triali) ;
           iO(triali, :) = allCat(triali,~h) ;
        end
        iO = [ind iO];
        iO = iO(:);

        clear i2O i2oR
        for i = 1:nDims
           i2O(:,i) = cX(iO, i); 
        end
        i2O = reshape (i2O, size(X, 1),[],  nDims);
        
        
        clear dis2c dis2
        for i = 1:6
            i2Otmp = squeeze(i2O(:,i,:));
            dis2c(:,i,:) =  sqrt (sum ( (X - i2Otmp) .^2 , 2, 'omitnan'));
        end
        
        
        % % % %normalize matrix z-score
        %dis2c = ( dis2c - mean(dis2c(:)) ) ./ std(dis2c(:));
        
        % % % % normalize matrix range
        %dis2c = dis2c - min(dis2c(:));
        %dis2c = dis2c ./ max(dis2c(:)); % normalize distances
        
        % %% mean distance to centroid for corresponding category
        d2u = dis2c(:,1);
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TR(subji, timei, 1, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 2, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 3, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 4, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 5, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
        
        % %% mean distance to centroid for non-corresponding categories
        d2u = squeeze(mean(dis2c(:,2:6),2));
        dis2c_CC = d2u; dis2c_CI = d2u; dis2c_I = d2u; dis2c_CCNI = d2u; 
        dis2c_CC(~id_CC) = nan;
        dis2c_CI(~id_CI) = nan;
        dis2c_I(~id_I) = nan;
        dis2c_CCNI(~id_CCNI) = nan; 
        
        md2_TR(subji, timei, 6, :) = accumarray(ind, dis2c_CC, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 7, :) = accumarray(ind, dis2c_CI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 8, :) = accumarray(ind, dis2c_I, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 9, :) = accumarray(ind, dis2c_CCNI, [6 1], @(x)mean(x,'omitnan')); 
        md2_TR(subji, timei, 10, :) = accumarray(ind, d2u, [6 1], @(x)mean(x,'omitnan')); % all trials
         
         
        
        %centroid area
%        areaCtr(subji, timei, :) = polyarea(cX(:, 1), cX(:, 2)); 
        
%         clear categs categsH areaCat hullCat
%         for i = 1:6 categs{i} = rdmMDS(ind==i,1:2); end 
%         for i = 1:6 if length(categs{i}) > 2 hullCat{i} = convhull(categs{i}); end ; end
%         if length(hullCat) < 6 hullCat{6} = NaN;end
%         for i = 1:6 categsH{i} = categs{i}(hullCat{i}(1:end-1),:,:); end 
%         for i = 1:6 areaCat(i) = polyarea(categsH{i}(:, 1), categsH{i}(:, 2)); end % only 2D area for now
%         md2ARE(subji, timei,:) = areaCat; 
              
        
  
    end
    
    
    
end

toc


































%% 