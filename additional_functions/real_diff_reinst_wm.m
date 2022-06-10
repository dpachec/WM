
function [out_real] = real_diff_reinst (real_cfg)

cfg = real_cfg;
all_cond1 = cfg.all_cond1; all_cond2 = cfg.all_cond2; 
out_real = [];
% > first determine real labels
n_subj = length(all_cond1);

clear junts labels meanReal_cond1 meanReal_cond2;
for si = 1:n_subj
    if cfg.pT == 1 % select random samples from the condition with more trials
            ntrlsC1 = size(all_cond1{1,si},1); ntrlsC2 = size(all_cond2{1,si},1);
            if ntrlsC1 > ntrlsC2
                %disp ('more trials in C1')
                x = randsample(size(all_cond1{1,si},1),size(all_cond2{1,si},1));
                all_cond1{1,si} = all_cond1{1,si}(x,:,:);
            else
                %disp ('more trials in C2')
                x = randsample(size(all_cond2{1,si},1),size(all_cond1{1,si},1));
                all_cond2{1,si} = all_cond2{1,si}(x,:,:);
            end
    end
    junts = cat(1, all_cond1{1,si}, all_cond2{1,si});
    labels = [zeros(1, size(all_cond1{1, si}, 1))  ones(1, size(all_cond2{1, si}, 1))];
    %real labels
    meanReal_cond1(si, :, :) = squeeze(mean(junts(labels==0, :, :), 1, 'omitnan')); 
    meanReal_cond2(si, :, :) = squeeze(mean(junts(labels==1, :, :), 1, 'omitnan')); 
    %meanReal_cond1(si, :, :) = squeeze(mean(junts(labels==0, :, :), 1)); 
    %meanReal_cond2(si, :, :) = squeeze(mean(junts(labels==1, :, :), 1)); 
end

%%real cluster, sum of t's, number of pixels, etc
%check for normality
%     for m = 1:size(meanReal_cond1, 2)
%         for n = 1:size(meanReal_cond1, 3)
%             %meanReal_cond1(:,m,n)
%             [hk1(:,m,n)] = kstest(meanReal_cond1(:,m,n));
%             [hk2(:,m,n)] = kstest(meanReal_cond1(:,m,n));
%         end
%     end    
%ttest
if strcmp (cfg.test2use, 'ttest')
    [sigMH_real sigMP_real sigMCI_real sigMT_real] = ttest (meanReal_cond1, meanReal_cond2, 'alpha', cfg.pval);
    sigMH_real = squeeze(sigMH_real);sigMP_real = squeeze(sigMP_real); 
    sigMT_real = squeeze (sigMT_real.tstat);
end

% %t-test in for loop, to check that it works for wilcox, use matrix in
% code above 
% for m = 1:size(meanReal_cond1, 2)
%     for n = 1:size(meanReal_cond1, 3)
%         %meanReal_cond1(:,m,n)
%         [sigMH_real(:,m,n) sigMP_real(:,m,n) ci t] = ...
%                 ttest (meanReal_cond1(:,m,n), meanReal_cond2(:,m,n), 'alpha', cfg.pval);
%         sigMT_real(:,m,n) = t.tstat;
%     end
% end
% sigMP_real = squeeze(sigMP_real ); sigMH_real  = squeeze(sigMH_real); 
% sigMT_real  = squeeze(sigMT_real);

% %wilcox
if strcmp (cfg.test2use, 'wilcox')
    for m = 1:size(meanReal_cond1, 2)
        for n = 1:size(meanReal_cond1, 3)
            %meanReal_cond1(:,m,n)
            [sigMP_real(:,m,n) sigMH_real(:,m,n) stats] = ...
                    ranksum(meanReal_cond1(:,m,n), meanReal_cond2(:,m,n), 'alpha', cfg.pval);
            sigMT_real(:,m,n) = stats.ranksum/100;
        end
    end
    sigMP_real = squeeze(sigMP_real ); sigMH_real  = squeeze(sigMH_real); 
    sigMT_real  = squeeze(sigMT_real);
end


sigMH_real(isnan(sigMH_real)) = 0; %for the half-matrix analysis
clustInfoReal = bwconncomp(sigMH_real, cfg.connectivity); 
numPixReal = cellfun(@numel,clustInfoReal.PixelIdxList);
[max_clust_infoR maxiR] = max([ cellfun(@numel,clustInfoReal.PixelIdxList) ]); % the zero accounts for empty maps


%%get all sums of t's
clear all_clust_tsum_real;
for ci = 1:length(clustInfoReal.PixelIdxList)
    %disp('entering');
    all_clust_tsum_real(ci) = sum (sigMT_real(clustInfoReal.PixelIdxList{ci}));
end
if exist ('all_clust_tsum_real')
    all_clust_tsum_real = all_clust_tsum_real';
    all_clust_tsum_real(:,2) =  numPixReal;

%max_clust_sum_real = sum (sigMT_real(clustInfoReal.PixelIdxList{maxiR}));
max_clust_sum_real = max(all_clust_tsum_real(:, 1));

clustInfoReal.PixelIdxList = clustInfoReal.PixelIdxList';
else
    max_clust_sum_real = [];
    clustInfoReal.PixelIdxList = [];
    all_clust_tsum_real = [];
end

%fill struct
out_real.max_clust_sum_real     =   max_clust_sum_real;
out_real.meanReal_cond1         =   meanReal_cond1;
out_real.meanReal_cond2         =   meanReal_cond2;
out_real.sigMH_real             =   sigMH_real;
out_real.sigMT_real             =   sigMT_real;
out_real.clustInfoReal          =   clustInfoReal;
out_real.numPixReal             =   numPixReal;
out_real.all_clust_tsum_real    =   all_clust_tsum_real;

end









