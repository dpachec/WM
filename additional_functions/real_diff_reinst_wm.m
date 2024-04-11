
function [out_real] = real_diff_reinst (cfg)

all_cond1 = cfg.all_cond1; all_cond2 = cfg.all_cond2; 
n_subj = size(all_cond1, 1);

for subji = 1:n_subj
    junts = cat(1, all_cond1{subji}, all_cond2{subji});
    labels = [zeros(1, size(all_cond1{subji}, 1))  ones(1, size(all_cond2{subji}, 1))];
    meanReal_cond1(subji, :, :) = squeeze(mean(junts(labels==0, :, :), 1, 'omitnan')); 
    meanReal_cond2(subji, :, :) = squeeze(mean(junts(labels==1, :, :), 1, 'omitnan')); 
end

[sigMH_real subjigMP_real subjigMCI_real sigMT_real] = ttest (meanReal_cond1, meanReal_cond2, 'alpha', cfg.alpha);
sigMH_real = squeeze(sigMH_real);
 if cfg.cond1(6:8) == 'EM1'
     sigMH_real(1:15, 1:5) = 0; % % % remove clusters before baseline
 end
subjigMP_real = squeeze(subjigMP_real); 
sigMT_real = squeeze (sigMT_real.tstat);

sigMH_real(isnan(sigMH_real)) = 0; %for the half-matrix analysubjis
clustInfoReal = bwconncomp(sigMH_real); 
numPixReal = cellfun(@numel,clustInfoReal.PixelIdxList);
[max_clust_infoR maxiR] = max([ cellfun(@numel,clustInfoReal.PixelIdxList) ]); % the zero accounts for empty maps


%%get sums of t's
for ci = 1:length(clustInfoReal.PixelIdxList)
    all_clust_tsum_real(ci) = sum (sigMT_real(clustInfoReal.PixelIdxList{ci}));
end
if exist ('all_clust_tsum_real')
    all_clust_tsum_real = all_clust_tsum_real';
    max_clust_sum_real = max(all_clust_tsum_real);
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









