
function [out_perm] = myPerm (cfg_perm)

    meanReal_cond1      =   cfg_perm.out_real.meanReal_cond1;
    meanReal_cond2      =   cfg_perm.out_real.meanReal_cond2;
    numPixReal          =   cfg_perm.out_real.numPixReal;
    max_clust_sum_real  =   cfg_perm.out_real.max_clust_sum_real;
    sigMH_real          =   cfg_perm.out_real.sigMH_real;
    clustInfoReal       =   cfg_perm.out_real.clustInfoReal;
    
    
    n_subj = size(meanReal_cond1, 1);
    rsa_perm = zeros (cfg_perm.n_perm, n_subj, 2, size(meanReal_cond1, 2), size(meanReal_cond1, 3)); 
    
    for permi = 1:cfg_perm.n_perm
        for subji = 1:n_subj
            junts = cat(1, meanReal_cond1(subji,:,:), meanReal_cond2(subji,:,:));
            labels = [zeros(1, size(meanReal_cond1(subji), 1))  ones(1, size(meanReal_cond2(subji), 1))];
            fake_condition_mapping = labels(randperm(size(labels, 2)));
            meanFake_cond1(subji, :, :) = squeeze(mean(junts(fake_condition_mapping==0, :, :), 1)); 
            meanFake_cond2(subji, :, :) = squeeze(mean(junts(fake_condition_mapping==1, :, :), 1)); 
        end
        rsa_perm(permi, :, 1, :, :) =  meanFake_cond1;
        rsa_perm(permi, :, 2, :, :) =  meanFake_cond2;
        
    end
    
    out_perm.m_fakeH= zeros (cfg_perm.n_perm, length(meanReal_cond1), length(meanReal_cond1));
    out_perm.m_fakeT= zeros (cfg_perm.n_perm, length(meanReal_cond1), length(meanReal_cond1));
    
    for permi = 1:cfg_perm.n_perm
        rsa_all_1 = rsa_perm(permi, :,1,:,:); 
        rsa1_perm = squeeze(rsa_all_1);
        rsa_all_2 = rsa_perm(permi, :,2,:,:); 
        rsa2_perm = squeeze(rsa_all_2);

        [sigMH_permi sigMP_permi sigMCI_permi sigMT_permi] = ttest (rsa1_perm, rsa2_perm, 'alpha', cfg_perm.pval);
        sigMT_permi = squeeze (sigMT_permi.tstat);
        
        % get number of elements in largest supra-threshold cluster
        sigMH_permi(isnan(sigMH_permi)) = 0; %for the half-matrix analysis
        clustinfo = bwconncomp(sigMH_permi);
        [numPixPermi(permi) maxi] = max([0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
        if numPixPermi(permi) > 0
            out_perm.max_clust_sum(permi) = sum (sigMT_permi(clustinfo.PixelIdxList{maxi-1}));
            for pixi = 1:length(clustinfo.PixelIdxList)
                out_perm.all_clust_sum{permi}(pixi) = sum (sigMT_permi(clustinfo.PixelIdxList{pixi}));
            end
            
        else
            %disp (['no significant cluster in permutation ' num2str(permi)]);
            out_perm.max_clust_sum(permi) = 0; 
        end
        
        

        
        %out_perm.allClustInfo{permi} =  clustinfo;
    end
   
    numPixPermi = numPixPermi';
    out_perm.max_clust_sum = out_perm.max_clust_sum';
    
    
    %%threshold to the sum of the significant cluster of t values
    lower_threshold = prctile(out_perm.max_clust_sum,    (cfg_perm.pval*100)/2);
    out_perm.max_clust_sum_perm = prctile(out_perm.max_clust_sum,100-((cfg_perm.pval*100)/2)); % /2 for a 2 sided test
    %out_perm.max_clust_sum_perm = prctile(out_perm.max_clust_sum,100-((cfg_perm.pval*100)));

    %numPixThres = prctile(abs(numPixPermi),100-cfg_perm.pvalC*100);
    %whichclusters2remove = find(numPixReal < numPixThres);
    
    %keep the biggest cluster
    %whichclusters2remove = find(cfg_perm.out_real.all_clust_tsum_real(:,1) < max_clust_sum_real);
    
    %keep based on summed t-values
    if ~isempty(cfg_perm.out_real.all_clust_tsum_real)
        whichclusters2remove = find(cfg_perm.out_real.all_clust_tsum_real(:,1) < out_perm.max_clust_sum_perm);
    else
        whichclusters2remove = [];
    end
    %whichclusters2remove = 1:31 %hack for NaN
    
    out_perm.sigMH_thres = sigMH_real;
    %remove small clusters
    for i=1:length(whichclusters2remove)
        out_perm.sigMH_thres(clustInfoReal.PixelIdxList{whichclusters2remove(i)})=0;
    end
    %manually delete all items
     %load ('clustInfoReal_plot');
     %out_perm.sigMH_thres(clustInfoReal.PixelIdxList{1})=1;
     %out_perm.sigMH_thres(clustInfoReal.PixelIdxList{16})=1;
     %out_perm.sigMH_thres(clustInfoReal.PixelIdxList{17})=0;
     
     
        fprintf('\n');
     
     
        all_clust_tsum_real             =   cfg_perm.out_real.all_clust_tsum_real;
        max_clust_sum_real              =   cfg_perm.out_real.max_clust_sum_real;
        max_clust_sum                   =   out_perm.max_clust_sum;
        max_clust_sum_perm              =   out_perm.max_clust_sum_perm;
        out_perm.max_clust_sum_real     =   cfg_perm.out_real.max_clust_sum_real;
        
        if cfg_perm.savePerm
            filename =  [num2str(cfg_perm.n_perm) 'p_' cfg_perm.cond1 '_' cfg_perm.cond2];
                save (filename, 'all_clust_tsum_real', 'out_perm');
        end

end






