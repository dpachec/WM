
function [out_perm] = myPerm (cfg_perm)

    cfg                 =   cfg_perm; 
    meanReal_cond1      =   cfg.out_real.meanReal_cond1;
    numPixReal          =   cfg.out_real.numPixReal;
    max_clust_sum_real  =   cfg.out_real.max_clust_sum_real;
    sigMH_real          =   cfg.out_real.sigMH_real;
    clustInfoReal       =   cfg.out_real.clustInfoReal;
    all_cond1           =   cfg.all_cond1; 
    all_cond2           =   cfg.all_cond2;
    
    
    out_perm        =       [];
    n_subj = length(all_cond1);
    clear fake_condition_mapping juntsFake_cong juntsFake_inc rsa_perm;
    %rsa_perm = zeros (n_perm, n_subj, 2, mlim(end), mlim(end)); 
    rsa_perm = zeros (cfg.n_perm, n_subj, 2, size(meanReal_cond1, 2), size(meanReal_cond1, 3)); 
    disp('Perm:      '); 
    for permi = 1:cfg.n_perm
        if (permi < 10) fprintf('\b'); fprintf('\b'); fprintf('\b'); 
        elseif (permi < 100) fprintf('\b'); fprintf('\b'); fprintf('\b'); fprintf('\b');  
        elseif (permi < 1000) fprintf('\b'); fprintf('\b'); fprintf('\b') ; fprintf('\b'); fprintf('\b'); 
        else fprintf('\b'); fprintf('\b'); fprintf('\b'); fprintf('\b') ; fprintf('\b'); fprintf('\b');  end
        fprintf('%d %s', permi, ' '); 
        
        % > Shuffle labels
        for si = 1:n_subj
            junts = zeros (size(all_cond1{1, si}, 1) + size(all_cond2{1,si},1), size(meanReal_cond1, 2), size(meanReal_cond1, 3));
            
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
            fake_condition_mapping = labels(randperm(size(labels, 2)));
            
            %real labels
            meanFake_cond1(si, :, :) = squeeze(mean(junts(fake_condition_mapping==0, :, :), 1)); 
            meanFake_cond2(si, :, :) = squeeze(mean(junts(fake_condition_mapping==1, :, :), 1)); 
        end
        %size(meanFake_cond1)
        rsa_perm(permi, :, 1, :, :) =  meanFake_cond1;
        rsa_perm(permi, :, 2, :, :) =  meanFake_cond2;
        
    end
    
    out_perm.m_fakeH= zeros (cfg.n_perm, length(meanReal_cond1), length(meanReal_cond1));
    out_perm.m_fakeT= zeros (cfg.n_perm, length(meanReal_cond1), length(meanReal_cond1));
    
    clear max_pixel_vals; clear numPixPermi; clear out_perm.max_clust_sum; 
    for permi = 1:cfg.n_perm
        %fprintf('%d %s', permi, ' '); if (rem(permi,20) == 0) fprintf('\n'); end
        rsa_all_C = rsa_perm(permi, :,1,:,:); rsaC_perm = squeeze(rsa_all_C);
        rsa_all_I = rsa_perm(permi, :,2,:,:); rsaI_perm = squeeze(rsa_all_I);

        [sigMH_permi sigMP_permi sigMCI_permi sigMT_permi] = ttest (rsaC_perm, rsaI_perm, 'alpha', cfg.pval);
        sigMT_permi = squeeze (sigMT_permi.tstat);
        % save maximum pixel values
        max_pixel_pvals(permi,:) = [ min(sigMT_permi(:)) max(sigMT_permi(:)) ];
        m_fakeH(permi, :,:) = sigMH_permi; m_fakeT(permi, :,:) = sigMT_permi;

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
            disp (['no significant cluster in permutation ' num2str(permi)]);
            out_perm.max_clust_sum(permi) = 0; 
        end
        
        

        
        %out_perm.allClustInfo{permi} =  clustinfo;
    end
   
    numPixPermi = numPixPermi';
    out_perm.max_clust_sum = out_perm.max_clust_sum';
    
    
    %%threshold to the sum of the significant cluster of t values
    lower_threshold = prctile(out_perm.max_clust_sum,    (cfg.pval*100)/2);
    out_perm.max_clust_sum_perm = prctile(out_perm.max_clust_sum,100-((cfg.pval*100)/2)); % /2 for a 2 sided test
    %out_perm.max_clust_sum_perm = prctile(out_perm.max_clust_sum,100-((cfg.pval*100)));

    %numPixThres = prctile(abs(numPixPermi),100-cfg.pvalC*100);
    %whichclusters2remove = find(numPixReal < numPixThres);
    
    %keep the biggest cluster
    %whichclusters2remove = find(cfg.out_real.all_clust_tsum_real(:,1) < max_clust_sum_real);
    
    %keep based on summed t-values
    if ~isempty(cfg.out_real.all_clust_tsum_real)
        whichclusters2remove = find(cfg.out_real.all_clust_tsum_real(:,1) < out_perm.max_clust_sum_perm);
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
     
     
        all_clust_tsum_real             =   cfg.out_real.all_clust_tsum_real;
        max_clust_sum_real              =   cfg.out_real.max_clust_sum_real;
        max_clust_sum                   =   out_perm.max_clust_sum;
        max_clust_sum_perm              =   out_perm.max_clust_sum_perm;
        out_perm.max_clust_sum_real     =   cfg.out_real.max_clust_sum_real;
        
        if cfg.savePerm
            filename =  [num2str(cfg.n_perm) 'p_' cfg.cond1 '_' cfg.cond2];
                save (filename, 'all_clust_tsum_real', 'out_perm');
        end

end






