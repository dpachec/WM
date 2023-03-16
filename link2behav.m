

%% extract trial level fits in PFC and compare between 3 different types of items 

clear , clc
f2sav = 'RNN_pfc_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1.mat'; 
cfg = getParams(f2sav); if strcmp(cfg.brainROI, 'vvs')sub2exc = [18 22]; elseif strcmp(cfg.brainROI, 'pfc')sub2exc = [1];end
paths = load_paths_WM(cfg.brainROI);
load([paths.results.trial_level f2sav]);


%%

foI = 14:26
toI = 7:12

for subji = 1:length(nnFit)

    ids = nnFit{subji, 2};
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0);
    
    idIL = cellfun(@(x) x(19), ids, 'un', 0); idIL = double(string(idIL));
    idCL = cellfun(@(x) x(20), ids, 'un', 0);idCL= double(string(idCL));
    id_Incorrect = idCL == 0 & idIL == 0;
    id_CorrIL = idIL == 1;
    id_CorrCL = idCL == 1;   
    id_CorrCNI = idIL == 0 & idCL == 1;
    disp(['Subj : ' num2str(subji) ' > CorrIL = ' num2str(sum(id_CorrIL))  ' > CorrCL = ' num2str(sum(id_CorrCL))  ' > Incorrect = ' ...
            num2str(sum(id_Incorrect)) ' > Corr category but not item = ' num2str(sum(id_CorrCNI)) ])
    nTrls(subji, 1) = sum(id_CorrIL);nTrls(subji, 2) = sum(id_CorrCL); nTrls(subji, 3) = sum(id_Incorrect);nTrls(subji, 4) = sum(id_CorrCNI);

    for layi = 1:size(nnFit{1}, 1)
       
        fit_CIL(subji, layi, :) = mean(atanh(nnFit{subji, 1}(layi,id_CorrIL,foI,toI)), 'all');
        fit_CCL(subji, layi, :) = mean(atanh(nnFit{subji, 1}(layi,id_CorrCL,foI,toI)), 'all');
        fit_Inc(subji, layi, :) = mean(atanh(nnFit{subji, 1}(layi,id_Incorrect,foI,toI)), 'all');
        fit_CINC(subji, layi, :) = mean(atanh(nnFit{subji, 1}(layi,id_CorrCNI,foI,toI)), 'all');
        
    end

end


% % %  stats 
[h p ci ts] = ttest(fit_Inc, fit_CINC); 
t = squeeze(ts.tstat); 
disp(['p = ' num2str(p)])
disp(['t = ' num2str(ts.tstat)])























