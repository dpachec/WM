function [all_r_Times] = fitModel_WM(neuralRDMs, networkRDMs, fitMode)


if ndims(neuralRDMs) == 4 % if frequency resolved

    nLays  = size(networkRDMs, 1);
    nFreqs = size(neuralRDMs, 3); 
    nTimes = size(neuralRDMs, 4); 

    if fitMode == 0 % trials or no trials
        all_r_Times = zeros(nLays, nFreqs, nTimes);

        for layi = 1:nLays
            M =  squeeze(networkRDMs(layi,:,:)); 
            M = vectorizeRDM(M);
            parfor freqi = 1:nFreqs
                for timei = 1:nTimes
                    rdm = squeeze(neuralRDMs(:, :, freqi, timei));
                    rdm = vectorizeRDM(rdm);
                    allTEst = corr(rdm', M', 'type', 's');
                    all_r_Times(layi,freqi,timei) = allTEst;  
                end
            end
        end

    else

        nTrials = size(neuralRDMs, 1);
        all_r_Times = zeros(nLays,nTrials, nFreqs, nTimes);

        for layi = 1:nLays
            for triali = 1:size(neuralRDMs, 1)
                M =  squeeze(networkRDMs(layi,triali,:)); 
                M(triali) = []; 
                for freqi = 1:nFreqs
                    parfor timei = 1:nTimes
                        rdm = squeeze(neuralRDMs(:, triali, freqi, timei));
                        rdm(triali) = []; 
                        allTEst = corr(rdm, M, 'type', 's');
                        all_r_Times(layi,triali, freqi,timei) = allTEst;  
                    end
                end
            end
        end    


    end
    
    
    
    
    
else % if only 3 dimensions (frequencies used as features)
    
    nLays  = size(networkRDMs, 1);
    nTimes = size(neuralRDMs, 3); 

    if fitMode == 0
        all_r_Times = zeros(nLays, nTimes);

        for layi = 1:nLays
            M =  squeeze(networkRDMs(layi,:,:)); 
            M = vectorizeRDM(M);
            parfor timei = 1:nTimes
                rdm = squeeze(neuralRDMs(:, :, timei));
                rdm = vectorizeRDM(rdm);
                allTEst = corr(rdm', M', 'type', 's');
                all_r_Times(layi,timei) = allTEst;  
            end
        end

    else

    nTrials = size(neuralRDMs, 1);
    all_r_Times = zeros(nLays,nTrials, nTimes);

    for layi = 1:nLays
        for triali = 1:size(neuralRDMs, 1)
            M =  squeeze(networkRDMs(layi,triali,:)); 
            M(triali) = []; 
            parfor timei = 1:nTimes
                rdm = squeeze(neuralRDMs(:, triali, timei));
                rdm(triali) = []; 
                allTEst = corr(rdm, M, 'type', 's');
                all_r_Times(layi,triali, timei) = allTEst;  
            end
        end
    end    


    
    all_r_Times = squeeze(all_r_Times);
    
    
    
    
    
    
    
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end
    
    

end