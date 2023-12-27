function [all_r_Times] = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, fitMode)


if ndims(neuralRDMs) == 4 % if frequency resolved

    nLays  = size(networkRDMs, 1);
    nFreqs = size(neuralRDMs, 3); 
    nTimes = size(neuralRDMs, 4); 

    if fitMode == 0 % trials or no trials
        all_r_Times = zeros(nLays, nFreqs, nTimes);

%         for layi = 1:nLays
%             Z = squeeze(networkRDMs(layi,:,:)); 
%             Zall(:, layi) = vectorizeRDM(Z);
%         end

        CM = squeeze(load_CATMODEL_activ(cfg_contrasts.oneListIds));
        CM = vectorizeRDM(CM); 

        for layi = 1:nLays
            M = squeeze(networkRDMs(layi,:,:)); 
            M = vectorizeRDM(M);
%             zAT = Zall; 
%             zAT(:, layi) = []; 
            
            parfor freqi = 1:nFreqs
                for timei = 1:nTimes
                    rdm = squeeze(neuralRDMs(:, :, freqi, timei));
                    rdm = vectorizeRDM(rdm);
                    %allTEst = partialcorr(rdm, M, zAT, 'type', 's');
                    allTEst = partialcorr(rdm, M, CM, 'type', 's');
                    all_r_Times(layi,freqi,timei) = allTEst;  
                end
            end
        end
        
    else
        

    end
    
    
    
    
    
else % if only 3 dimensions (frequencies used as features)
    
    
end
    


end