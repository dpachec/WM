function [all_r_Times] = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, CoI)


if ndims(neuralRDMs) == 4 % if frequency resolved

    nLays  = size(networkRDMs, 1);
    nFreqs = size(neuralRDMs, 3); 
    nTimes = size(neuralRDMs, 4); 

    
    all_r_Times = zeros(nLays, nFreqs, nTimes);

%         for layi = 1:nLays
%             Z = squeeze(networkRDMs(layi,:,:)); 
%             Zall(:, layi) = vectorizeRDM(Z);
%         end

    if strcmp(CoI, 'C')
       CM = squeeze(load_CATMODEL_activ(cfg_contrasts.oneListIds));
       CM = vectorizeRDM(CM); 
       nids = isnan(CM); 
       CM(nids) = []; 

    elseif strcmp(CoI, 'I')
       CM = squeeze(load_ITMODEL_activ(cfg_contrasts.oneListIds));
       CM = vectorizeRDM(CM); 
       nids = isnan(CM); 
       CM(nids) = []; 
    end
    % 
    % for layi = 1:nLays
    %     M = squeeze(networkRDMs(layi,:,:)); 
    %     M = vectorizeRDM(M);
    %     M(nids) = []; 
    %     parfor freqi = 1:nFreqs
    %         for timei = 1:nTimes
    %             rdm = squeeze(neuralRDMs(:, :, freqi, timei));
    %             rdm = vectorizeRDM(rdm);
    %             rdm(nids) = []; 
    %             allTEst = partialcorr(rdm, M, CM, 'type', 's');
    %             all_r_Times(layi,freqi,timei) = allTEst;  
    %         end
    %     end
    % end
        
      
        % % same as code above in vectorized form
        % parfor layi = 1:nLays
        %     M =  squeeze(networkRDMs(layi,:,:)); 
        %     M = vectorizeRDM(M);
        % 
        %     rdm = neuralRDMs;
        %     rdm = vectorizeRDM(rdm);
        % 
        %     allTEst = partialcorr(rdm, M, CM, 'type', 's');
        %     all_r_Times(layi,:, :) = reshape(allTEst, nFreqs, nTimes);;  
        % 
        % end    

        % % % The layer loop can be avoided as well
        all_r_Times = zeros(nLays, nFreqs, nTimes);
        Ms =  vectorizeRDM(networkRDMs)';
        nanIds = isnan(Ms(:, 1)); 
        nanIds = nanIds|nids; 
        Ms(any(nanIds, 2), :) = [];
        rdm = neuralRDMs;
        rdm = vectorizeRDM(rdm);
        rdm(any(nanIds, 2), :) = [];
        allTEst = partialcorr(rdm, Ms, CM, 'type', 's');
        all_r_Times = reshape(allTEst', nLays, nFreqs, nTimes);
    
    
    
    
    
    
else % if only 3 dimensions (frequencies used as features)
    
    
end
    


end