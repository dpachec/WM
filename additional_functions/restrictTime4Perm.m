
function [neuralRDMs] = restrictTime4Perm(cfg, neuralRDMs)

% % % restrict time for permutation data
    if strcmp(cfg.period(1), 'M')
        if ndims(neuralRDMs) == 4
            neuralRDMs = neuralRDMs(:,:,:,6:40); %frequency-resolved
        else
            neuralRDMs = neuralRDMs(:,:,6:40); %band analysis
        end
    else
        if ndims(neuralRDMs) == 4
            neuralRDMs = neuralRDMs(:,:,:,6:15); %frequency-resolved
        else
            neuralRDMs = neuralRDMs(:,:,6:15); %band analysis
        end
    end

end