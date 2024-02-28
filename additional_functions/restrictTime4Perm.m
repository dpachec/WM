
function [neuralRDMs] = restrictTime4Perm(cfg, neuralRDMs)

% % % restrict time for permutation data
    if strcmp(cfg.period(1), 'M')
        if strcmp(cfg.period, 'M123')
            if ndims(neuralRDMs) == 4
                neuralRDMs = neuralRDMs(:,:,:,6:40); %frequency-resolved
            else
                neuralRDMs = neuralRDMs(:,:,6:40); %band analysis
            end
        else % restrict time differently for the M1 period
            if ndims(neuralRDMs) == 4
                neuralRDMs = neuralRDMs(:,:,:,6+8:6+8+20); %frequency-resolved
            else
                neuralRDMs = neuralRDMs(:,:,6+8:6+8+20); %band analysis
            end
        end
    else
        if ndims(neuralRDMs) == 4
            neuralRDMs = neuralRDMs(:,:,:,6:15); %frequency-resolved
        else
            neuralRDMs = neuralRDMs(:,:,6:15); %band analysis
        end
    end

end