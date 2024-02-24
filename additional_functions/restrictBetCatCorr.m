function [networkRDMs] = restrictBetCatCorr(cfg_contrasts, networkRDMs)

ITM                         = squeeze(load_ITMODEL_activ(cfg_contrasts.oneListIds));
ITMR                        = repmat(ITM, 1, 1, size(networkRDMs,1));
ITMR                        = permute(ITMR, [3 1 2]);
networkRDMs(isnan(ITMR))    = nan; 


end