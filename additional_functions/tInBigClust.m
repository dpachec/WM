function [tPerm] = tInBigClust(h, t); 
    clustinfo = bwconncomp(h);
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
