function [data] = normalize_within_cond_WM(data)

    mT = mean(data,1, 'omitnan');
    stdT = std(data,[], 1, 'omitnan');
    data = bsxfun(@rdivide, data - mT, stdT);  
 
end
























                