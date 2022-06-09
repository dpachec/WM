function [data] = norm_trials(data)
    mT = mean(data, 1);
    stdT = std(data, [], 1);
    data = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);
end