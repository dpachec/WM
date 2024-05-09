
function [CCI] = compute_CCI_WM(mRDM, ind)

    %M = kron(eye(6), ones(10)); %equivalent to
    M = zeros(length(ind));
    [ids1 ids2] = unique(ind); 
    [ids1 ids3] = unique(ind, 'Last'); 
    for i = 1:6
        M(ids2(i):ids3(i), ids2(i):ids3(i)) = 1; 
    end

    mRDM(mRDM==0) = 2000; 
    mRDM = tril(mRDM, -1); mRDM(mRDM==0) = nan; 
    mRDM(mRDM==2000) = 0; 
    mWithin = mean(mRDM(M == 1), 'all', 'omitnan');
    mAcross = mean(mRDM(M == 0), 'all', 'omitnan');
    CCI = mWithin-mAcross;

end
