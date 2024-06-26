
function [ACT nRep] = load_ITMODEL_activ(oneListIDs)


    ids = cellfun(@(x) strsplit(string(x)), oneListIDs, 'UniformOutput', false);
    ids0 = double(string(cellfun(@(x) x(3), ids, 'UniformOutput', false)));

    [a,b] = histc(ids0,unique(ids0));
    nRep = a(b);

    M = zeros (length(ids0));
    
    for i = 1:length(M)
        for j = 1:length(M)
            x = ids0(i); 
            y = ids0(j); 
            if x == y 
                M(i, j) = 1;
            elseif floor(x/100) ~= floor(y/100)
                M(i, j) = nan;
            end
        end 
    end    
    M(M==0) = 2; M = tril(M, -1); M(M==0) = 3; M(M==2) = 0; 
    ACT(1,:,:) = M; 
      
end
