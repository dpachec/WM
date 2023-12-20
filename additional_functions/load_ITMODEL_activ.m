
function [ACT] = load_ITMODEL_activ(oneListIDs)

    ids = cellfun(@(x) strsplit(string(x)), oneListIDs, 'UniformOutput', false);
    ids0 = double(string(cellfun(@(x) x(3), ids, 'UniformOutput', false)));

    M = zeros (length(ids0));
    
    for i = 1:length(M)
        for j = 1:length(M)
            x = ids0(i); 
            y = ids0(j); 
            if x == y 
                M(i, j) = 1;
            end
        end 
    end    
    ACT(1,:,:) = M; 
      
end
