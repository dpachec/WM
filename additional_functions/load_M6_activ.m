
function [ACT] = load_M6_activ(oneListIDs)


    ids = cellfun(@(x) strsplit(string(x)), oneListIDs, 'UniformOutput', false);
    ids0 = cellfun(@(x) x(3), ids, 'UniformOutput', false);

    M = zeros (length(ids0));
    
    for i = 1:length(M)
        for j = 1:length(M)
            x = char(ids0{i}); 
            y = char(ids0{j}); 
            if strcmp(x(1), y(1)) & ~strcmp(x, y)
                M(i, j) = x(1)-48;
            else
                M(i, j) = (x(1)-48)*-1;
            end
        end 
    end    

    %M(M==0) = 2; M = tril(M, -1); M(M==0) = 3; M(M==2) = 0; 
    M(logical(eye(size(M, 1)))) = nan; 
    ACT(1,:,:) = M; 
    
      
end
