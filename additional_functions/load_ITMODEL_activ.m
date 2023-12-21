
function [ACT] = load_ITMODEL_activ(oneListIDs)

if ~strcmp(oneListIDs{1}(1), '7')
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

else

    ids = cellfun(@(x) strsplit(x), oneListIDs, 'UniformOutput', false);
    cues = double(string(cellfun(@(x) x(2), ids, 'UniformOutput', false)));
    
    
    clear idsF
    for cuei = 1:length(cues)
        idsF(cuei,:) = double(string(ids{cuei}(12+cues(cuei))));
    end

    M = zeros (length(idsF));
    
    for i = 1:length(M)
        for j = 1:length(M)
            x = idsF(i); 
            y = idsF(j); 
            if x == y 
                M(i, j) = 1;
            end
        end 
    end    
    ACT(1,:,:) = M; 


end
      
end
