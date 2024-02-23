function[vectorizedRDM] = vectorizeRDM(rdm)

    if ndims(rdm)== 2
        rows = size(rdm, 1);
        vectorizedRDM = rdm(tril(true(rows), -1));
    elseif ndims(rdm)== 3
        rows = size(rdm, 2);
        rdmTemplate = zeros(rows); 
        rdmTemplate(tril(true(rows), -1)) = 1;
        size1 = size(rdm, 1); 
        rdmTemplateMovie = repmat(rdmTemplate, 1, 1, size1); 
        rdmTemplateMovie = permute(rdmTemplateMovie, [3 1 2]); 
        rdm(rdmTemplateMovie==0) = []; 
        vectorizedRDM = reshape(rdm, size1, []); 

    elseif ndims(rdm)== 4
        rows = size(rdm, 1);
        rdm = reshape(rdm, rows, rows, []);
        rdmTemplate = zeros(rows); 
        rdmTemplate(tril(true(rows), -1)) = 1; 
        size3 = size(rdm, 3); 
        rdmTemplateMovie = repmat(rdmTemplate, 1, 1, size3); 
        rdm(rdmTemplateMovie==0) = []; 
        vectorizedRDM = reshape(rdm, [], size3); 
    end
end


