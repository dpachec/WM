function[vectorizedRDM] = vectorizeRDM(rdm)
    rows = size(rdm, 1);
    vectorizedRDM = rdm(tril(true(rows), -1));
end


