function[rdm] = vectorizeRDM(rdm)
        rdm(rdm ==1) = 1000;rdm(rdm ==0) = 2000;
        rdm = tril(rdm, -1);rdm(rdm==0) =[];
        rdm(rdm==1000) = 1;rdm(rdm==2000) = 0;

end