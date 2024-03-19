function [trlinfo] = getPresentedItems(trlinfo)


for triali = 1:size(trlinfo, 1)
    randvec = trlinfo(triali, 16:21);
    colvec = randvec;
    colvec(randvec == 1) = 4;
    colvec(randvec == 2) = 5;
    colvec(randvec == 3) = 6;
    colvec(randvec == 4) = 12;
    colvec(randvec == 5) = 13;
    colvec(randvec == 6) = 14;
    
    
    n28 = trlinfo(triali,28);
    n29 = trlinfo(triali,29);
    n30 = trlinfo(triali,30);
    
    switch n28
        case 1
            cid = colvec(1);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        case 2
            cid = colvec(2);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        case 3
            cid = colvec(3);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        case 4
            cid = colvec(4);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        case 5
            cid = colvec(5);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        case 6
            cid = colvec(6);
            trlinfo(triali, 45) = trlinfo(triali,cid);
        otherwise
            trlinfo(triali, 45) = NaN;
    end
       
    switch n29
        case 1
            cid = colvec(1);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        case 2
            cid = colvec(2);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        case 3
            cid = colvec(3);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        case 4
            cid = colvec(4);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        case 5
            cid = colvec(5);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        case 6
            cid = colvec(6);
            trlinfo(triali, 46) = trlinfo(triali,cid);
        otherwise
           trlinfo(triali, 45) = NaN;
    end
       
    switch n30
        case 1
            cid = colvec(1);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        case 2
            cid = colvec(2);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        case 3
            cid = colvec(3);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        case 4
            cid = colvec(4);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        case 5
            cid = colvec(5);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        case 6
            cid = colvec(6);
            trlinfo(triali, 47) = trlinfo(triali,cid);
        otherwise
            trlinfo(triali, 45) = NaN;
    end
       
    
end



