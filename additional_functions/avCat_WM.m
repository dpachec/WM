function[wCatSim] = avCat_WM(rdm)


    rdm(logical(eye(size(rdm, 1)))) = nan; 

    allCM= [mean(rdm(1:10, 1:10), 'all', 'omitnan') mean(rdm(1:10, 11:20), 'all', 'omitnan') mean(rdm(1:10, 21:30), 'all', 'omitnan') mean(rdm(1:10, 31:40), 'all', 'omitnan') mean(rdm(1:10, 41:50), 'all', 'omitnan') mean(rdm(1:10, 51:60), 'all', 'omitnan') ;...
        mean(rdm(11:20, 1:10), 'all', 'omitnan') mean(rdm(11:20, 11:20), 'all', 'omitnan') mean(rdm(11:20, 21:30), 'all', 'omitnan') mean(rdm(11:20, 31:40), 'all', 'omitnan') mean(rdm(11:20, 41:50), 'all', 'omitnan') mean(rdm(11:20, 51:60), 'all', 'omitnan') ;...
        mean(rdm(21:30, 1:10), 'all', 'omitnan') mean(rdm(21:30, 11:20), 'all', 'omitnan') mean(rdm(21:30, 21:30), 'all', 'omitnan') mean(rdm(21:30, 31:40), 'all', 'omitnan') mean(rdm(21:30, 41:50), 'all', 'omitnan') mean(rdm(21:30, 51:60), 'all', 'omitnan') ;...
        mean(rdm(31:40, 1:10), 'all', 'omitnan') mean(rdm(31:40, 11:20), 'all', 'omitnan') mean(rdm(31:40, 21:30), 'all', 'omitnan') mean(rdm(31:40, 31:40), 'all', 'omitnan') mean(rdm(31:40, 41:50), 'all', 'omitnan') mean(rdm(31:40, 51:60), 'all', 'omitnan') ;...
        mean(rdm(41:50, 1:10), 'all', 'omitnan') mean(rdm(41:50, 11:20), 'all', 'omitnan') mean(rdm(41:50, 21:30), 'all', 'omitnan') mean(rdm(41:50, 31:40), 'all', 'omitnan') mean(rdm(41:50, 41:50), 'all', 'omitnan') mean(rdm(41:50, 51:60), 'all', 'omitnan') ;...
        mean(rdm(51:60, 1:10), 'all', 'omitnan') mean(rdm(51:60, 11:20), 'all', 'omitnan') mean(rdm(51:60, 21:30), 'all', 'omitnan') mean(rdm(51:60, 31:40), 'all', 'omitnan') mean(rdm(51:60, 41:50), 'all', 'omitnan') mean(rdm(51:60, 51:60), 'all', 'omitnan') ;...
        ];


    wCatSim = allCM (logical(eye(size(allCM)))); 





end