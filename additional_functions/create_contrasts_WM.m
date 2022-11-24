function [out_contrasts] = create_contrasts_WM (cfg_contrasts)
 
if isfield(cfg_contrasts, 'oneListIds_enc')
    oneListIds              =       [cfg_contrasts.oneListIds_enc ; cfg_contrasts.oneListIds_maint];
    oneListPow              =       cat(1, cfg_contrasts.oneListPow_enc, cfg_contrasts.oneListPow_maint);
else
    oneListIds              =       cfg_contrasts.oneListIds_c;
    oneListPow              =       cfg_contrasts.oneListPow;
end
chanNames               =       cfg_contrasts.chanNames;
contr2save              =       cfg_contrasts.contr2save;
n2s                     =       cfg_contrasts.n2s;
loadSurr                =       cfg_contrasts.loadSurr; 
batch_bin               =       cfg_contrasts.batch_bin;

disp ('>>>>> creating contrasts');
 
allContrasts = []; allContrastIds = [];

for ci = 1:length(contr2save)
    eval(['count' contr2save{ci} '= 1;'])
    eval([contr2save{ci} '= [];'])
end

cr2c = cellfun(@(x) strsplit(x,'_'), contr2save, 'un', 0);
c2c = unique(cellfun(@(x) x{2}, cr2c, 'un', 0));

for i = 1:length(oneListIds)
    
    evei = strsplit(oneListIds{i});
    posEnc = evei(1);
    

    

% % % % % % 3 classics
    
% % % % % % Encoding - Encoding
if ~isempty(intersect(c2c, 'EE'))
    if strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*') %* is for the average data
        for j = 1:length(oneListIds) % repetitions are needed so should not start at i (we are only saving half of the matrix)
           evej = strsplit(oneListIds{j});
           if strcmp(evej(1), '1') | strcmp(evej(1), '3') | strcmp(evej(1), '5') | strcmp(evej(1), '*') 
            trli = string(evei(12)); trlj = string(evej(12));
            trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
            trlij = [trli trlj];
            
            if length(trlij) == length(unique(trlij)) % all averaged items were presented in different trials
                  it_id_Enc = evei(3); it_id_Ret = evej(3); 
                  cat_id_Enc = it_id_Enc{1}(1); cat_id_Ret = it_id_Ret{1}(1);
                  if strcmp(it_id_Enc,it_id_Ret) 
                     %disp (['SISC_EE> ' oneListIds{i} '//' oneListIds{j}]);   
                     if exist('countSISC_EE')
                        new_sisc_ee{countSISC_EE} = [i, j];
                        countSISC_EE = countSISC_EE+1;
                     end
                  elseif ~strcmp(it_id_Enc,it_id_Ret) & strcmp(cat_id_Enc,cat_id_Ret)
                     %disp (['DISC_EE > ' trlij]);    
                     if exist('countDISC_EE')
                         new_disc_ee{countDISC_EE} = [i, j];
                         countDISC_EE = countDISC_EE+1;
                     end
                  elseif ~strcmp(cat_id_Enc,cat_id_Ret)
                     %disp (['DIDC_EE > ' oneListIds{i} '//' oneListIds{j}]);     
                     if exist('countDIDC_EE')
                        new_didc_ee{countDIDC_EE} = [i, j];
                        countDIDC_EE = countDIDC_EE+1;
                     end
                   end
               end
           %end 
        end
        end
    end
    end
    
    
    
   % % %  %ENCODING - MAINTENANCE 2
   if ~isempty(intersect(c2c, 'EM2'))
    if strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*')
        for j = 1:length(oneListIds)
            evej = strsplit(oneListIds{j});
            if ( strcmp(evej(1), '7') | strcmp(evej(1), '7*') | strcmp(evej(1), '-7') ) & ~strcmp(evej(2), '4')
               trli = string(evei(12)); trlj = string(evej(12));
               trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
               trlij = [trli trlj];
            
               if length(trlij) == length(unique(trlij)) % there are no repetitions
                  idEnc = evei(3);
                  cueRet = str2double(evej{2});
                  idx = 12 + cueRet;
                  idRet = evej(idx);
                  oneIdEnc = idEnc{1}(1); oneIdRet = idRet{1}(1);
                  if strcmp(idEnc,idRet) 
                     %disp (['SISC_EM2> ' oneListIds{i} '//' oneListIds{j}]);  
                     if exist('countSISC_EM2')
                        new_sisc_em2{countSISC_EM2} = [i, j];
                        countSISC_EM2 = countSISC_EM2+1;
                     end
                  end
                   if ~strcmp(idEnc,idRet) & strcmp(oneIdEnc, oneIdRet)
                        %disp (['DISC_EM2 > ' trlij]);    
                        if exist('countDISC_EM2')
                            new_disc_em2{countDISC_EM2} = [i, j];
                            countDISC_EM2 = countDISC_EM2+1;
                        end
                   end
                   if ~strcmp(oneIdEnc, oneIdRet) 
                        %disp (['DIDC_EM2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        if exist('countDIDC_EM2')
                            new_didc_em2{countDIDC_EM2} = [i, j];
                            countDIDC_EM2 = countDIDC_EM2+1;
                        end
                   end
              end 
           end
        end
   % end
    end
   end
    
    
    
    % % % % % 1) M2M2 for 123 Cued items (only category)
    if ~isempty(intersect(c2c, 'M2M2'))
      if ( strcmp(evei(1), '7') | strcmp(evei(1), '7*') | strcmp(evei(1), '-7') | strcmp(evei(1), '9'))  & ~strcmp(evei(2), '4')
        cueRetI = str2double(evei{2});
        idx = 12 + cueRetI;
        idRet_I = double(string(evei(idx)));
        cuedCI = floor(idRet_I / 100);
        
        for j = i:length(oneListIds)
           evej = strsplit(oneListIds{j});
            trli = string(evei(12)); trlj = string(evej(12));
           trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
           trlij = [trli trlj];
            
           if length(trlij) == length(unique(trlij)) % all averaged instances are from different trials
            if ( strcmp(evej(1), '7') | strcmp(evej(1), '7*') | strcmp(evej(1), '-7') | strcmp(evei(1), '9')  )  & ~strcmp(evei(2), '4')
                  cueRetJ = str2double(evej{2});idx = 12 + cueRetJ;
                  idRet_J = double(string(evej(idx)));
                  cuedCJ = floor(idRet_J / 100);
                  if idRet_I == idRet_J
                      if exist('countSISC_M2M2')
                        %disp (['SISC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_sisc_m2m2{countSISC_M2M2} = [i, j];
                        countSISC_M2M2 = countSISC_M2M2+1; 
                      end 
                  end
                  if cuedCI == cuedCJ
                    if exist('countDISC_M2M2')
                        %disp (['DISC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_disc_m2m2{countDISC_M2M2} = [i, j];
                        countDISC_M2M2 = countDISC_M2M2+1;
                    end
                  else
                    if exist('countDIDC_M2M2')
                        %disp (['DIDC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_didc_m2m2{countDIDC_M2M2} = [i, j];
                        countDIDC_M2M2 = countDIDC_M2M2+1;
                    end
                end
              end 
           end
        end
      end
    end
    
    
    
    
    
    
% % % % % % % % % % M2M2 123 cued trials V1
    if ~isempty(intersect(c2c, 'M2123V1'))
      if strcmp(evei(1), '7') & ~strcmp(evei(2), '4')
        cuedIt_I = double(string(char(evei(12 + str2double(evei{2})))));
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, cuedIt_I);
        cuedCat_I = floor(cuedIt_I / 100); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') & ~strcmp(evej(2), '4')
                    cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                    cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                    nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if cuedCat_I == cuedCat_J & length(unique(nonCuedCatIJ)) == 4 & length(unique(allItIJ)) == 6 
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                        if exist('countDISC_M2123V1')
                            new_disc_m2123v1{countDISC_M2123V1} = [i, j];
                            countDISC_M2123V1 = countDISC_M2123V1+1;
                        end
                    end
                    if length(unique(allCatIJ)) == 6
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2123V1')
                            new_didc_m2123v1{countDIDC_M2123V1} = [i, j];
                            countDIDC_M2123V1 = countDIDC_M2123V1+1;
                        end
                    end
               end
           end
        end
      end
    end

      
% % % % % % % % % % M2M2 123 cued trials V2
if ~isempty(intersect(c2c, 'M2123V2'))
      if strcmp(evei(1), '7') & ~strcmp(evei(2), '4')
        cuedIt_I = double(string(char(evei(12 + str2double(evei{2})))));
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)]))));
        nonCuedIt_I = setdiff(allIt_I, cuedIt_I);
        cuedCat_I = floor(cuedIt_I / 100); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           allItIJ = [allIt_I ; allIt_J]; 
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') & ~strcmp(evej(2), '4')
                    cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)]))));
                    nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                    cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                    nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                    allItIJ = [allIt_I ; allIt_J]; 
                    if cuedCat_I == cuedCat_J & length(unique(allItIJ)) == 6 
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                       if exist('countDISC_M2123V2')
                                new_disc_m2123v2{countDISC_M2123V2} = [i, j];
                                countDISC_M2123V2 = countDISC_M2123V2+1;
                        end
                    end
                    if cuedCat_I ~= cuedCat_J & length(unique(allItIJ)) == 6; 
                        %disp (['All different categories > ' string(trlijIt)]);     
                       if exist('countDIDC_M2123V2')
                                new_didc_m2123v2{countDIDC_M2123V2} = [i, j];
                                countDIDC_M2123V2 = countDIDC_M2123V2+1;
                        end
                    end
               end
           end
        end
      end
end
    
          
% % % % % % % % % % M2M2 123 cued and non-cued trials V1
if ~isempty(intersect(c2c, 'M2123CNCV1'))
      if strcmp(evei(1), '7') & ~strcmp(evei(2), '4')
        cuedIt_I = double(string(char(evei(12 + str2double(evei{2})))));
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, cuedIt_I);
        cuedCat_I = floor(cuedIt_I / 100); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') & ~strcmp(evej(2), '4')
                    cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                    cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                    nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(allCatIJ)) == 3 & length(unique(allItIJ)) == 6 & cuedCat_I ~= cuedCat_J
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                        if exist('countDISC_M2123CNCV1')
                            new_disc_m2123cncv1{countDISC_M2123CNCV1} = [i, j];
                            countDISC_M2123CNCV1 = countDISC_M2123CNCV1+1;
                        end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2123CNCV1')
                            new_didc_m2123cncv1{countDIDC_M2123CNCV1} = [i, j];
                            countDIDC_M2123CNCV1 = countDIDC_M2123CNCV1+1;
                        end
                    end
               end
           end
        end
      end
end
   
      % % % % % % % % % % M2M2 123 cued and non-cued trials V2
  if ~isempty(intersect(c2c, 'M2123CNCV2'))
      if strcmp(evei(1), '7') & ~strcmp(evei(2), '4')
        cuedIt_I = double(string(char(evei(12 + str2double(evei{2})))));
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, cuedIt_I);
        cuedCat_I = floor(cuedIt_I / 100); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') & ~strcmp(evej(2), '4')
                    cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                    cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                    nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(allCatIJ)) == 3 & length(unique(allItIJ)) == 6 
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                        if exist('countDISC_M2123CNCV2')
                            new_disc_m2123cncv2{countDISC_M2123CNCV2} = [i, j];
                            countDISC_M2123CNCV2 = countDISC_M2123CNCV2+1;
                        end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2123CNCV2')
                            new_didc_m2123cncv2{countDIDC_M2123CNCV2} = [i, j];
                            countDIDC_M2123CNCV2 = countDIDC_M2123CNCV2+1;
                        end
                    end
               end
           end
        end
      end
  end
      
      
      % % % % % % % % % % M2M2 123 non-cued 
  if ~isempty(intersect(c2c, 'M2123NC'))

      if (strcmp(evei(1), '7') | strcmp(evei(1), '7*')) & ~strcmp(evei(2), '4')
        cuedIt_I = double(string(char(evei(12 + str2double(evei{2})))));
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, cuedIt_I);
        cuedCat_I = floor(cuedIt_I / 100); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
            trli = string(evei(12)); trlj = string(evej(12));
           trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
           trlij = [trli trlj];
                          
           if length(trlij) == length(unique(trlij)) % all averaged trials are from different trials
               if (strcmp(evej(1), '7') | strcmp(evej(1), '7*')) & ~strcmp(evej(2), '4')
                    cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                    cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                    nonCuedItIJ = [nonCuedIt_I ; nonCuedIt_J]; 
                    nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(nonCuedCatIJ)) == 2 & length(unique(nonCuedItIJ)) == 4 & cuedCat_I ~= cuedCat_J
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                        if exist('countDISC_M2123NC')
                            new_disc_m2123nc{countDISC_M2123NC} = [i, j];
                            countDISC_M2123NC = countDISC_M2123NC+1;
                        end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2123NC')
                            new_didc_m2123nc{countDIDC_M2123NC} = [i, j];
                            countDIDC_M2123NC = countDIDC_M2123NC+1;
                        end
                    end
               end
           end
        end
      end
  end
      
if ~isempty(intersect(c2c, 'EM2UV1'))
      if strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*')
        idEnc = double(string(evei{3})); idCEnc = floor(idEnc/100); 
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, idEnc); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
               if ( strcmp(evej(1), '7') | strcmp(evej(1), '7*') | strcmp(evej(1), '-7') ) & ~strcmp(evej(2), '4')
                   trli = string(evei(12)); trlj = string(evej(12));
                   trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
                   trlij = [trli trlj];
                                  
                   if length(trlij) == length(unique(trlij)) % all averaged trials are from different trials
                        cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                        allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                        nonCuedIt_J = setdiff(allIt_J, cuedIt_J);
                        cuedCat_J = floor(cuedIt_J / 100); nonCuedCat_J = floor(nonCuedIt_J / 100); 
                        nonCuedItIJ = [nonCuedIt_I ; nonCuedIt_J]; 
                        nonCuedCatIJ = [nonCuedCat_I ; nonCuedCat_J]; 
                        allItIJ = [allIt_I ; allIt_J]'; 
                        allCatIJ = [allCat_I ; allCat_J]'; 
                        if length(unique(nonCuedCatIJ)) == 2 & length(unique(nonCuedItIJ)) == 4 & idCEnc ~= cuedCat_J
                            %disp (['EM2U DISC V1 > '  string(allItIJ)]);     % oneListIds{i} '//' oneListIds{j}]);     
                            if exist('countDISC_EM2UV1')
                                new_disc_em2uv1{countDISC_EM2UV1} = [i, j];
                                countDISC_EM2UV1 = countDISC_EM2UV1+1;
                            end
                        end
                        if length(unique(allCatIJ)) == 6 
                            %disp (['EM2U DIDC V1 > '  string(allItIJ)]);    
                            if exist('countDIDC_EM2UV1')
                                new_didc_em2uv1{countDIDC_EM2UV1} = [i, j];
                                countDIDC_EM2UV1 = countDIDC_EM2UV1+1;
                            end
                        end
               end
           end
        end
      end
  end


    if ~isempty(intersect(c2c, 'EM2UV2')) %Encoding maintenance uncued
      if strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*')
        idEnc = double(string(evei{3})); idCEnc = floor(idEnc/100); 
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        nonCuedIt_I = setdiff(allIt_I, idEnc); nonCuedCat_I = floor(nonCuedIt_I / 100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
               if ( strcmp(evej(1), '7') | strcmp(evej(1), '7*') | strcmp(evej(1), '-7') ) & ~strcmp(evej(2), '4')
                   trli = string(evei(12)); trlj = string(evej(12));
                   trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
                   trlij = [trli trlj];
                   
                   
                   if length(trlij) == length(unique(trlij)) % all averaged trials are from different trials
                       cuedIt_J = double(string(char(evej(12 + str2double(evej{2})))));
                       allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                       nonCuedIt_J = setdiff(allIt_J, cuedIt_J);cuedCat_J = floor(cuedIt_J / 100); 
                       nonCuedCat_J = floor(nonCuedIt_J / 100); 
                       allItIJ = [allIt_I ; allIt_J]; 
                       allCatIJ = [allCat_I ; allCat_J]; 
    
                        cncI = [nonCuedIt_I(:) ; nonCuedIt_J(:)]';
                        cncC = [nonCuedCat_I(:) ; nonCuedCat_J(:)]';
                        
                        if length(unique(cncI)) == 3  & idEnc ~= cuedIt_J% 
                            %disp (['EM2UV2 SISC > '  string(cncI) num2str(idCEnc) '-' num2str(cuedCat_J)]);     % oneListIds{i} '//' oneListIds{j}]);     
                            if exist('countSISC_EM2UV2')
                                new_sisc_em2uv2{countSISC_EM2UV2} = [i, j];
                                countSISC_EM2UV2 = countSISC_EM2UV2+1;
                            end
                        end
                        if length(unique(cncC)) == 3 & length(unique(cncI)) == 4 & idCEnc ~= cuedCat_J  & isempty(intersect(cuedCat_J, cncC))
                            %disp (['EM2UV2 DISC > '  string(cncI) num2str(idCEnc) '-' num2str(cuedCat_J)]);     % oneListIds{i} '//' oneListIds{j}]);     
                            if exist('countDISC_EM2UV2')
                                new_disc_em2uv2{countDISC_EM2UV2} = [i, j];
                                countDISC_EM2UV2 = countDISC_EM2UV2+1;
                            end
                        end
                        if length(unique(cncC)) == 4 & idCEnc ~= cuedCat_J & isempty(intersect(cuedCat_J, cncC))
                            %disp (['EM2UV2 DIDC > '  string(cncI) num2str(idCEnc) '-' num2str(cuedCat_J)]);     % oneListIds{i} '//' oneListIds{j}]);     
                            if exist('countDIDC_EM2UV2')
                                new_didc_em2uv2{countDIDC_EM2UV2} = [i, j];
                                countDIDC_EM2UV2 = countDIDC_EM2UV2+1;
                            end
                        end
                   end
               end
            end
      end
  end  
       % % % % % % % % % % M2A All trials 
   if ~isempty(intersect(c2c, 'M2A'))
      if ( strcmp(evei(1), '7') | strcmp(evei(1), '7*') | strcmp(evei(1), '-7') ) & strcmp(evei(2), '4')
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           trli = string(evei(12)); trlj = string(evej(12));
           trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
           trlij = [trli trlj];
            
            if length(trlij) == length(unique(trlij)) % all averaged items were presented in different trials
               if ( strcmp(evej(1), '7') | strcmp(evej(1), '7*') | strcmp(evej(1), '-7') ) & strcmp(evej(2), '4')
                   %disp('hola')
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(allCatIJ)) == 3 & length(unique(allItIJ)) == 6
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                       if exist('countDISC_M2A')
                            new_disc_m2a{countDISC_M2A} = [i, j];
                            countDISC_M2A = countDISC_M2A+1;
                       end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2A')
                            new_didc_m2a{countDIDC_M2A} = [i, j];
                            countDIDC_M2A = countDIDC_M2A+1;
                        end
                    end
               end
           end
        end
      end
   end
      
   % % % % % % % % % % M2M2 all + 123 trials (all trials + 123 trials (cued and non-cued))
   if ~isempty(intersect(c2c, 'M2A123'))
      if strcmp(evei(1), '7')
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') 
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(allCatIJ)) == 3 & length(unique(allItIJ)) == 6
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                       if exist('countDISC_M2A123')
                                new_disc_m2a123{countDISC_M2A123} = [i, j];
                                countDISC_M2A123 = countDISC_M2A123+1;
                       end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M2A123')
                            new_didc_m2a123{countDIDC_M2A123} = [i, j];
                            countDIDC_M2A123 = countDIDC_M2A123+1;
                        end
                    end
               end
           end
        end
      end       
   end
    
    
    
    
    
    % % % % %  %ENCODING CUED CORRECT-INCORRECT
    if ~isempty(intersect(c2c, 'CI')) |  ~isempty(intersect(c2c, 'CC'))  | ~isempty(intersect(c2c, 'CCNI'))  |  ~isempty(intersect(c2c, 'I')) 
        if ( strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*') ) & ~strcmp(evei(2), '4') & ...
            strcmp( evei(3) ,  evei (12 + str2double(evei{2})) )
            for j = 1:length(oneListIds)
               evej = strsplit(oneListIds{j});
               if strcmp(evej(1), '1') | strcmp(evej(1), '3') | strcmp(evej(1), '5') | strcmp(evej(1), '*') & ~strcmp(evej(2), '4') & ...
                strcmp( evej(3) ,  evej (12 + str2double(evej{2})) )
                    trli = string(evei(12)); trlj = string(evej(12));
                    trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
                    trli = double(trli); trlj = double(trlj);
                    trlij = [trli trlj];
                    if length(trlij) == length(unique(trlij)) % different trials
                        %if strcmp(idEnc,idRet) 
                          idEnc = evei(3); idRet = evej(3); 
                          oneIdEnc = idEnc{1}(1); oneIdRet = idRet{1}(1);
                          if strcmp(evei(19), '1') & strcmp(evej(19), '1') 
                             %disp (['CI > ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countCI')
                                new_ci{countCI} = [i, j];
                                countCI = countCI+1;
                             end
                          end
                          if strcmp(evei(20), '1')  & strcmp(evej(20), '1') 
                             %disp (['CC > ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countCC')
                                new_cc{countCC} = [i, j];
                                countCC = countCC+1;
                             end
                          end
                          if strcmp(evei(19), '0') & strcmp(evei(20), '1') &  strcmp(evej(19), '0') & strcmp(evej(20), '1') 
                             %disp (['CCNI > ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countCCNI')
                                new_ccni{countCCNI} = [i, j];
                                countCCNI = countCCNI+1;
                             end
                          end
                          if strcmp(evei(19), '0') & strcmp(evei(20), '0') & strcmp(evej(19), '0') & strcmp(evej(20), '0') 
                             %disp (['I > ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countI')
                                new_i{countI} = [i, j];
                                countI = countI+1;
                             end
                          end
                        %end
                     end
                 end 
            end
        end
    end
    
    
        
% % % % %   %M2M2 for 123 non cued items
if ~isempty(intersect(c2c, '123M2A'))
      if strcmp(evei(1), '7') & ~strcmp(evei(2), '4')
        cueRetI = str2double(evei{2});
        idx = 12 + cueRetI;
        idRet_I = evei(idx);
        oneIdM2_I = idRet_I{1}(1);
        allItEnc = [evei{13} '_' evei{14} '_' evei{15}];
        cuedCI = double(string(char(oneIdM2_I)));
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') & ~strcmp(evej(2), '4')
                  cueRetJ = str2double(evej{2});
                  idx = 12 + cueRetJ;
                  idRet_J = evej(idx);
                  oneIdM2_J = idRet_J{1}(1); 
                  allItMaint = [evej{13} '_' evej{14} '_' evej{15}];
                  cuedCJ = double(string(char(oneIdM2_J)));
                      
                  trliIt = string(allItEnc); trljIt = string(allItMaint);
                  trliIt = strsplit(trliIt, '_'); trljIt = strsplit(trljIt, '_');
                  trliIt = double(trliIt); trljIt = double(trljIt);
                  trlijIt = [trliIt trljIt];

                  trlijCat = floor(trlijIt./100);
                  
                  if cuedCI ~= cuedCJ
                      if length(unique(trlijCat)) <= 3 & length(unique(trlijIt)) == 6 % all diff items and category repetitions
                        %disp (['Diff item same category > ' string(trlijIt)]);     
                        if exist('countDISC_123M2A')
                            new_disc_123m2a{countDISC_123M2A} = [i, j];
                            countDISC_123M2A = countDISC_123M2A+1;
                        end
                      end
                      if length(unique(trlijCat)) == 6 % all different categories
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_123M2A')
                            new_didc_123m2a{countDIDC_123M2A} = [i, j];
                            countDIDC_123M2A = countDIDC_123M2A+1;
                        end
                      end
                  end
              end 
           end
        end
      end
 end


% %  % %  M2 All and 123 combined
if ~isempty(intersect(c2c, '123M2A'))
     if strcmp(evei(1), '7') 
        allItEnc = [evei{13} '_' evei{14} '_' evei{15}];
        cueRetI = str2double(evei{2});
        idx = 12 + cueRetI;
        idRet_I = evei(idx);
        oneIdM2_I = idRet_I{1}(1);
        cuedCI = double(string(char(oneIdM2_I)));
        for j = i:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '7') 
                  cueRetJ = str2double(evej{2});
                  idx = 12 + cueRetJ;
                  idRet_J = evej(idx);
                  oneIdM2_J = idRet_J{1}(1); 
                  cuedCJ = double(string(char(oneIdM2_J)));

                  allItMaint = [evej{13} '_' evej{14} '_' evej{15}];

                  trliIt = string(allItEnc); trljIt = string(allItMaint);
                  trliIt = strsplit(trliIt, '_'); trljIt = strsplit(trljIt, '_');
                  trliIt = double(trliIt); trljIt = double(trljIt);
                  trlijIt = [trliIt trljIt];

                  trlijCat = floor(trlijIt./100);

                  %if cuedCI == cuedCJ
                      if length(unique(trlijCat)) == 3 & length(unique(trlijIt)) == 6 % all diff items and category repetitions
                        %disp (['DISC M2A > ' string(trlijIt)]);     
                        if exist('countDISC_123M2A')
                            new_disc_m2a{countDISC_M2A} = [i, j];
                            countDISC_M2A = countDISC_M2A+1;
                        end
                      end
                  %end
                      if length(unique(trlijCat)) == 6 % all different categories
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_123M2A')
                            new_didc_m2a{countDIDC_M2A} = [i, j];
                            countDIDC_M2A = countDIDC_M2A+1;
                        end
                      end
                end 
           end
        end
    end
end

     
      

% % % % % % Encoding - Encoding CORRECT
if ~isempty(intersect(c2c, 'EEC'))
    if (strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5') ) & ... %* is for the average data
        strcmp(evei{8}, '1') & ( strcmp(evei{1}, evei{2})  ) % correct at the item level
        for j = 1:length(oneListIds) % repetitions are needed so should not start at i (we are only saving half of the matrix)
           evej = strsplit(oneListIds{j});
           if (strcmp(evej(1), '1') | strcmp(evej(1), '3') | strcmp(evej(1), '5') | strcmp(evej(1), '*') ) & ...
               strcmp(evej{8}, '1') & ( strcmp(evej{1}, evej{2})  )
                trli = string(evei(12)); trlj = string(evej(12));
                trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
                trlij = [trli trlj];
            
            if length(trlij) == length(unique(trlij)) % all averaged items were presented in different trials
                  it_id_Enc = evei(3); it_id_Ret = evej(3); 
                  cat_id_Enc = it_id_Enc{1}(1); cat_id_Ret = it_id_Ret{1}(1);
                  if strcmp(it_id_Enc,it_id_Ret) 
                     %disp (['SISC_EEC> ' oneListIds{i} '//' oneListIds{j}]);   
                     if exist('countSISC_EEC')
                        new_sisc_eec{countSISC_EEC} = [i, j];
                        countSISC_EEC = countSISC_EEC+1;
                     end
                  elseif ~strcmp(it_id_Enc,it_id_Ret) & strcmp(cat_id_Enc,cat_id_Ret)
                     %disp (['DISC_EEC > ' trlij]);    
                     if exist('countDISC_EEC')
                         new_disc_eec{countDISC_EEC} = [i, j];
                         countDISC_EEC = countDISC_EEC+1;
                     end
                  elseif ~strcmp(cat_id_Enc,cat_id_Ret)
                     %disp (['DIDC_EE > ' oneListIds{i} '//' oneListIds{j}]);     
                     if exist('countDIDC_EEC')
                        new_didc_eec{countDIDC_EEC} = [i, j];
                        countDIDC_EEC = countDIDC_EEC+1;
                     end
                   end
               end
           %end 
        end
        end
    end
    end

      
          
end

   


% restrict number of different items
%length(contr2save)
for ci = 1:length(contr2save)
    
    % % select only a subset of combinations
     eval(['if exist(''new_' lower(contr2save{ci}) ''') ' newline ...
           '   disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
           '   if ~loadSurr' newline ... % & isempty(allPermIds(' num2str(ci) '))'
           '     if length(new_' lower(contr2save{ci}) ') > n2s ' newline ...
           '        id2ns' lower(contr2save{ci}) ' = randsample(length(new_' lower(contr2save{ci}) '), n2s); ' newline ...
           '     else' newline ...
           '        id2ns' lower(contr2save{ci}) ' = 1:length(new_' lower(contr2save{ci}) '); ' newline ...
           '     end' newline ... 
           '     allPermIds_N{' num2str(ci) '} = id2ns' lower(contr2save{ci}) ';' newline ...
           '     new_' lower(contr2save{ci}) ' =  new_' lower(contr2save{ci}) '(id2ns' lower(contr2save{ci}) '); ' newline ...
           '   else' newline ...
           '        new_' lower(contr2save{ci}) ' =  new_' lower(contr2save{ci}) '(allPermIds{ci}); ' newline ...
           '   end' newline ... 
           'else' newline ...
           '        %allPermIds_N{' num2str(ci) '} = [];' newline ...
           'end' ]);

end


for ci = 1:length(contr2save)
       
    eval(['if exist(''new_' lower(contr2save{ci}) ''')    new_' lower(contr2save{ci}) ' =  new_' lower(contr2save{ci}) ''' ; '  ...
          'new_' lower(contr2save{ci}) ' = vertcat(new_' lower(contr2save{ci}) '{:}); end' ])   

end 


for ci = 1:length(contr2save)
eval(  ['if exist(''new_' lower(contr2save{ci}) ''') & any(strcmp(contr2save, ''' contr2save{ci} ''')) ...' newline ...
           'disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
           ' nTrials = length(new_' lower(contr2save{ci}) '); ' newline ...
           ' if nTrials > batch_bin ' newline ...
               ' nBatches = ceil(nTrials/batch_bin)' newline ... 
               contr2save{ci} '{nBatches} = [];' newline ...
               ' for batchi = 1:nBatches ' newline ... 
                   'if strcmp(''' contr2save{ci}(end-1:end) ''', ''EE'')' newline ... 
                       ' if batchi == 1 ' newline ... 
                            't = batchi:batchi*batch_bin;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,1), :, :, 1:25); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,2), :, :, 1:25); ...' newline ... 
                        ' elseif (batchi)*batch_bin < nTrials ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,1), :, :, 1:25); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,2), :, :, 1:25); ...' newline ... 
                            ' else ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 1), :, :, 1:25); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 2), :, :, 1:25); ...' newline ... 
                            'end' newline ... 
                 'else' newline ... 
                      ' if batchi == 1 ' newline ... 
                        't = batchi:batchi*batch_bin;' newline ... 
                        contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,1), :, :, :); ...' newline ... 
                        contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,2), :, :, :); ...' newline ... 
                    ' elseif (batchi)*batch_bin < nTrials ' newline ... 
                            't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,1), :, :, :); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,2), :, :, :); ...' newline ... 
                        ' else ' newline ... 
                            't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 1), :, :, :); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 2), :, :, :); ...' newline ... 
                        'end' newline ... 
               'end' newline...               
            'end' newline ... 
            'else' newline ...
                'if strcmp(''' contr2save{ci}(end-1:end) ''', ''EE'')' newline ... 
                    contr2save{ci} '{1}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,1), :, :, 1:25); ...' newline ... 
                    contr2save{ci} '{1}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,2), :, :, 1:25); ...' newline ... 
                'else' newline ... 
                    contr2save{ci} '{1}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,1), :, :, :); ...' newline ... 
                    contr2save{ci} '{1}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,2), :, :, :); ...' newline ... 
                'end' newline ...
            'end' newline ... 
            contr2save{ci} '_ID = [oneListIds(new_' lower(contr2save{ci}) '(:,1)) oneListIds(new_' lower(contr2save{ci}) '(:,2))];' newline ...
            'end'    ]      );
       
end


        


for ci = 1:length(contr2save)
    
    
allContrastIds = contr2save'; 
for i = 1:length(allContrastIds)
    if exist( allContrastIds{i} ) 
        allContrasts{i} = eval([allContrastIds{i}]);
    end
end
allContrasts = allContrasts';
 
for i = 1:length(allContrastIds)
    %[allContrastIds{i} '_ID']
    % exist( [allContrastIds{i} '_ID']) 
    if exist( [allContrastIds{i} '_ID']) 
        allIDs{i} = eval([allContrastIds{i} '_ID']);
    else
        allIDs{i} = [];
    end
    if exist( [allContrastIds{i} '_ids']) 
        all_ids{i} = eval([allContrastIds{i} '_ids']);
    else
        all_ids{i} = []; 
    end
    
end
allIDs = allIDs';
%all_ids = all_ids';
 
 
out_contrasts.allContrasts              = allContrasts;
out_contrasts.allContrastIds            = allContrastIds ;
out_contrasts.chanNames                 = chanNames; 
out_contrasts.allIDs                    = allIDs; 
 
 
end
 
disp ('>> conditions extracted');
 
 
%%
 
 
 

