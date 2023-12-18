function [out_contrasts] = create_contrasts_WM (cfg_contrasts)
 
chanNames               =       cfg_contrasts.chanNames;
contr2save              =       cfg_contrasts.contr2save;
n2s                     =       cfg_contrasts.n2s;
loadSurr                =       cfg_contrasts.loadSurr; 
batch_bin               =       cfg_contrasts.batch_bin;
oneListIds              =       cfg_contrasts.oneListIds;
oneListPow              =       cfg_contrasts.oneListPow;

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
    
% % % % % % Encoding - Encoding
if ~isempty(intersect(c2c, 'EE'))
    if strcmp(posEnc, '1') | strcmp(posEnc, '3') | strcmp(posEnc, '5')| strcmp(posEnc, '*') %* is for the average data
        for j = 1:length(oneListIds) % repetitions are needed so should not start at i (we are only saving half of the matrix)
           evej = strsplit(oneListIds{j});
           if strcmp(evej(1), '1') | strcmp(evej(1), '3') | strcmp(evej(1), '5') | strcmp(evej(1), '*') 
            trli = string(evei(12)); trlj = string(evej(12));
            trli = strsplit(trli, '_'); trlj = strsplit(trlj, '_');
            trlij = [trli trlj];
            
            if length(trlij) == length(unique(trlij)) % all averaged items were presented in different trials. 
                    
                  it_id_Enc = evei(3); it_id_Ret = evej(3); 
                  cat_id_Enc = it_id_Enc{1}(1); cat_id_Ret = it_id_Ret{1}(1);

                  if exist('countALL_EE')
                    new_all_ee{countALL_EE} = [i, j];
                    countALL_EE = countALL_EE+1;
                  end
                  if strcmp(it_id_Enc,it_id_Ret) 
                     %disp(['SISC_EE> ' oneListIds{i} '//' oneListIds{j}]);   
                     if exist('countSISC_EE')
                        new_sisc_ee{countSISC_EE} = [i, j];
                        countSISC_EE = countSISC_EE+1;
                     end
                  elseif ~strcmp(it_id_Enc,it_id_Ret) & strcmp(cat_id_Enc,cat_id_Ret)
                     %disp(['DISC_EE > ' trlij]);    
                     if exist('countDISC_EE')
                         new_disc_ee{countDISC_EE} = [i, j];
                         countDISC_EE = countDISC_EE+1;
                     end
                  elseif ~strcmp(cat_id_Enc,cat_id_Ret)
                     %disp(['DIDC_EE > ' oneListIds{i} '//' oneListIds{j}]);     
                     if exist('countDIDC_EE')
                        new_didc_ee{countDIDC_EE} = [i, j];
                        countDIDC_EE = countDIDC_EE+1;
                     end
                   end
            else
                it_id_Enc = evei(3); it_id_Ret = evej(3); 
                cat_id_Enc = it_id_Enc{1}(1); cat_id_Ret = it_id_Ret{1}(1);
                %if ~strcmp(it_id_Enc,it_id_Ret) & strcmp(cat_id_Enc,cat_id_Ret)
                %    disp(['Trial Overlap > ' it_id_Enc ' > ' it_id_Ret]);
                %end
                 if ~strcmp(cat_id_Enc,cat_id_Ret)
                     disp(['Trial Overlap > ' it_id_Enc ' > ' it_id_Ret]);
                 end 
            end


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
                  if exist('countALL_EM2')
                    new_all_em2{countALL_EM2} = [i, j];
                    countALL_EM2 = countALL_EM2+1;
                  end
                  if strcmp(idEnc,idRet) 
                     %disp(['SISC_EM2> ' oneListIds{i} '//' oneListIds{j}]);  
                     if exist('countSISC_EM2')
                        new_sisc_em2{countSISC_EM2} = [i, j];
                        countSISC_EM2 = countSISC_EM2+1;
                     end
                  end
                   if ~strcmp(idEnc,idRet) & strcmp(oneIdEnc, oneIdRet)
                        %disp(['DISC_EM2 > ' trlij]);    
                        if exist('countDISC_EM2')
                            new_disc_em2{countDISC_EM2} = [i, j];
                            countDISC_EM2 = countDISC_EM2+1;
                        end
                   end
                   if ~strcmp(oneIdEnc, oneIdRet) 
                        %disp(['DIDC_EM2 > ' oneListIds{i} '//' oneListIds{j}]);     
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

   if ~isempty(intersect(c2c, 'M1A'))
      if strcmp(evei(1), '5')
        allIt_I = double(string(char(([evei(13) evei(14) evei(15)])))); allCat_I = floor(allIt_I/100); 
        for j = 1:length(oneListIds)
           evej = strsplit(oneListIds{j});
           if ~strcmp(evei(12), evej(12)) %not same trial
               if strcmp(evej(1), '5')
                    allIt_J = double(string(char(([evej(13) evej(14) evej(15)])))); allCat_J = floor(allIt_J/100); 
                    allItIJ = [allIt_I ; allIt_J]; 
                    allCatIJ = [allCat_I ; allCat_J]; 
                    if length(unique(allCatIJ)) == 3 & length(unique(allItIJ)) == 6
                        %disp (['M2M2 DISC > '  string(trlijIt)]);     % oneListIds{i} '//' oneListIds{j}]);     
                       if exist('countDISC_M1A')
                            new_disc_m1a{countDISC_M1A} = [i, j];
                            countDISC_M1A = countDISC_M1A+1;
                       end
                    end
                    if length(unique(allCatIJ)) == 6 
                        %disp (['All different categories > ' string(trlijIt)]);     
                        if exist('countDIDC_M1A')
                            new_didc_m1a{countDIDC_M1A} = [i, j];
                            countDIDC_M1A = countDIDC_M1A+1;
                        end
                    end
               end
           end
        end
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
                        %disp(['SISC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_sisc_m2m2{countSISC_M2M2} = [i, j];
                        countSISC_M2M2 = countSISC_M2M2+1; 
                      end 
                  end
                  if cuedCI == cuedCJ
                    if exist('countDISC_M2M2')
                        %disp(['DISC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_disc_m2m2{countDISC_M2M2} = [i, j];
                        countDISC_M2M2 = countDISC_M2M2+1;
                    end
                  else
                    if exist('countDIDC_M2M2')
                        %disp(['DIDC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_didc_m2m2{countDIDC_M2M2} = [i, j];
                        countDIDC_M2M2 = countDIDC_M2M2+1;
                    end
                end
              end 
           end
        end
      end
    end



    % % % % % 1) M2M2 for SCSP SCDP
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
                  
                  if cuedCI == cuedCJ & cueRetI == cueRetJ
                    if exist('countSCSP_M2M2')
                        %disp(['SCSP_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_scsp_m2m2{countSCSP_M2M2} = [i, j];
                        countSCSP_M2M2 = countSCSP_M2M2+1;
                    end
                  else
                    if cuedCI ~= cuedCJ & cueRetI == cueRetJ
                        %disp(['DIDC_M2M2 > ' oneListIds{i} '//' oneListIds{j}]);     
                        new_scdp_m2m2{countSCDP_M2M2} = [i, j];
                        countSCDP_M2M2 = countSCDP_M2M2+1;
                    end
                end
              end 
           end
        end
      end
    end



    
end


for ci = 1:length(contr2save)
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
 
 
 

