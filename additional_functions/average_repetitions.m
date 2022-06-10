function [cfg_contrasts] = average_repetitions(cfg_contrasts)

oneListIds = cfg_contrasts.oneListIds_c;
oneListPow = cfg_contrasts.oneListPow;


% % % first average only 1 2 3 trials at encoding

ids0 = cellfun(@(x) strsplit(string(x)), oneListIds, 'un', 0);
x3c =  cell2mat(cellfun(@(x) str2double(x{3}), ids0, 'un', 0));
fidr = isnan(x3c); 
x3c(fidr) = [];
oneListPow(fidr,:,:,:) = [];
oneListIds(fidr) = [];

[C t] = sort(x3c);
Cpow = oneListPow(t,:,:,:);
Cids = oneListIds(t);

clear dataF oneListPowAv oneListIdsAv
[C1, iaF, icF] = unique(C,'first');   [C1, iaL, icL] = unique(C,'last');
for i = 1:length(C1)
   dataF{i, 1} = C1(i);dataF{i, 2} = iaF(i):iaL(i); 
end           

if exist('dataF')
    for dfi = 1:length(dataF)

       if  length(dataF{dfi,2}) > 1
        oneListPowAv(dfi, :, :, :) = mean(Cpow(dataF{dfi,2}, :, :, :));   
       else
        oneListPowAv(dfi, :, :, :) = Cpow(dataF{dfi,2}, :, :, :);   
       end

       %oneListIdsAv(dfi,:) =  Cids(dataF{dfi,2}(1));
       allRepIds = Cids(dataF{dfi,2});
       if length(allRepIds) > 1
           clear tri iup
           for insi = 1:length(allRepIds)
              iup = allRepIds{insi};
              iup = strsplit(iup);
              tri{insi}= ['t' iup{11} 's' iup{end}];
           end

           sr2t = strsplit(string(Cids(dataF{dfi,2}(1))));
           sr2t{1} = '*';
           sr2t{2} = '*';
           %disp([ 'Length tri >>  ' num2str(length(tri))]);
           tri2 = join(tri,'_');
           sr2t{12}  = char(tri2);
           oneListIdsAv{dfi,:} =  join(sr2t, '  ');

       else 
           str1j = strsplit(string(allRepIds{1})); %store session in id as well for the non-repeated items
           str1j{12} =  ['t' str1j{12} 's' str1j{end}];
           oneListIdsAv{dfi,:} = join(str1j, '  ');
       end
    end

    if size(Cpow, 2) == 1 %one electrode squeezes and creates problems 
        oneListPowAv = permute(oneListPowAv, [1 4 2 3]);
    end

end


% %  % % 2) Average cued items during maintenance

oneListIds = cfg_contrasts.oneListIds_c;
oneListPow = cfg_contrasts.oneListPow;

for i = 1:length(oneListIds)
    eve = strsplit(oneListIds{i});
     %if isnan(str2double(eve{3}))
     if ( strcmp(eve{1}, '7') | strcmp(eve{1}, '9') ) & ~strcmp(eve{2}, '4')
        idx12p = 12 + str2double(eve{2});
        eve{3} = eve{idx12p} ;
        oneListIds{i} = join(eve, '  ');
     else
        eve{3} = 'NaN';
        oneListIds{i} = join(eve, '  '); 
     end
end

ids0 = cellfun(@(x) strsplit(string(x)), oneListIds, 'un', 0);
x3c =  cell2mat(cellfun(@(x) str2double(x{3}), ids0, 'un', 0));

oneListIds(isnan(x3c)) = []; 
oneListPow(isnan(x3c),:,:,:) = []; 
x3c(isnan(x3c)) = []; 

clear C t
[C t] = sort(x3c);
Cpow = oneListPow(t,:,:,:);
Cids = oneListIds(t);


clear dataF oneListPowAv2 oneListIdsAv2
[C1, iaF, icF] = unique(C,'first');   [C1, iaL, icL] = unique(C,'last');
for i = 1:length(C1)
   dataF{i, 1} = C1(i);dataF{i, 2} = iaF(i):iaL(i); 
end      


for dfi = 1:length(dataF)

       if  length(dataF{dfi,2}) > 1
        oneListPowAv2(dfi, :, :, :) = mean(Cpow(dataF{dfi,2}, :, :, :));   
       else
        oneListPowAv2(dfi, :, :, :) = Cpow(dataF{dfi,2}, :, :, :);   
       end

   %oneListIdsAv(dfi,:) =  Cids(dataF{dfi,2}(1));
   allRepIds = Cids(dataF{dfi,2});
   if length(allRepIds) > 1
       clear tri iup
       for insi = 1:length(allRepIds)
          iup = allRepIds{insi};
          iup = strsplit(string(iup));
          tri{insi}= ['t' iup{11} 's' iup{end}];
       end
       sr2t = strsplit(string(Cids(dataF{dfi,2}(1))));
       sr2t{1} = '7*';
       %disp([ 'Length tri >>  ' num2str(length(tri))]);
       tri2 = join(tri,'_');
       sr2t{12}  = char(tri2);
       oneListIdsAv2{dfi,:} =  join(sr2t, '  ');
   else 
       str1j = strsplit(string(allRepIds{1})); %store session in id as well for the non-repeated items
       str1j{12} =  ['t' str1j{12} 's' str1j{end}]; %t= trial; s= session
       oneListIdsAv2{dfi,:} = join(str1j, '  ');
   end
end

if size(Cpow, 2) == 1 %one electrode squeezes and creates problems 
    oneListPowAv2 = permute(oneListPowAv2, [1 4 2 3]);
end


% % % combine 123 and 7 trials
if exist('oneListPowAv')
    oneListPowAv(61:end, :, :, :) = []; 
    oneListPowAv2(61:end, :, :, :) = []; 
    oneListIdsAv(61:end) = []; 
    oneListIdsAv2(61:end) = []; 
    cfg_contrasts = rmfield(cfg_contrasts, 'oneListIds_c'); 
    cfg_contrasts = rmfield(cfg_contrasts, 'oneListPow'); 
    cfg_contrasts.oneListIds_enc = oneListIdsAv;
    cfg_contrasts.oneListIds_maint = oneListIdsAv2;
    cfg_contrasts.oneListPow_enc = oneListPowAv; 
    cfg_contrasts.oneListPow_maint = oneListPowAv2;
else %analysis locked to the probe
   cfg_contrasts.oneListIds_c = oneListIdsAv2;
   cfg_contrasts.oneListPow = oneListPowAv2; 
end





