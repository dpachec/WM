%% RSA GLOBAL CREATE FOLDERS
% run in folder with gOBO_rsa files from all subj together
% first create the folders based on the conditions's names
sublist = dir('*_rsa.mat'); sublist = {sublist.name};
fname_tmp = 's0';
for filei=1:length(sublist)
   str = sublist{filei};
   fname = str(5:end-13);
   %fname = str(5:end-7);
   if ~strcmp(fname, fname_tmp)
      mkdir(fname)
      movefile(str, fname)
   else
      movefile(str, fname)
   end
   str_tmp = sublist{filei};
   fname_tmp = str_tmp (5:end-13);
   %fname_tmp = str_tmp (5:end-7);
end

%%GO THROUGH FOLDERS RECURSIVELY
% run in folder with subfolders SI_C, SI_F... SI_R
folders = dir(); dirs = find(vertcat(folders.isdir));
folders = folders(dirs);
 
 
for foldi = 3:length(folders) %start at 3 cause 1 and 2 are . and ...
    direct = folders(foldi);
    cd (direct.name)
    sublist = dir('*_rsa.mat');
    sublist = {sublist.name};
    disp (['measurements -> ' num2str(length(sublist))]);
 
    for subji=1:length(sublist)
        load(sublist{subji});
        %chan2plot(subji)
        all{subji,1} = rsaZ; %when more than 1 electrode
        if exist('allIDs')
            all_IDs{subji, :} = allIDs; 
        end
        %all{subji,1} = rsaZ;
 
    end
 
    if exist ('timeBins')
        timeBins1 = timeBins (:, [1, end]);
    end
 
 
    cd .. % goes up one directory
    
    if ismac
        rmdir(direct.name, 's')
    else
        cmd_rmdir(direct.name, 's')
    end

    
    filename = [sublist{subji}(5:end-13)];
    %filename = [sublist{subji}(5:end-12)];
    eval([filename '= all;']);
    if exist('all_IDs')
        save (filename, filename, 'all_IDs', '-v7.3');
    else
        save (filename, filename, '-v7.3');
    end
    %save (filename, filename);
    
    
end












