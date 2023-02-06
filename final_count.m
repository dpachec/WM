
cLen1 = sum(cellfun(@(x) size(x, 1), chanNames_all_VVS))  ;
cLen2 = sum(cellfun(@(x) size(x, 1), chanNames_all_PFC));

cLen1+cLen2


%% 
[6 7 9 13 18 19 20 21 22 23 27 28] 
[1 2 3  5  9 10 11 12 13 14 15 16]

[4 6 7 8 ]


c = [23 6 13 14 8 18 34 19 17 22 27 22 19 19 15 25 28 6 24 36 48 3 41 7 14 5 8 6 3 19 2 37]

length(c)
sum(c)
mean(c)
std(c)

%% 
clear

cd D:\_WM\analysis\out_contrasts

load vvs_elec
VVS = chanNames_all; 

count = 1; 
for subji = 1:length(VVS)
    
    chS = VVS{subji}; 
    for chani = 1:size(chS, 1)
        chNames(count, :) = chS{:, 5}; 
        count = count+1;
    end

end

chNamVVS = unique(chNames)


%% check number of channels 
clear

cd D:\_WM\analysis\out_contrasts

load pfc_elec
PFC = chanNames_all; 

count = 1; 
for subji = 1:length(PFC)
    
    chS = PFC{subji}; 
    for chani = 1:size(chS, 1)
        chNames(count, :) = chS{:, 5}; 
        count = count+1;
    end

end

chNamPFC = unique(chNames)





%% check number of trials 
clear 
%cd D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs
cd D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc

sublist = dir('*mat'); sublist = {sublist.name}'

for subji = 1:length(sublist)

    load (sublist{subji})
    ids = cfg_contrasts.oneListIds_c; 
    
    ids0 = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    x3c =  cell2mat(cellfun(@(x) str2double(x{11}), ids0, 'un', 0));
    
    nTrials(subji, :) = length(x3c);

end




%% numbner of sessions and blocks 
clear
cd ('D:\_WM\data\iEEG\set - monopolar')
sessList = dir('*mat'); sessList = {sessList.name}'

subN = 0; 
sessN = 0; 
bloN = 0; 

subStrOld = 's00';
for sessi = 1:length(sessList)

    subStrNew = sessList{sessi}(1:3);

    if ~strcmp(subStrOld, subStrNew)
        disp ('new!')
        subN = subN+1;
        subStrOld = subStrNew ; 
    end


    
end

%% count number of blocks
clearvars -except sessList
for i = 1:length(sessList)
    currentFile = sessList{i};
    
    subjects{i, 1} = currentFile(1:3);
    subjects{i, 2} = currentFile(7:9);
    subjects{i, 3} = currentFile(14);

end


[C, iaF, icF] = unique(subjects(:, 1),'first');
[C, iaL, icL] = unique(subjects(:, 1),'last');

for i = 1:length(C)
    dataF{i, 1} = C(i);
    dataF{i, 2} = length(iaF(i):iaL(i));%this works because unique sorts the items by name
end   

mean(double([dataF{:,2}]))

%% count number of sessions
clearvars -except sessList
for i = 1:length(sessList)
    currentFile = sessList{i};
    %count sessions
    subjects{i, 1} = currentFile(1:9);
end


[C, iaF, icF] = unique(subjects(:, 1),'first');
[C, iaL, icL] = unique(subjects(:, 1),'last');
for i = 1:length(C)
    dataF{i, 1} = C(i);
    dataF{i, 2} = length(iaF(i):iaL(i));%this works because unique sorts the items by name
end   


for i = 1:length(dataF)
    newSubjects(i, :) = string(dataF{i, 1}{1}(1:3));

end

clear dataF
[C, iaF, icF] = unique(newSubjects,'first');
[C, iaL, icL] = unique(newSubjects,'last');
for i = 1:length(C)
    dataF{i, 1} = C(i);
    dataF{i, 2} = length(iaF(i):iaL(i));%this works because unique sorts the items by name
end   



mean(double([dataF{:,2}]))
std(double([dataF{:,2}]))

%% 
clear 
cd D:\_WM\analysis\out_contrasts\raw_traces\allTrials

load untitled
load untitled2
load untitled3

dataF{1 , 3} = nTrialsVVS([1])
dataF{2 , 3} = nTrialsVVS([2])
dataF{3 , 3} = nTrialsVVS([3])
dataF{4 , 3} = nTrialsVVS([4])
dataF{5 , 3} = nTrialsVVS([5])
dataF{6 , 3} = nTrialsVVS([6])
dataF{7 , 3} = nTrialsVVS([7])
dataF{8 , 3} = nTrialsVVS([8])
dataF{9 , 3} = nTrialsVVS([9])
dataF{11 , 3} = nTrialsVVS([10])
dataF{12 , 3} = nTrialsVVS([11])
dataF{13 , 3} = nTrialsVVS([12])
dataF{14 , 3} = nTrialsVVS([13])
dataF{15 , 3} = nTrialsVVS([14])
dataF{16 , 3} = nTrialsVVS([15])
dataF{18 , 3} = nTrialsVVS([16])
dataF{21 , 3} = nTrialsVVS([17])
dataF{22 , 3} = nTrialsVVS([18])
dataF{23 , 3} = nTrialsVVS([19])
dataF{24 , 3} = nTrialsVVS([20])
dataF{25 , 3} = nTrialsVVS([21])
dataF{26 , 3} = nTrialsVVS([22])
dataF{27 , 3} = nTrialsVVS([23])
dataF{28 , 3} = nTrialsVVS([24])
dataF{29 , 3} = nTrialsVVS([25])
dataF{30 , 3} = nTrialsVVS([26])
dataF{31 , 3} = nTrialsVVS([27])
dataF{32 , 3} = nTrialsVVS([28])





dataF{10 , 3} = nTrials_PFC([4])
dataF{17 , 3} = nTrials_PFC([6])
dataF{19 , 3} = nTrials_PFC([7])
dataF{20 , 3} = nTrials_PFC([8])







x = [dataF{:, 3}] ./ [dataF{:, 2}]
x = x'


x1 = 240-x; 
mean(x1)
std(x1)




%%



fileNames = sessList; 



% Get the number of subjects
subjects = unique(cellfun(@(x) x(1:3), fileNames, 'UniformOutput', false));
numSubjects = length(subjects);
finalTable = subjects; 

% Get the number of sessions for each subject
sessionCount = cell(numSubjects, 1);
for i = 1:numSubjects
    sessions = unique(cellfun(@(x) x(5:9), fileNames(strcmp(cellfun(@(x) x(1:3), fileNames, 'un', 0), subjects{i})), 'un', 0))
    sessionCount{i} = length(sessions);   
end
finalTable(:, 2) = sessionCount; 

% Get the number of blocks for each session for each subject
blockCount = cell(numSubjects,1);
for i = 1:numSubjects
    blocks = unique(cellfun(@(x) x(5:14), fileNames(strcmp(cellfun(@(x) x(1:3), fileNames, 'un', 0), subjects{i})), 'un', 0))
    blockCount{i} = length(blocks); 
end

finalTable(:, 3) = blockCount; 





























