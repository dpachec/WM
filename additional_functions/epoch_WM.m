function [oneListTraces oneListIds oneListMarkers] = epoch_WM (EEG, eLim)

disp ('>>>>> epoching...');

countE =1; 
EEGM = EEG;EEGM.data = EEG.markers;

for i = 1:length(EEG.event)
    eve = strsplit(EEG.event(i).type);
    
    if  strcmp(eve(1), '1') | strcmp(eve(1), '3')  | strcmp(eve(1), '5')  | strcmp(eve(1), '7') 
        EEG_b = pop_epoch( EEG, {EEG.event(i).type}, eLim, 'newname', ...
            'verbose', 'no', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
         EEGM_b = pop_epoch( EEGM, {EEG.event(i).type}, eLim, 'newname', ...
             'verbose', 'no', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
        traces{countE} = EEG_b.data;
        oneListIds{countE,:} = EEG.event(i).type;
        markers{countE} = EEGM_b.data;
        countE = countE+1;
               
    end
end

%find empty epochs at the beginning and end of the exp
ids = find(cellfun(@isempty,traces));
oneListIds(ids) = [];
traces(ids) = []; 
markers(ids) = []; 

oneListTraces = cat(3,traces{:});  
oneListMarkers = cat(3, markers{:});  


disp ('>>>>> epoching done');


