function EEG = re_reference_WM(EEG, montage, loc)

disp ('>>>>> re-referencing...');

%first remove EEG and POL for chinese data
chs = {EEG.chanlocs(:).labels}';
chs = erase(chs, 'POL'); %deletes over the whole array
chs = erase(chs, 'EEG');
chs = erase(chs, '-Ref');
chs = erase(chs, 'G_');
chs = erase(chs, '-');

% then re-order alphabetically to avoid jumps
[x ids] = natsortfiles(chs);
EEG.chanlocs = EEG.chanlocs(ids);
EEG.data = EEG.data(ids,:);

if strcmp (montage, 'bipo');
    disp ('>> BIPOLAR REFERENCE...');

        for i = 1:length(EEG.chanlocs)-1
            currChan = x{i};nextChan = x{i+1};
            index = find(isletter(currChan));currChan = currChan(index);
            index = find(isletter(nextChan));nextChan = nextChan(index);
            currChanB = x{i};nextChanB = x{i+1};
            indexNum = find(~isletter(currChanB)); currNum = str2num(currChanB(indexNum));
            indexNumNext = find(~isletter(nextChanB)); nextNum = str2num(nextChanB(indexNumNext));

            
            disp([currChan ' ' nextChan]);
            disp(['number >> ' num2str(nextNum) ' ' num2str(currNum)]);

            if  strcmp(currChan, nextChan) & (nextNum - currNum == 1)
                disp(['Equal: ' EEG.chanlocs(i).labels ' = ' EEG.chanlocs(i+1).labels ' - ' EEG.chanlocs(i).labels]);
                EEG.data(i,:) = EEG.data(i+1,:,:)  - EEG.data(i,:,:);
                EEG.chanlocs(i).labels = [EEG.chanlocs(i).labels ' - ' EEG.chanlocs(i+1).labels];
                EEG.chanlocs(i).X = ( EEG.chanlocs(i).X + EEG.chanlocs(i+1).X ) / 2;
                EEG.chanlocs(i).Y = ( EEG.chanlocs(i).Y + EEG.chanlocs(i+1).Y ) / 2;
                EEG.chanlocs(i).Z = ( EEG.chanlocs(i).Z + EEG.chanlocs(i+1).Z ) / 2;
                
            else
                disp('Different');
                ids2rem(i) = 1;
            end 
        end
        EEG.data(end,:) = []; %remove last channel (not bipolarized)
        EEG.chanlocs(end) = [];
        ids2rem = logical (ids2rem);
        EEG.chanlocs(ids2rem) = []; 
        EEG.data(ids2rem, :) = [];
        
        EEG.chanlocs = rmfield(EEG.chanlocs, {'ref', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', 'type', 'urchan'});
    
        
        disp (' >>>> bipolar reference all electrodes');
        
    end
    
% average reference
if strcmp (montage, 'aver');
    disp ('>> AVERAGE REFERENCE...');
    dataRef = EEG.data; %dataRef(chans2exc1, :) = []; EEG.chanlocs(chans2exc1, :) = [];
    EEG_average = mean(dataRef, 1);
    EEG.data = dataRef - EEG_average;
end

%%end function