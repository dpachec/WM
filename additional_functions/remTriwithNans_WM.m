 
function [tr2exc]  = remTriwithNans_WM (oneListTraces, oneListIds, oneListMarkers, xlimE,eLim);
   
                                                
nChans = size(oneListTraces, 1);
nPnts = size(oneListTraces, 2);
nTrials = size(oneListTraces, 3);
tr2exc = zeros(nChans, nTrials);

times = (eLim(1)*1000)+1:1:(eLim(2)*1000);

for chani = 1:nChans
    for ti = 1:nTrials
        if ~isempty(oneListMarkers(:,:, ti))
               time_s = dsearchn(times',xlimE(1));
               time_e = dsearchn(times',xlimE(2));
            data = squeeze(oneListMarkers(chani,time_s:time_e, ti));
            if any(isnan(data))
               % disp (['trial with NaN: Chan ' num2str(chani) ' Trial: ' num2str(ti)]);
                tr2exc(chani, ti) = 1;
            end
        else
            %empty trials
        end
        
    end
end
 
tr2exc1 = logical(sum (tr2exc,1));
disp ([num2str(sum(tr2exc1)) ' of ' num2str(nTrials) ' trials were excluded']);
tr2exc = tr2exc1;
 

 
 
 
 
 
 
 
 
 
 

