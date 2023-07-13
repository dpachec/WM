function [mySurro] =  rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, sessi, TG, aVTime)
    

currentContrast = out_contrasts.allContrasts;
currentIds = out_contrasts.allContrastIds;


for i = 1:length(currentContrast)
    if ~isempty(out_contrasts.allContrasts{i})
        nTimepoints = size (out_contrasts.allContrasts{i}{end}, 5); %%all2all file is stored within this cell array
        aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
    end
end

for coni = 1:length(currentContrast)
    bins = aBins(coni);
    allIDs = out_contrasts.allIDs{coni};
    id = currentIds{coni};
    clear allRSA %critical to not merge conditions
    for batchi = 1:length(currentContrast{coni}) % for every batch
        
        if ~isempty(currentContrast{coni}) 
            all2all = currentContrast{coni}{batchi};
        end
    
        trialN = size(all2all, 1);    
        chanN = size(all2all, 3);
    
    
        %disp (['Cond ' id '   ' num2str(size(all2all, 1)) ' trials']);
    
    
        if meanInTime
            if meanInFreq
                xM = zeros (trialN, bins,  chanN );
                yM = zeros (trialN, bins,  chanN );
            else 
                xM = zeros (trialN, bins,  chanN * length(f));
                yM = zeros (trialN, bins,  chanN * length(f));
            end
    
        else 
            if meanInFreq
                xM = zeros (trialN, bins,  chanN * win_width);
                yM = zeros (trialN, bins,  chanN * win_width);
            else
                xM = zeros (trialN, bins,  chanN * length(f) * win_width);
                yM = zeros (trialN, bins,  chanN * length(f) * win_width);
            end
    
        end
                
        for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            if meanInTime
                if meanInFreq
                    x = mean(all2all(:, 1,:,f, timeBins), 5);
                    x = mean(x, 4);
                    %size(x)
                    y = mean(all2all(:, 2,:,f,timeBins), 5);
                    y = mean(y, 4);
                else
                    x = mean(all2all(:, 1,:,f, timeBins), 5, 'omitnan');
                    x = reshape (x, [trialN, chanN * length(f)]);
                    y = mean(all2all(:, 2,:,f,timeBins), 5, 'omitnan');
                    y = reshape (y, [trialN, chanN * length(f)]);
                end
            else
                if meanInFreq
                    x = all2all(:, 1,:,f,timeBins);
                    x = squeeze(x);
                    x = mean(x, 3);
                    x = reshape (x, [trialN, chanN * win_width]);
    
                    y = all2all(:, 2,:,f,timeBins);
                    y = squeeze(y);
                    y = mean(y, 3);
                    y = reshape (y, [trialN, chanN * win_width]);
                else
                    x = all2all(:, 1,:,f,timeBins);
                    x = reshape (x, [trialN, chanN * length(f)* win_width]);
    
                    y = all2all(:, 2,:,f,timeBins);
                    y = reshape (y, [trialN, chanN * length(f)* win_width]);
                end
            end
    
            xM(:, timei, :) =  x;
            yM(:, timei, :) =  y;
            %disp(['size xM >>    ' num2str(size(xM))])
            
        end
    
        
    
        rsaZ = zeros (trialN, bins, bins);
        %fprintf('\n'); fprintf('trial correlation:          '); 
        for triali = 1:trialN
            mX= squeeze(xM(triali,:,:));
            mY= squeeze(yM(triali,:,:));
            r = corr (mX', mY','Type', 's'); 
            idC = strsplit(id, '_');
            if ~strcmp(idC{2}, 'EM2') 
                id0 = find(r==0);
                r = tril(squeeze(r)); %symmetric so only half is saved
                r(r == 0) = nan;r(id0)=0;
            end
            rsaZ(triali, :, :) = atanh(r);
        end
        
        rsaZ(isinf(rsaZ)) = nan;
      
        allRSA{batchi} = rsaZ; 
     
        
    end
    
 
 


    if TG 
        filename = ['s' num2str(sessi, '%02.f') '_' id '_gOBO'   '_rsa.mat'];
        rsaZ = cat(1, allRSA{:});
        if ~isempty(allIDs) & ndims(rsaZ) == 3
            if aVTime
                rsaZ = squeeze(mean(rsaZ(:,6:15,6:15), 'all', 'omitnan'));  
                save (filename, 'rsaZ'); %, 'timeBins'
            else
                save (filename, 'rsaZ', 'allIDs'); %, 'timeBins'
            end
        else 
            rsaZ = []
            save (filename, 'rsaZ'); %, 'timeBins'
        end
        
    else %only store the diagonal 
        rsaZ = cat(1, allRSA{:});
        parfor triali = 1:size(rsaZ, 1)
            rsaN(triali, :) = diag(squeeze(rsaZ(triali, :, :)));
        end
        rsaZ = rsaN;
        filename = ['s' num2str(sessi, '%02.f') '_' id '_dOBO'   '_rsa.mat'];
        save (filename, 'rsaZ', 'allIDs'); %, 'timeBins'
    end



    
    end 
end


 
 
 
 

