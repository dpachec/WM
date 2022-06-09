function [mySurro] =  rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, takeElec, ...
                                takeFreq, idxCH, idxF, subji, TG, aVTime)
    

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
    for bini = 1:length(currentContrast{coni})
        %fprintf('\n');
    if ~isempty(currentContrast{coni}) 
        if takeElec & ~takeFreq
            iCH = idxCH(:,subji); 
            all2all = currentContrast{coni}{bini}(:,:,iCH,:,:);
            disp('selected electrodes');
        elseif takeFreq & ~takeElec
            disp('only frequencies');
            iF = idxF(:,subji);
            all2all = currentContrast{coni}{bini}(:,:,:,iF,:);
            f = 1:size(all2all, 4);
        elseif takeElec & takeFreq
            disp('>>>> elec and freqs');
            iCH = idxCH(:,subji);
            iF = idxF(:,subji);
            all2all = currentContrast{coni}{bini}(:,:,iCH,iF,:);
            f = 1:size(all2all, 4);
        else
            all2all = currentContrast{coni}{bini};
            %size(all2all)
        end

        n2s = size(all2all, 1);
        trialN = size(all2all, 1);    
        chanN = size(all2all, 3);


        disp (['Cond ' id '   ' num2str(size(all2all, 1)) ' trials']);
        
      
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

            fprintf('\n');

            rsaZ = zeros (trialN, bins, bins);
            %fprintf('\n'); fprintf('trial correlation:          '); 
            parfor triali = 1:trialN
                mX= squeeze(xM(triali,:,:));
                mY= squeeze(yM(triali,:,:));
                r = corr (mX', mY','Type', 's', 'Rows', 'pairwise'); 
                if (~aVTime)
                    r(r==0) = 10000;
                    r = tril(squeeze(r)); %symmetric so only half saved
                    r(r == 0) = nan;r(r==10000)=0;
                end
                rsaZ(triali, :, :) = atanh(r);
            end
            
           rsaZ(isinf(rsaZ)) = nan;
            
            
          
         allRSA{bini} = rsaZ; 
         


% 
% 
% 
%            zM = zeros (trialN, bins,  chanN * length(f), 2);
%             
%             for timei = 1:bins
%                 %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
%                 timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
%                 
%                 x = mean(all2all(:, 1,:,f, timeBins), 5, 'omitnan');
%                 x = reshape (x, [trialN, chanN * length(f)]);
%                 y = mean(all2all(:, 2,:,f,timeBins), 5, 'omitnan');
%                 y = reshape (y, [trialN, chanN * length(f)]);
% 
%                 zM(:, timei, :,1) = x; 
%                 zM(:, timei, :,2) = y; 
%             end
%             
%             rsaZT = arrayfun(@(i)tril(atanh(corr(squeeze(zM(i, :,:,1))', squeeze(zM(i, :,:,2))','Type', 's'))), 1:size(zM,1), 'un', 0);
%             rsaZT = cat(3, rsaZT {:}); rsaZT = permute(rsaZT, [3 1 2]);
%             rsaZT(rsaZT==0) = nan; rsaZT(isinf(rsaZT)) = nan;
%             allRSA{bini} = rsaZT; 
% 
%        end
    

        
    end
    
 end 
 


 if TG %avoiding the loop over trials is the fastest, but we only store the diagonal 
            filename = ['s' num2str(subji, '%02.f') '_' id '_gOBO'   '_rsa.mat'];
            rsaZ = cat(1, allRSA{:});
            if ~isempty(allIDs) & ndims(rsaZ) == 3
                if aVTime
                    if strcmp(id(end-2:end), 'M2A') | strcmp(id(end-2:end), '2M2')
                        %disp('averaging M2A');
                        %rsaZ = squeeze(mean(rsaZ(:,6:45,6:45), 'all', 'omitnan'));  
                        rsaZ = squeeze(mean(mean(rsaZ(:,6:45,6:45), 2, 'omitnan'), 3, 'omitnan'));  
                        %rsaZ = squeeze(mean(mean(rsaZ(:,6:10,6:10), 2, 'omitnan'), 3, 'omitnan'));  
                        %rsaZ = squeeze(mean(mean(rsaZ(:,11:45,11:45), 2, 'omitnan'), 3, 'omitnan'));  
                    elseif strcmp(id(end-2:end), 'EM2')
                        rsaZ = squeeze(mean(rsaZ(:,11:45,11:45), 'all', 'omitnan'));  
                        %rsaZ = squeeze(mean(mean(rsaZ(:,6:13,6:45), 2, 'omitnan'), 3, 'omitnan'));  
                        %rsaZ = squeeze(mean(mean(rsaZ(:,6:10,6:10), 2, 'omitnan'), 3, 'omitnan'));  
                        %rsaZ = squeeze(mean(mean(rsaZ(:,11:45,11:45), 2, 'omitnan'), 3, 'omitnan'));  
                    else
                        rsaZ = squeeze(mean(rsaZ(:,6:13,6:13), 'all', 'omitnan'));  
                    end
                    save (filename, 'rsaZ'); %, 'timeBins'
                else
                    save (filename, 'rsaZ', 'allIDs'); %, 'timeBins'
                end
            else 
                rsaZ = []
                save (filename, 'rsaZ'); %, 'timeBins'
            end
            
            
           % fprintf('\n');
         
        else
            disp('hola')
            rsaZ = cat(1, allRSA{:});
            parfor triali = 1:size(rsaZ, 1)
                rsaN(triali, :) = diag(squeeze(rsaZ(triali, :, :)));
            end
            rsaZ = rsaN;
            filename = ['s' num2str(subji, '%02.f') '_' id '_dOBO'   '_rsa.mat'];
            save (filename, 'rsaZ', 'allIDs'); %, 'timeBins'
        end
    
 end


 
 
 
 

