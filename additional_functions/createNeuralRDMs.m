function[neuralRDMs] = createNeuralRDMs(cfg, cfg_contrasts)

oneListPow = cfg_contrasts.oneListPow;
freqs2test = cfg.freqs; 
win_width = cfg.win_width; 
mf = cfg.mf; 
fR = cfg.fR; 
avTFV = cfg.avTFV;


if fR
% %     for freqi = 1:length(freqs2test)
% %         f  = freqs2test(freqi);
% %         pow = oneListPow(:, :, f, :);
% %         nTrials = size(pow, 1); nChans = size(pow, 2); nTimes = size(pow, 4); 
% %         ALLER = pow;
% %         bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
% % 
% %         parfor timei = 1:bins 
% %             %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
% %             timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
% %             if avTFV
% %                 x = mean(ALLER(:, :, timeBins), 3, 'omitnan');
% %             else
% %                 x = ALLER(:, :, timeBins);
% %                 x = reshape(x, nTrials, []); 
% %             end
% %             neuralRDMs(:, :, freqi, timei) = corr(x', 'type', 's');
% %         end
% % 
% %     end
        
        % Pre-allocate output array with appropriate dimensions
        pow = oneListPow(:, :, freqs2test, :);
        nTrials = size(pow, 1); nChans = size(pow, 2); nTimes = size(pow, 4); nFreqs = length(freqs2test); 
        bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
        neuralRDMs = zeros(nTrials, nTrials, length(freqs2test), bins);
        timeBins = (1:bins)' * mf - (mf-1) + (0:win_width-1);
        x = pow(:, :, :, timeBins);
        x = reshape (x,  nTrials, size(timeBins, 1)* size(timeBins, 2)* nChans* nFreqs); 
        

        % Calculate RDMs using vectorized correlation
        myRho = corr(x, 'type', 's');
        neuralRDMs = reshape(myRho, nTrials, nTrials, nFreqs, nTimes);
    

else
        f = freqs2test; 
        pow = oneListPow(:, :, f, :);
        nTrials = size(pow, 1); nChans = size(pow, 2); nFreq = size(pow, 3); nTimes = size(pow, 4); 
        ALLER = reshape (pow, [nTrials, nChans*nFreq, nTimes]);
        bins  =  floor ( (nTimes/mf)- win_width/mf+1 );

        for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            if avTFV
                x = mean(ALLER(:, :, timeBins), 3, 'omitnan');
            else
                x = ALLER(:, :, timeBins);
                x = reshape(x, nTrials, []); 
            end
            neuralRDMs(:, :, timei) = corr(x', 'type', 's');
        end

    
end

   
end