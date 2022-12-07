function[neuralRDMs] = createNeuralRDMs(cfg, oneListPow)

freqs2test = cfg.freqs; 
win_width = cfg.win_width; 
mf = cfg.mf; 
fR = cfg.fR; 
avTFV = cfg.avTFV;

if fR
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        pow = oneListPow(:, :, f, :);
        nTrials = size(pow, 1); nChans = size(pow, 2); nFreq = size(pow, 3); nTimes = size(pow, 4); 
        ALLER = reshape (pow, [nTrials, nChans*nFreq, nTimes]);
        bins  =  floor ( (nTimes/mf)- win_width/mf+1 );

        parfor timei = 1:bins %parfor possible
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            if avTFV
                x = mean(ALLER(:, :, timeBins), 3, 'omitnan');
            else
                x = ALLER(:, :, timeBins);
                x = reshape(x, nTrials, []); 
            end
            neuralRDMs(:, :, freqi, timei) = corr(x', 'type', 's');
        end

    end
    
else
        f = freqs2test; 
        pow = oneListPow(:, :, f, :);
        nTrials = size(pow, 1); nChans = size(pow, 2); nFreq = size(pow, 3); nTimes = size(pow, 4); 
        ALLER = reshape (pow, [nTrials, nChans*nFreq, nTimes]);
        bins  =  floor ( (nTimes/mf)- win_width/mf+1 );

        parfor timei = 1:bins %parfor possible
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