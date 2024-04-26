function [all_events allTrlInfo] = loadLogsWM ()
    
    currentFolder = pwd;

    all_pat_folder='D:\_WM\data\iEEG\';
    cd(all_pat_folder);
    all_pat = dir('sub*');all_pat = {all_pat.name};

    %%read in logfiles and create trialinfo for each session
    % loop sub, sess, logfiles
    for sub=1:size(all_pat,2)
        sel_sub=char(all_pat(sub));
        fname1 = string(strcat(all_pat_folder,sel_sub));
        cd(fname1);
        all_sess=dir('Sess*');all_sess = {all_sess.name};
        
        for ses=1:size(all_sess,2)
            sel_sess=all_sess(ses);
            fname = string(strcat(all_pat_folder,sel_sub,'/',sel_sess,'/Beh'));
            cd(fname);
            all_logs = dir('*log');all_logs = {all_logs.name};
            for log=1:size(all_logs, 2)        
                sel_log=char(all_logs(log)); 
                sel_log(isspace(sel_log))=[]; % delete blanks for better indei
                [trlinfo,datelog]=log2trlinfo(sel_log);
                fname2 = strcat(sel_sub(1:5),'_',sel_log(1:end-4), '_trialinfo');
                save(fname2,'trlinfo','datelog');
            end
        end
    end
    cd(all_pat_folder);

    
    
for sub=1:size(all_pat,2)
    sel_sub=all_pat(sub);
    fname = char(strcat(all_pat_folder,sel_sub));
    cd(fname);
    all_sess=dir('Sess*');all_sess = {all_sess.name};
    countWM = 0;
    for ses=1:size(all_sess, 2)
        sel_sess=all_sess(ses);
        fname = string(strcat(all_pat_folder,sel_sub,'/',sel_sess,'/Beh'));
        cd(fname);
        trl_files=dir('*trialinfo.mat'); trl_files = {trl_files.name};
        % load trialinfo files
        for blo=1:size(trl_files,2)
            
            sel_log=char(trl_files(blo));
            load(sel_log);
            % construct trigger channels from trlinfo
            % 1: block_type 2=WM
            % 2: WM block
            % 3: trialnumber
            % 4:9 enco item & times
            % 10:11cue & times
            % 12:15: test item token + time ( at test always 6 items are presented, the three enco items plus three added items (12:14)
            % 16:21 which response codes for which item (randvec: position of x in randvec codes for item: ie randvec(4)=2 means response 4 means enco item 2, also item shown at position 4 on test screen)
            % 22-27 recode randvec (column 22: correct response for item 1, 25: correct response for item 2 etc...)
            % 28-30: responses
            % 31-33: lure pos (there are always some some similar lures in the response display, this codes on which position similar lures to the respective enco item 1-3 are shown, should make sense when you check the token vs type correct)
            % 34-36. all correct (token correct)
            % 37:39: correct on type level
            % 40: completely correct all token
            % 41: completely correct all types
            % 42-44: RT for response relative to test onset
            % trlinfo(:, 48) = m1_time';
            % trlinfo(:, 49) = m2_time';
            % trlinfo(:, 50) = fix_time(:,1); 
            % trlinfo(:, 51) = fix_time(:,2); 
            
            
            if ~isempty(strfind(sel_log,'wmblock'));
                
                countWM = countWM+1;
                subsessblo(sub, 1) = countWM;
        
                
                trlinfo = getPresentedItems(trlinfo);
                log_type='wmblock';
                
                clear c1 c2 c3 c4 c5 c6 c7 c8 c91011 event_info
                c1 = [trlinfo(:,7);trlinfo(:,8);trlinfo(:,9);trlinfo(:,11); ...
                                   trlinfo(:,15); ... 
                                   trlinfo(:,15) + trlinfo(:,42); ...
                                   trlinfo(:,15) + trlinfo(:,43); ...
                                   trlinfo(:,15) + trlinfo(:,44); ...
                                   trlinfo(:,48); trlinfo(:,49); ...
                                   trlinfo(:,50); trlinfo(:,51)];
               n= 12;
               %c2 = trial type
               c2 = 1:n; % events in a trial
               c2 = repmat(c2,1,60)';   
               c3 = trlinfo(:, 10); % column 10 and not n ! 
               c3 = repelem(c3,n,1);

               
               c4 = [trlinfo(:,4);trlinfo(:,5);trlinfo(:,6);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,45); ...
                                   trlinfo(:,46); ...
                                   trlinfo(:,47); ...
                                   trlinfo(:,6); NaN(60,1); ... 
                                   trlinfo(:,4); trlinfo(:,5)]; % start of m1 is end of pic3, trlinfo(:,6) 
                                   %shoould go after trlinfo(:,5) when done properly 
               c5 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,28); ...
                                   trlinfo(:,29); ...
                                   trlinfo(:,30); ...
                                   NaN(60,1); NaN(60,1); ... 
                                   NaN(60,1); NaN(60,1)];

               c6 = [ trlinfo(:,22);trlinfo(:,23);trlinfo(:,24);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,22); ...
                                   trlinfo(:,23); ...
                                   trlinfo(:,24); ...
                                   NaN(60,1); NaN(60,1); ... 
                                   NaN(60,1); NaN(60,1)];
                c7 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,34); ...
                                   trlinfo(:,35); ...
                                   trlinfo(:,36); ...
                                   NaN(60,1); NaN(60,1); ... 
                                   NaN(60,1); NaN(60,1)];
                c8 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,37); ...
                                   trlinfo(:,38); ...
                                   trlinfo(:,39); ...
                                   NaN(60,1); NaN(60,1); ... 
                                   NaN(60,1); NaN(60,1)];
                               
                c9 = trlinfo(:, 40);
                c9 = repelem(c9,n,1);
                c10 = trlinfo(:, 41); 
                c10 = repelem(c10,n,1);
                
                c11 = trlinfo(:, 1); 
                c11 = repelem(c11,n,1);
                c12 = trlinfo(:, 3);  %or (1:60)
                c12 = repelem(c12,n,1);
                c13 = trlinfo(:, 3); 
                c13 = repelem(c13,n,1);
                 
                c141516 = trlinfo(:, 4:6);
                c141516 = repelem(c141516,n,1);
                c171819 = trlinfo(:, 12:14);
                c171819 = repelem(c171819,n,1); 
                
                
                c21 = trlinfo(:, 40);
                c21 = repelem(c21,n,1);
                c22 = trlinfo(:, 41);
                c22 = repelem(c22,n,1);

                ctrlinf1630 = trlinfo(:, 22:30);
                ctrlinf1630 = repelem(ctrlinf1630, n, 1);
                
                
                [c1 idx]=sort(c1);
                c4 = c4(idx);
                c5 = c5(idx);
                c6 = c6(idx);
                c7 = c7(idx);
                c8 = c8(idx);
                                
                event_info = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c141516 c171819 c21 c22 ctrlinf1630];
                
                [uni idx] = unique(c1);
                event_info = event_info(idx, :);
                
                %%event by type
                %3 = cue; 
                %4= id(1-3 = of pres, 6-8 of respnse)
                %5 = 6-8position response at test
                %6 = position of item at test (repeated in 1-3 just in case)
                %7 = correct / incorrect item
                %8 = correct / incorrect type
                %9 = all correct type
                %10 = all correct token
                %11 = block type (2 = WM)
                %12 = block number
                %13 = trial number
                %14-19 = all presented items (14-16 = presented, 16-19 new)
                %20 = latency of response with respect to test onset
                
                
                % extract events
                %disp ('building events');

                %disp (['sub ' num2str(sub),' ses ' num2str(ses), ' blo ' num2str(blo)])
                
                % build diff_12
                clear diff_12_prf
                diff_12_prf(1, 1, 4) = 24172;       %> 1     > sub01_100918
                diff_12_prf(2, 1, 5) = 13040;       %> 2     > sub03_260918
                diff_12_prf(2, 1, 6) = 29020;       %> 3     > sub03_260918
                diff_12_prf(3, 1, 4) = 20362;       %> 4     > sub04_091018
                diff_12_prf(3, 2, 2) = 40400;       %> 5     > sub04_091018
                diff_12_prf(4, 1, 4) = 73482;       %> 6     > sub05_201018
                diff_12_prf(4, 2, 3) = 6664;        %> 7     > sub05_201018
                diff_12_prf(5, 1, 4) = 16834;       %> 8     > sub06_191118
                diff_12_prf(5, 1, 5) = 19785;       %> 9     > sub06_191118
                diff_12_prf(6, 1, 5) = 21036;       %> 10    > sub07_301118
                diff_12_prf(6, 1, 6) = 7854;        %> 11    > sub07_301118
                diff_12_prf(6, 1, 7) = 11929;       %> 12    > sub07_301118
                diff_12_prf(6, 2, 3) = 7835;        %> 13    > sub07_301118
                diff_12_prf(6, 2, 4) = 5933;        %> 14    > sub07_301118
                diff_12_prf(7, 1, 4) = 13846;       %> 15    > sub11_030219
                diff_12_prf(7, 1, 5) = 49861;       %> 16    > sub11_030219
                diff_12_prf(7, 2, 3) = 7616.8;      %> 17    > sub11_030219
                diff_12_prf(7, 2, 4) = 6168.8;      %> 18    > sub11_030219
                diff_12_prf(8, 1, 4) = 61497;       %> 19    > sub12_100219
                diff_12_prf(8, 1, 5) = 25824;       %> 20    > sub12_100219
                diff_12_prf(9, 1, 4) = 6957.3;      %> 21    > sub13_210219
                diff_12_prf(9, 1, 5) = 4406;        %> 22    > sub13_210219
                diff_12_prf(10, 1, 4) = 38012;      %> 23    > sub14_040319
                diff_12_prf(10, 1, 5) = 36989;      %> 24    > sub14_040319
                diff_12_prf(11, 1, 4) = 17572;      %> 25    > sub18_031520
                diff_12_prf(11, 1, 5) = -32909;     %> 26    > sub18_031520
                diff_12_prf(12, 1, 5) = 22498;      %> 27    > sub19_230520
                diff_12_prf(12, 1, 6) = 2101.2;     %> 28    > sub19_230520
                diff_12_prf(13, 1, 5) = 1568.6;     %> 29    > sub20_100720
                diff_12_prf(13, 1, 6) = 59345;      %> 30    > sub20_100720
                diff_12_prf(14, 1, 5) = 2490.4;     %> 31    > sub21_160720
                diff_12_prf(14, 1, 6) = 464110;     %> 32    > sub21_160720
                diff_12_prf(15, 1, 6) = 10988;      %> 33    > sub22_061120
                diff_12_prf(15, 1, 7) = 26472;      %> 34    > sub22_061120
                diff_12_prf(15, 1, 8) = -11351;     %> 35    > sub22_061120
                diff_12_prf(15, 2, 4) = 1003;       %> 36    > sub22_061120
                diff_12_prf(15, 2, 5) = 6621.6;     %> 37    > sub22_061120
                diff_12_prf(16, 1, 4) = 25509;      %> 38    > sub23_061220
                diff_12_prf(16, 1, 5) = 13282;      %> 39    > sub23_061220
                diff_12_prf(17, 1, 6) = 51033;      %> 40    > sub24_110321
                diff_12_prf(17, 1, 7) = 2.3242e+04; %> 41    > sub24_110321
                diff_12_prf(17, 1, 8) = 1.2168e+04; %> 42    > sub24_110321
                diff_12_prf(18, 1, 5) = 1.7303e+04; %> 43    > sub25_220321
                diff_12_prf(18, 1, 6) = 1.0182e+04; %> 44    > sub25_220321
                diff_12_prf(18, 2, 4) = 6.4681e+04; %> 45    > sub25_220321
                diff_12_prf(18, 2, 5) = 8076;       %> 46    > sub25_220321
                
                diff_12_prf(19, 1, 4) = 27278;      %> 47    > sub31_
                diff_12_prf(19, 1, 5) = 17857;      %> 48    > sub31_
                diff_12_prf(20, 1, 5) = 12628;      %> 49    > sub32_
                diff_12_prf(20, 1, 6) = 64796;      %> 50    > sub32_
                diff_12_prf(21, 1, 5) = 15068;      %> 51    > sub33_
                diff_12_prf(21, 1, 6) = 7418.8;     %> 52    > sub33_
                diff_12_prf(22, 1, 5) = 29531;      %> 53    > sub34_
                diff_12_prf(22, 1, 6) = 25955;      %> 54    > sub34_
                diff_12_prf(23, 1, 5) = 25851;      %> 55    > sub35_
                diff_12_prf(23, 1, 6) = 15358;      %> 56    > sub35_
                diff_12_prf(24, 1, 5) = 32825;      %> 57    > sub36_
                diff_12_prf(24, 1, 6) = 46942;      %> 58    > sub36_
                diff_12_prf(25, 1, 5) = -3900.4;    %> 59    > sub37_
                diff_12_prf(25, 1, 6) = -11319;     %> 60    > sub37_                
                diff_12_prf(26, 1, 5) = -14026;     %> 61    > sub38_
                diff_12_prf(26, 1, 6) = -4045.6;    %> 62    > sub38_
                diff_12_prf(27, 1, 5) = -9210.6;    %> 63    > sub39_
                diff_12_prf(27, 1, 6) = 13969;      %> 64    > sub39_
                
                diff_12_prf(28, 1, 7) = -25189;     %> 65    > sub51_
                diff_12_prf(28, 1, 8) = -23465;     %> 66    > sub51_
                diff_12_prf(28, 1, 9) = -37847;     %> 67    > sub51_
                diff_12_prf(28, 1, 10) = -47219;    %> 68    > sub51_
                diff_12_prf(29, 1, 6) = -23451;     %> 69    > sub52_
                diff_12_prf(29, 1, 7) = 306270;     %> 70    > sub52_
                diff_12_prf(29, 1, 8) = -13927;     %> 71    > sub52_
                diff_12_prf(30, 1, 7) = -15666;     %> 72    > sub53_
                diff_12_prf(30, 1, 8) = 69429;      %> 73    > sub53_
                diff_12_prf(30, 1, 9) = -59014;     %> 74    > sub53_
                diff_12_prf(30, 1, 10) = -63686;    %> 75    > sub53_
                diff_12_prf(31, 1, 7) = -67823;     %> 76    > sub54_
                diff_12_prf(31, 1, 8) = -12403;     %> 77    > sub54_
                diff_12_prf(31, 1, 9) = -26549;     %> 78    > sub54_
                diff_12_prf(31, 1, 10) = -70942;    %> 79    > sub54_
                diff_12_prf(32, 1, 7) = -24069;     %> 80    > sub55_
                diff_12_prf(32, 1, 8) = -43360;     %> 81    > sub55_
                diff_12_prf(32, 1, 9) = -91856;     %> 82    > sub55_
                diff_12_prf(32, 1, 10) = -34050;    %> 83    > sub55_
                
               
                diff_12 =  diff_12_prf(sub, ses, blo); 
                
                EEG.event = struct('latency', [], 'type', '', 'urevent', []);
                for i = 1:length(event_info)
                    EEG.event(length(EEG.event) + 1) = struct('latency', event_info(i)  /5 - diff_12, 'type', ...
                        num2str(event_info(i,2:end)),'urevent', 0);
                end
                EEG.event(1) = [];
                

                EEG.event = nestedSortStruct(EEG.event, 'latency', 'type');
                 
                x = ([EEG.event.latency] / 2)'; %2 for 1000, 4 for 500Hz
                x = num2cell(x');

                EEG.event = struct('latency',x, 'type', {EEG.event.type});

                allEventInfo{sub}{ses}{blo} = event_info;
                alleeg_event{sub}{ses}{blo} = EEG.event;
                allTrlInfo{sub}{ses}{blo} = trlinfo;
                
            elseif ~isempty(strfind(sel_log,'practice'));    
                log_type='practice';
                event_info = [];
            elseif ~isempty(strfind(sel_log,'decoding'));
                %log_type='decoding';
                %trig_log=trlinfo(:,8);
                %event_info =[];
            else
            end
          
        end
             
    end
end

%UNNEST alleeg_event
count = 1;
for subji = 1:length(alleeg_event)
    allevh = alleeg_event{subji};
    for sesi = 1:length (allevh)
        allevsess = allevh{sesi};
        for bloi = 1:length(allevsess) 
           %disp (['sub ' num2str(subji) ' sess ' num2str(sesi) ' block ' num2str(bloi)]);
           all_events{count,:} = allevsess{bloi}; 
           count = count +1;
        end 
    end  
end

all_events = all_events(~cellfun('isempty',all_events));


cd (currentFolder)


disp('logs imported');





