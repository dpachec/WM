function [trlinfo, datelog]=log2trlinfo(sel_log)       

loginfo=log2mat(sel_log);
log_tmp=textread(sel_log,'%s');
logtime=log_tmp(7:8);
datelog = datetime(strcat(logtime{1},logtime{2}),'InputFormat','MM/dd/yyyyHH:mm:ss');
clear logtime log_tmp


if ~isempty(strfind(sel_log,'wmblock'));         
log_type='wmblock';
elseif ~isempty(strfind(sel_log,'practice'));    
log_type='practice';
elseif ~isempty(strfind(sel_log,'decoding'));
log_type='decoding';
else
end



switch log_type
     case 'decoding'

        trlinfo =  [];

    case 'wmblock'
        % which info is needed: 
        % 1: block_type 2=WM
        % 2: WM block
        % 3: trialnumber
        % 4:9 enco item & times
        % 10:11cue & times
        % 12:15: test item token + time
        % 16:21 which response codes for which item (randvec: position of x in randvec codes for item: 
        % randvec(4)=2 means response 4 means enco item 2) 
        % 22-27 recode randvec (column 21: correct response for item 1, 22:correct response for item 2 etc...)
        % 28-30: responses 
        % 31-33: lure pos
        % 34-36. all correct (token correct)
        % 37-39: correct on type level
        % 40: completely correct all token
        % 41: completely correct all types
        % 42-44: RT for response relative to test onset
            
        % check which randomization
        %sel_log
        %sel_log(1)
        rand_sub=str2double(sel_log(1));
        % check which block
        block_num=str2double(sel_log(end-4));
        all_block_def=[1:60;61:120;121:180;181:240];
        block_def=all_block_def(block_num,:);
        clear all_block_def
        % load item_mat
        %load('E:\WM_paradigm\iEEG\WM_Freiburg\WM_paradigm_freiburg\WM_seq\item_mat_iEEG1.mat')
        
%         if strcmp(chi, 'chi')
%             %disp('loading data from china');
%             load('D:\_WM\data\item_mat_iEEGchina1.mat');
%         else
%             load('D:\_WM\data\item_mat_iEEG1.mat')
%         end
%         
        % get relevant part of item_mat
        %item_seq=squeeze(item_mat(rand_sub,block_def,:));
        %cue_seq=squeeze(allsubs_trials_token_mat(rand_sub,block_def,5))';
        
        %item_seq=squeeze(item_mat(rand_sub,block_def,:));
        %item_seq_before =squeeze(item_mat(rand_sub,block_def,:));
        
        clear allsubs_trials_token_mat item_mat block_def 
        
        logs=loginfo;
        clear loginfo;
         loginfo=logs{1,3};
         loginfo_res=logs{1,2};
         loginfo_time=logs{1,4};
         
        % get rid of enter responses
    enter_ind=find(strncmp('7',loginfo,1));
    loginfo(enter_ind)=[];
    loginfo_res(enter_ind)=[];
    loginfo_time(enter_ind)=[];
    
   % find start of every trial
    trial_ind=find(strncmp('trial',loginfo,5));
    trial_ind(end+1)=length(loginfo); %add filler at the end

    res123=zeros(1,numel(trial_ind)-1);
    resall=zeros(numel(trial_ind)-1,3);
    rtall=zeros(numel(trial_ind)-1,3);
    rt123=zeros(numel(trial_ind)-1,3);
% get responses and score them
    for t=1:(numel(trial_ind)-1)
         %get all trial indexes
         borders=trial_ind(t):trial_ind(t+1);
         rel_trial=loginfo(borders);
         % get responses in trial
         res_trial=loginfo_res(borders);
         ind_res=find(strncmp('Response',res_trial,8));
         time_trial=loginfo_time(borders);
         % get rand vec (randomization of test items on screen
                % randvec are 6 number, 1-3 code of position of enco item, 4 is
                % new type, 5-6 similar lures
         rand_vec_str=rel_trial(find(strncmp('m2randseq',rel_trial,9)));
             for v=1:6; rand_vec(t,v)=str2num(rand_vec_str{1}(10+(v-1))); end

         % cue type
         cue_str=rel_trial(find(strncmp('cue', rel_trial,3)));
         %cue_type=find(strncmp(cue_str,{'cue1','cue2','cue3','cueall'},4));
         cue_type = double(string(cue_str{1}(4)));
         if cue_type ~= 1 & cue_type ~= 2 & cue_type ~= 3  
             cue_type = 4;
         end
         trial_type(t)=cue_type;
         pic_ind=find(strncmp('pic', rel_trial,3));
         fix_ind = find(strncmp('fix', rel_trial,3));
         test_onset=pic_ind(4);
         cue_time(t)=time_trial((strncmp('cue', rel_trial,3)));
         m1_time(t)=time_trial((strncmp('m1', rel_trial,3)));
         m2_time(t)=time_trial((strncmp('m2', rel_trial,2)));
         enco_time(t,:)=time_trial(pic_ind(1:3));
         lengH = cellfun(@length, rel_trial);
         [u2h u2c] = max(lengH);
         item_secB = rel_trial(u2c); % sometimes answered are logged, but max size is always the line we have to pick
         item_secBB = item_secB{1}([5:6 ' ' 13:14 ' ' 21:22 ' ' 29:30 ' ' 37:38 ' ' 45:46 ]);
         item_secBB1 = string([item_secBB(1) '0' item_secBB(2) ' ' item_secBB(4) '0' item_secBB(5) ' '  item_secBB(7) '0' item_secBB(8) ' ' ...
                        item_secBB(10) '0' item_secBB(11) ' ' item_secBB(13) '0' item_secBB(14) ' '  item_secBB(16) '0' item_secBB(17) ]);
         item_secBBB = double(strsplit(item_secBB1, ' '));
         item_secBBB(logical(~mod(item_secBBB, 100))) = item_secBBB(logical(~mod(item_secBBB, 100))) + 10;
         item_sec4B(t, :) = item_secBBB;
         
         cueB = rel_trial((strncmp('cue', rel_trial,3)));
         cueBB = double(string(cueB{1}(4))); 
         if cueBB > 3 | isnan(cueBB) cueBB = 4; end
         cue_seq(t,:) = cueBB;

         fix_time(t,:)=time_trial(fix_ind(1:2));
         test_time(t)=time_trial(test_onset);
         % number of responses in trial
         num_res=numel(ind_res);    
         if cue_type < 4 
             if num_res==1
                res123(t)=str2num(rel_trial{ind_res});
                rt123(t,cue_type)=time_trial(ind_res)-time_trial(test_onset);
             elseif num_res>=1
                 % count too much response
                 tmr123(t)=1;
                    %use first response
                  res123(t)=str2num(rel_trial{ind_res(1)});
                  res2_123(t)=str2num(rel_trial{ind_res(2)});
                  rt123(t,cue_type)=time_trial(ind_res(1))-time_trial(test_onset);

             elseif num_res==1
                 % count no response trials
                 res123(t)=0;  
                 rt123(t)=NaN;
             end 
         elseif cue_type==4
             if num_res==3
                % keep the three responses
                 resall(t,:)=cellfun(@str2num,rel_trial(ind_res))'; 
                 rtall(t,:)=time_trial(ind_res)-repmat(time_trial(test_onset),3,1);

             elseif num_res<=2 & num_res>0
                 % less than three responses means one reponse was missed
                 % missed responses are characterized by no response between
                 % feedback (wrong-wrong)
                 right=find(strncmp('ok',rel_trial,2));
                 wrong=find(strncmp('wrong',rel_trial,5));
                 fb_ind=sort([right;wrong]);
                 % when difference between fb_ind=1, then there is missed
                 % response, if no diff ==1, first response was missed
                    if ind_res(end)>fb_ind(end)
                       ind_res(end)=[];
                       num_res=numel(ind_res);
                    end

                 miss_res23=find(diff(fb_ind)==1);
                     if isempty(miss_res23) && num_res==2
                        nomiss_resind=[0,1,1];
                      elseif numel(miss_res23)==1 && num_res==2
                        nomiss_resind=ones(1,3);
                        nomiss_resind(miss_res23+1)=0;
                      elseif numel(miss_res23)==1 && num_res==1      
                        nomiss_resind=ones(1,3);
                        nomiss_resind(miss_res23+1)=0;
                        nomiss_resind(1)=0;
                      elseif numel(miss_res23)==2 && num_res==1      
                        nomiss_resind=ones(1,3);
                        nomiss_resind(miss_res23+1)=0;
                      elseif num_res==0;
                         nomiss_resind=zeros(1,3);
                      else
                     end

                     if num_res==0
                         nomiss_resind=zeros(1,3);
                     end


             res_ind=find(nomiss_resind);
             resall(t,res_ind)=cellfun(@str2num,rel_trial(ind_res))';
              rtall(t,res_ind)=time_trial(ind_res)-repmat(time_trial(test_onset),numel(res_ind),1);

              elseif num_res>3 %for now keep the first 3 only
                 resall(t,:)=cellfun(@str2num,rel_trial(ind_res(1:3)))'; 
                 rtall(t,:)=time_trial(ind_res(1:3))-repmat(time_trial(test_onset),3,1);
             end
         end
    end
    
    
    clear res_trial res_ind time_trial ind_res block_def borders cue_str cue_type enter_ind fb_ind miss_res23 nomiss_resind num_res rand_vec_str trial_ind v wrong right t test_onset
     
      
    
    % recode rand vec
for r=1:6
    [i,j]=ind2sub(size(rand_vec),find(rand_vec==r));
    tmp=sortrows([i,j],1);
    rand_vec_recoded(:,r)=tmp(:,2);
end

for r = 1:60
    item_seq(r, :) = item_sec4B(r, rand_vec_recoded(r,:));
end
        
% get info whether item_seq(:,4:6) are new types or similar lure
%item_sec4B has shape 410 207 610 108 602 402
    type_seq=floor(item_seq./100);
        
    % item_seq(:,4) always a new type
    for t_item=1:3
        check_lure(:,t_item)=(type_seq(:,t_item+3)==type_seq(:,1)|type_seq(:,t_item+3)==type_seq(:,2) |type_seq(:,t_item+3)==type_seq(:,3));        
        for tr=1:size(item_sec4B,1)
            if sum(type_seq(tr,t_item+3)==type_seq(tr,1:3))
                lure_pos(tr,t_item)=find(type_seq(tr,t_item+3)==type_seq(tr,1:3));
            else
                lure_pos(tr,t_item)=NaN;
            end
        end
    end
     



clear t_item r tmp

% get responses
ind123=find(cue_seq<4);
ind1=find(cue_seq==1);
ind2=find(cue_seq==2);
ind3=find(cue_seq==3);
indall=find(cue_seq==4);

resall=resall(indall,:);
res1=res123(ind1)';
res2=res123(ind2)';
res3=res123(ind3)';

responses=NaN(size(item_seq,1),3);
responses(ind1,1)=res1;
responses(ind2,2)=res2;
responses(ind3,3)=res3;
responses(indall,:)=resall;

% don't forget to uncooment this later
clear  ind1 ind2 ind3   res1 res2 res3 resall

% score token correctness
corr_res=double(responses==rand_vec_recoded(:,1:3));
corr_res(isnan(responses))=NaN;

% score type correct
% first get type for lure response
%%important > lure_resforpos indicates from item_seq 4:6, to which of the
%%presented items the lure corresponds to 
lure_response=rand_vec_recoded(:,4:6);
for t=1:size(rand_vec,1)
    tmp_ind=find(~isnan(lure_pos(t,:)));
    lure_resforpos(t,lure_pos(t,tmp_ind))=lure_response(t,tmp_ind);
end

clear tmp_ind
corr_type=double(responses==rand_vec_recoded(:,1:3)| responses==lure_resforpos);
corr_type(isnan(responses))=NaN;

% score correct all token
corr_alltoken=nanmean(corr_res,2)==1;

% score correct all type
corr_alltype=nanmean(corr_type,2)==1;

% RTs
RTs_all=zeros(size(rand_vec,1),3);
RTs_all(ind123,:)=rt123(ind123,:);
RTs_all(indall,:)=rtall(indall,:);

clear pic_ind i j indall ind123 rel_trial t tr loginfo loginfo_res loginfo_time
% which info is needed: 
        % 1: block_type 2=WM
        % 2: WM block
        % 3: trialnumber
        % 4:9 enco item & times
        % 10:11cue & times
        % 12:15: test item token + time
        % 16:23 which response codes for which item 
        % (randvec: position of x in randvec codes for item: ie randvec(4)=2 means response 4 means enco item 2) 
        % 24-29 recode randvec (column 24: correct response for item 1, 25:correct response for item 2 etc...)
        % 30-32: responses 
        % 33-35: lure pos
        % 36-38. all correct (token correct)
        % 39-41: correct on type level
        % 42: completely correct all token
        % 43: completely correct all types
        % 44-46: RT for response relative to test onset
        
        %randvec indices indicate the order of appearence in the sequence
        %randvecrecoded indicates the correct order of button presses
trlinfo=[ones(size(cue_seq)).*2,ones(size(cue_seq)).*block_num,(1:numel(cue_seq))',...
    item_seq(:,1:3),enco_time,cue_seq,cue_time',item_seq(:,4:6),test_time'...
    rand_vec,rand_vec_recoded,responses,lure_pos,corr_res,corr_type...
    corr_alltoken,corr_alltype,RTs_all];

 trlinfo(:, 48) = m1_time';
 trlinfo(:, 49) = m2_time';
 trlinfo(:, 50) = fix_time(:,1); 
 trlinfo(:, 51) = fix_time(:,2); 



    case 'practice'
    
        % which info is needed: 
        % 1: block_type 0 practice
        % 2: pic
        % 3: cue
        % 4: triggertime
        % 5: testtime
        

         logs=loginfo;
        clear loginfo;
         loginfo=logs{1,3};
         loginfo_res=logs{1,2};
         loginfo_time=logs{1,4};
         % cue type
         cue_ind=find(strncmp('cue', loginfo,3));
         pic_ind=find(strncmp('prapic',  loginfo,3)&cellfun(@numel,loginfo)<20);
         
         test_ind=find(cellfun(@numel,loginfo)>20);
         test_time_tmp=loginfo_time(test_ind);
         all_trig_ind= sort([cue_ind;pic_ind]);
         all_trig_info=loginfo(all_trig_ind);
         cue=(strncmp('cue', all_trig_info,3));
         pic=(strncmp('prapic',  all_trig_info,3));
         all_trig_time=loginfo_time(all_trig_ind);
         test_time=zeros(size(all_trig_time));
         test_time(cue)=test_time_tmp;
         
         trlinfo=[zeros(size(all_trig_time)), pic, cue, all_trig_time, test_time]; 
        
end




end











