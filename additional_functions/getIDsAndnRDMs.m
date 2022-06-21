function [allS ids actCH2 actFR2 cM] = getIDsAndnRDMs(rdm_prev, actCH, actFR)


% get ids
ids = rdm_prev(1).ids; % take first freq only to get ids
ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
ids0 = cellfun(@(x) strsplit(x(3), '_'), ids, 'un', 0); 
ids1 = vertcat(ids0{:}); 
ids2 = cellfun(@(x) x([1 3]), ids1, 'un', 0);
ids3 = double(string(ids2));
idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);
ids = ids4; 


cM = zeros (length(ids0));
ids0 = string(ids0);
for i = 1:length(cM)
    for j = 1:length(cM)
        if strcmp(ids0{i}(1), ids0{j}(1)) 
            cM(i, j) = 1;
        end
    end 
end    

% format neural data
if ndims(rdm_prev(1).rdm) == 3
    allS = cat(4, rdm_prev.rdm);
    allS2 = allS;allS2(allS2==1) = 1000; allS2(allS2==0) = 2000;
    allS3 = arrayfun(@(j)  (arrayfun(@(i)tril(squeeze(allS2(:,:,i,j)), -1), 1:size(allS2,3), 'un', 0)), 1:size(allS2,4), 'un', 0);
    allS4 = vertcat(allS3{:});
    allS7 = cellfun(@(x) x(x~=0), allS4, 'un', 0);     
    allS8 = cat(2, allS7 {:});
    allS8 (allS8 ==1000) = 1;allS8 (allS8 ==2000) = 0;
    allS = allS8; 
else
    allS = cat(3, rdm_prev.rdm);
    allS2 = allS;allS2(allS2==1) = 1000; allS2(allS2==0) = 2000;
    allS3 = arrayfun(@(i)tril(squeeze(allS2(:,:,i)), -1), 1:size(allS2,3), 'un', 0);
    allS7 = cellfun(@(x) x(x~=0), allS3, 'un', 0);     
    allS8 = cat(2, allS7 {:});
    allS8 (allS8 ==1000) = 1;allS8 (allS8 ==2000) = 0;
    allS = allS8; 
    
end



% build average matrices if 'all' analysis
if size(ids4, 2) == 3
    
    %mean in DNN 
    all3CH(:, :, :, 1) = actCH(:, ids4(:,1), ids4(:,1)); 
    all3CH(:, :, :, 2) = actCH(:, ids4(:,2), ids4(:,2)); 
    all3CH(:, :, :, 3) = actCH(:, ids4(:,3), ids4(:,3)); 

    all3FR(:, :, :, 1) = actFR(:, ids4(:,1), ids4(:,1)); 
    all3FR(:, :, :, 2) = actFR(:, ids4(:,2), ids4(:,2)); 
    all3FR(:, :, :, 3) = actFR(:, ids4(:,3), ids4(:,3));         

    actCH2 = mean(all3CH, 4);
    actFR2 = mean(all3FR, 4);

    % just one (as category model) if one of 3 is present (not useful)
%     M = zeros (size(ids4, 1));
%     for i = 1:length(M)
%         for j = 1:length(M)
%             junts = [ids1(i,:)  ids1(j, :)];
%             if length(unique(junts)) < 6 %at least one repetition
%                 M(i, j) = 1;
%             end
%         end 
%     end    
%     actCH2 = repmat(M, [1, 1, 56]);
%     actCH2 = permute(actCH2, [3 2 1]);
%     actFR2 = actCH2;

else
    
    actCH2 = actCH(:, ids4, ids4); 
    actFR2 = actFR(:, ids4, ids4); 
end

    
   