
function [x] = averSub2 (c, region)

disp('averaging sessions and subjects');
    switch region
        case 'vvs'
            x(1,:) = mean([c{1}], 'omitnan');
            x(2,:) = mean([c{2:3}], 'omitnan');
            x(3,:) = mean([c{4:5}], 'omitnan');
            x(4,:) = mean([c{6:7}], 'omitnan');
            x(5,:) = mean([c{8:9}], 'omitnan');
            x(6,:) = mean([c{10:14}], 'omitnan');
            x(7,:) = mean([c{15:18}], 'omitnan');
            x(8,:) = mean([c{19:20}], 'omitnan');
            x(9,:) = mean([c{21:22}], 'omitnan');
            x(10,:) = mean([c{23:24}], 'omitnan');
            x(11,:) = mean([c{25:26}], 'omitnan');
            x(12,:) = mean([c{27:28}], 'omitnan');
            x(13,:) = mean([c{29:30}], 'omitnan');
            x(14,:) = mean([c{31:35}], 'omitnan');
            x(15,:) = mean([c{36:37}], 'omitnan');
            x(16,:) = mean([c{38:41}], 'omitnan');
            x(17,:) = mean([c{42:43}], 'omitnan');
            x(18,:) = mean([c{44:45}], 'omitnan');
            x(19,:) = mean([c{46:47}], 'omitnan');
            x(20,:) = mean([c{48:49}], 'omitnan');
            x(21,:) = mean([c{50:51}], 'omitnan');
            x(22,:) = mean([c{52:53}], 'omitnan');
            x(23,:) = mean([c{54:55}], 'omitnan');
            x(24,:) = mean([c{56:59}], 'omitnan');
            x(25,:) = mean([c{60:62}], 'omitnan');
            x(26,:) = mean([c{63:66}], 'omitnan');
            x(27,:) = mean([c{67:70}], 'omitnan');
            x(28,:) = mean([c{71:74}], 'omitnan');
            
      case 'pfc'
            x(1,:) = mean([c{1:5}], 'omitnan');
            x(2,:) = mean([c{6:9}], 'omitnan');
            x(3,:) = mean([c{10:11}], 'omitnan');
            x(4,:) = mean([c{12:13}], 'omitnan');
            x(5,:) = mean([c{14:15}], 'omitnan');
            x(6,:) = mean([c{16:18}], 'omitnan');
            x(7,:) = mean([c{19:20}], 'omitnan');
            x(8,:) = mean([c{21:22}], 'omitnan');
            x(9,:) = mean([c{23:24}], 'omitnan');
            x(10,:) = mean([c{25:26}], 'omitnan');
            x(11,:) = mean([c{27:28}], 'omitnan');
            x(12,:) = mean([c{29:30}], 'omitnan');
            x(13,:) = mean([c{31:32}], 'omitnan');
            x(14,:) = mean([c{33:34}], 'omitnan');
            x(15,:) = mean([c{35:38}], 'omitnan');
            x(16,:) = mean([c{39:42}], 'omitnan');

            
        case 'hipp'
            eval ([     'x{1}  = '  c{i}  '{1};' ])
            eval ([     'x{2}  = [' c{i} '{2};'  c{i} '{3}];' ])
            eval ([     'x{3}  = [' c{i} '{4};'  c{i} '{5}];' ])
            eval ([     'x{4}  = [' c{i} '{6};'  c{i} '{7}];' ])
            eval ([     'x{5}  = [' c{i} '{8};'  c{i} '{9}];' ])
            eval ([     'x{6}  = [' c{i} '{10};' c{i} '{11};' c{i} '{12};' c{i} '{13};' c{i} '{14}];' ])
            eval ([     'x{7}  = [' c{i} '{15};' c{i} '{16}];'])
            eval ([     'x{8}  = [' c{i} '{17};' c{i} '{18}];'])
            eval ([     'x{9}  = [' c{i} '{19};' c{i} '{20}];' ])
            eval ([     'x{10} = [' c{i} '{21};' c{i} '{22}];' ])
            eval ([     'x{11} = [' c{i} '{23};' c{i} '{24}];' ])
            eval ([     'x{12} = [' c{i} '{25};' c{i} '{26};' c{i} '{27};' c{i} '{28}];'])
            eval ([     'x{13} = [' c{i} '{29};' c{i} '{30};'  c{i} '{31}];'])
            eval ([     'x{14} = [' c{i} '{32};' c{i} '{33};' c{i} '{34};' c{i} '{35}];'])
            eval ([     'x{15} = [' c{i} '{36};' c{i} '{37};' c{i} '{38};' c{i} '{39}];'])
            eval ([     'x{16} = [' c{i} '{40};' c{i} '{41};' c{i} '{42};' c{i} '{43}];'])
            
           
    end


end



%%















