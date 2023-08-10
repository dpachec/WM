function [gM1] = average_xGM (gM, region)

clear gM1
    
switch region
    
    case 'vvs'
        gM1(1,:,:,:) = gM(1,:,:,:);
        gM1(2,:,:,:) = mean(gM(2:3,:,:,:), 'omitnan');
        gM1(3,:,:,:) = mean(gM(4:5,:,:,:), 'omitnan');
        gM1(4,:,:,:) = mean(gM(6:7,:,:,:), 'omitnan');
        gM1(5,:,:,:) = mean(gM(8:9,:,:,:), 'omitnan');
        gM1(6,:,:,:) = mean(gM(10:14,:,:,:), 'omitnan');
        gM1(7,:,:,:) = mean(gM(15:18,:,:,:), 'omitnan');
        gM1(8,:,:,:) = mean(gM(19:20,:,:,:), 'omitnan');
        gM1(9,:,:,:) = mean(gM(21:22,:,:,:), 'omitnan');
        gM1(10,:,:,:) = mean(gM(23:24,:,:,:), 'omitnan');
        gM1(11,:,:,:) = mean(gM(25:26,:,:,:), 'omitnan');
        gM1(12,:,:,:) = mean(gM(27:28,:,:,:), 'omitnan');
        gM1(13,:,:,:) = mean(gM(29:30,:,:,:), 'omitnan');
        gM1(14,:,:,:) = mean(gM(31:35,:,:,:), 'omitnan');
        gM1(15,:,:,:) = mean(gM(36:37,:,:,:), 'omitnan');
        gM1(16,:,:,:) = mean(gM(38:41,:,:,:), 'omitnan');
        gM1(17,:,:,:) = mean(gM(42:43,:,:,:), 'omitnan');
        gM1(18,:,:,:) = mean(gM(44:45,:,:,:), 'omitnan');
        gM1(19,:,:,:) = mean(gM(46:47,:,:,:), 'omitnan');
        gM1(20,:,:,:) = mean(gM(48:49,:,:,:), 'omitnan');
        gM1(21,:,:,:) = mean(gM(50:51,:,:,:), 'omitnan');
        gM1(22,:,:,:) = mean(gM(52:53,:,:,:), 'omitnan');
        gM1(23,:,:,:) = mean(gM(54:55,:,:,:), 'omitnan');
        gM1(24,:,:,:) = mean(gM(56:59,:,:,:), 'omitnan');
        gM1(25,:,:,:) = mean(gM(60:62,:,:,:), 'omitnan');
        gM1(26,:,:,:) = mean(gM(63:66,:,:,:), 'omitnan');
        gM1(27,:,:,:) = mean(gM(67:70,:,:,:), 'omitnan');
        gM1(28,:,:,:) = mean(gM(71:74,:,:,:), 'omitnan');
        
    case 'pfc'
        gM1(1,:,:,:) = mean(gM(1:5,:,:,:), 'omitnan');
        gM1(2,:,:,:) = mean(gM(7:9,:,:,:), 'omitnan');
        gM1(3,:,:,:) = mean(gM(10:11,:,:,:), 'omitnan');
        gM1(4,:,:,:) = mean(gM(12:13,:,:,:), 'omitnan');
        gM1(5,:,:,:) = mean(gM(14:15,:,:,:), 'omitnan');
        gM1(6,:,:,:) = mean(gM(16:18,:,:,:), 'omitnan');
        gM1(7,:,:,:) = mean(gM(19:20,:,:,:), 'omitnan');
        gM1(8,:,:,:) = mean(gM(21:22,:,:,:), 'omitnan');
        gM1(9,:,:,:) = mean(gM(23:24,:,:,:), 'omitnan');
        gM1(10,:,:,:) = mean(gM(25:26,:,:,:), 'omitnan');
        gM1(11,:,:,:) = mean(gM(27:28,:,:,:), 'omitnan');
        gM1(12,:,:,:) = mean(gM(29:30,:,:,:), 'omitnan');
        gM1(13,:,:,:) = mean(gM(31:32,:,:,:), 'omitnan');
        gM1(14,:,:,:) = mean(gM(33:34,:,:,:), 'omitnan');
        gM1(15,:,:,:) = mean(gM(35:38,:,:,:), 'omitnan');
        gM1(16,:,:,:) = mean(gM(39:42,:,:,:), 'omitnan');
        
        
        
      case 'all'
          
        gM1(1,:,:,:) = gM(1,:,:,:);
        gM1(2,:,:,:) = mean(gM(2:3,:,:,:), 'omitnan');
        gM1(3,:,:,:) = mean(gM(4:5,:,:,:), 'omitnan');
        gM1(4,:,:,:) = mean(gM(6:7,:,:,:), 'omitnan');
        gM1(5,:,:,:) = mean(gM(8:9,:,:,:), 'omitnan');
        gM1(6,:,:,:) = mean(gM(10:14,:,:,:), 'omitnan');
        gM1(7,:,:,:) = mean(gM(15:18,:,:,:), 'omitnan');
        gM1(8,:,:,:) = mean(gM(19:20,:,:,:), 'omitnan');
        gM1(9,:,:,:) = mean(gM(21:22,:,:,:), 'omitnan');
        gM1(10,:,:,:) = mean(gM(23:24,:,:,:), 'omitnan');
        gM1(11,:,:,:) = mean(gM(25:26,:,:,:), 'omitnan');
        gM1(12,:,:,:) = mean(gM(27:28,:,:,:), 'omitnan');
        gM1(13,:,:,:) = mean(gM(29:30,:,:,:), 'omitnan');
        gM1(14,:,:,:) = mean(gM(31:32,:,:,:), 'omitnan');
        gM1(15,:,:,:) = mean(gM(33:37,:,:,:), 'omitnan');
        gM1(16,:,:,:) = mean(gM(38:39,:,:,:), 'omitnan');
        gM1(17,:,:,:) = mean(gM(40:42,:,:,:), 'omitnan');
        gM1(18,:,:,:) = mean(gM(43:46,:,:,:), 'omitnan');
        gM1(19,:,:,:) = mean(gM(47:48,:,:,:), 'omitnan');
        gM1(20,:,:,:) = mean(gM(49:50,:,:,:), 'omitnan');
        gM1(21,:,:,:) = mean(gM(51:52,:,:,:), 'omitnan');
        gM1(22,:,:,:) = mean(gM(53:54,:,:,:), 'omitnan');
        gM1(23,:,:,:) = mean(gM(55:56,:,:,:), 'omitnan');
        gM1(24,:,:,:) = mean(gM(57:58,:,:,:), 'omitnan');
        gM1(25,:,:,:) = mean(gM(59:60,:,:,:), 'omitnan');
        gM1(26,:,:,:) = mean(gM(61:62,:,:,:), 'omitnan');
        gM1(27,:,:,:) = mean(gM(63:64,:,:,:), 'omitnan');
        gM1(28,:,:,:) = mean(gM(65:68,:,:,:), 'omitnan');
        gM1(29,:,:,:) = mean(gM(69:71,:,:,:), 'omitnan');
        gM1(30,:,:,:) = mean(gM(72:75,:,:,:), 'omitnan');
        gM1(31,:,:,:) = mean(gM(76:79,:,:,:), 'omitnan');
        gM1(32,:,:,:) = mean(gM(80:83,:,:,:), 'omitnan');


        
        
end


% 
% 

% 
