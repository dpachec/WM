
function [out_c ] = averageSub_WM (c, d, contrData, idData, region, noAv)

disp('averaging sessions and subjects');

if ~noAv

    for i = 1:length(contrData) 
        eval([c{i} ' = contrData{i};']);
        %eval([d{i} ' = idData{i};']);

        switch region
            case 'vvs'
                eval ([     'x{1,:} = ' c{i} '{1};' ])
                eval ([     'x{2,:} = [' c{i} '{2};' c{i} '{3}];' ])
                eval ([     'x{3,:} = [' c{i} '{4};' c{i} '{5}];' ])
                eval ([     'x{4,:} = [' c{i} '{6};' c{i} '{7}];' ])
                eval ([     'x{5,:} = [' c{i} '{8};' c{i} '{9}];' ])
                eval ([     'x{6,:} = [' c{i} '{10};' c{i} '{11};' c{i} '{12};' c{i} '{13};' c{i} '{14}];' ])
                eval ([     'x{7,:} = [' c{i} '{15};' c{i} '{16};' c{i} '{17};' c{i} '{18}];'])
                eval ([     'x{8,:} = [' c{i} '{19};' c{i} '{20}];' ])
                eval ([     'x{9,:} = [' c{i} '{21};' c{i} '{22}];' ])
                eval ([     'x{10,:} = [' c{i} '{23};' c{i} '{24}];' ])
                eval ([     'x{11,:} = [' c{i} '{25};' c{i} '{26}];' ])
                eval ([     'x{12,:} = [' c{i} '{27};' c{i} '{28}];' ])
                eval ([     'x{13,:} = [' c{i} '{29};' c{i} '{30}];' ])
                eval ([     'x{14,:} = [' c{i} '{31};' c{i} '{32};' c{i} '{33};' c{i} '{34};' c{i} '{35}];' ])
                eval ([     'x{15,:} = [' c{i} '{36};' c{i} '{37}];' ])
                eval ([     'x{16,:} = [' c{i} '{38};' c{i} '{39};' c{i} '{40};' c{i} '{41}];' ])
                eval ([     'x{17,:} = [' c{i} '{42};' c{i} '{43}];' ])
                eval ([     'x{18,:} = [' c{i} '{44};' c{i} '{45}];' ])
                eval ([     'x{19,:} = [' c{i} '{46};' c{i} '{47}];' ])
                eval ([     'x{20,:} = [' c{i} '{48};' c{i} '{49}];'])
                eval ([     'x{21,:} = [' c{i} '{50};' c{i} '{51}];'])
                eval ([     'x{22,:} = [' c{i} '{52};' c{i} '{53}];'])
                eval ([     'x{23,:} = [' c{i} '{54};' c{i} '{55}];'])
                eval ([     'x{24,:} = [' c{i} '{56};' c{i} '{57};' c{i} '{58};' c{i} '{59}];'])
                eval ([     'x{25,:} = [' c{i} '{60};' c{i} '{61};' c{i} '{62}];' ])
                eval ([     'x{26,:} = [' c{i} '{63};' c{i} '{64};' c{i} '{65};' c{i} '{66}];'])
                eval ([     'x{27,:} = [' c{i} '{67};' c{i} '{68};' c{i} '{69};' c{i} '{70}];'])
                eval ([     'x{28,:} = [' c{i} '{71};' c{i} '{72};' c{i} '{73};' c{i} '{74}];'])
                
%                 eval ([     'y{1,:} = ' d{i} '{1};' ])
%                 eval ([     'y{2,:} = [' d{i} '{2};' d{i} '{3}];' ])
%                 eval ([     'y{3,:} = [' d{i} '{4};' d{i} '{5}];' ])
%                 eval ([     'y{4,:} = [' d{i} '{6};' d{i} '{7}];' ])
%                 eval ([     'y{5,:} = [' d{i} '{8};' d{i} '{9}];' ])
%                 eval ([     'y{6,:} = [' d{i} '{10};' d{i} '{11};' d{i} '{12};' d{i} '{13};' d{i} '{14}];' ])
%                 eval ([     'y{7,:} = [' d{i} '{15};' d{i} '{16};' d{i} '{17};' d{i} '{18}];'])
%                 eval ([     'y{8,:} = [' d{i} '{19};' d{i} '{20}];' ])
%                 eval ([     'y{9,:} = [' d{i} '{21};' d{i} '{22}];' ])
%                 eval ([     'y{10,:} = [' d{i} '{23};' d{i} '{24}];' ])
%                 eval ([     'y{11,:} = [' d{i} '{25};' d{i} '{26}];' ])
%                 eval ([     'y{12,:} = [' d{i} '{27};' d{i} '{28}];' ])
%                 eval ([     'y{13,:} = [' d{i} '{29};' d{i} '{30}];' ])
%                 eval ([     'y{14,:} = [' d{i} '{31};' d{i} '{32};' d{i} '{33};' d{i} '{34};' d{i} '{35}];' ])
%                 eval ([     'y{15,:} = [' d{i} '{36};' d{i} '{37}];' ])
%                 eval ([     'y{16,:} = [' d{i} '{38};' d{i} '{39};' d{i} '{40};' d{i} '{41}];' ])
%                 eval ([     'y{17,:} = [' d{i} '{42};' d{i} '{43}];' ])
%                 eval ([     'y{18,:} = [' d{i} '{44};' d{i} '{45}];' ])
%                 eval ([     'y{19,:} = [' d{i} '{46};' d{i} '{47}];' ])
%                 eval ([     'y{20,:} = [' d{i} '{48};' d{i} '{49}];'])
%                 eval ([     'y{21,:} = [' d{i} '{50};' d{i} '{51}];'])
%                 eval ([     'y{22,:} = [' d{i} '{52};' d{i} '{53}];'])
%                 eval ([     'y{23,:} = [' d{i} '{54};' d{i} '{55}];'])
%                 eval ([     'y{24,:} = [' d{i} '{56};' d{i} '{57};' d{i} '{58};' d{i} '{59}];'])
%                 eval ([     'y{25,:} = [' d{i} '{60};' d{i} '{61};' d{i} '{62}];' ])
%                 eval ([     'y{26,:} = [' d{i} '{63};' d{i} '{64};' d{i} '{65};' d{i} '{66}];'])
%                 eval ([     'y{27,:} = [' d{i} '{67};' d{i} '{68};' d{i} '{69};' d{i} '{70}];'])
%                 eval ([     'y{28,:} = [' d{i} '{71};' d{i} '{72};' d{i} '{73};' d{i} '{74}];'])


          case 'pfc'
                eval ([     'x{1} = [' c{i} '{1};' c{i} '{2};' c{i} '{3};' c{i} '{4};' c{i} '{5}];' ])
                eval ([     'x{2} = [' c{i} '{6};' c{i} '{7};' c{i} '{8};' c{i} '{9}];'])
                eval ([     'x{3} = [' c{i} '{10};' c{i} '{11}];' ])
                eval ([     'x{4} = [' c{i} '{12};' c{i} '{13}];' ])
                eval ([     'x{5} = [' c{i} '{14};' c{i} '{15}];' ])
                eval ([     'x{6} = [' c{i} '{16};' c{i} '{17};' c{i} '{18}];'  ])
                eval ([     'x{7} = [' c{i} '{19};' c{i} '{20}];' ])
                eval ([     'x{8} = [' c{i} '{21};' c{i} '{22}];' ])
                eval ([     'x{9} = [' c{i} '{23};' c{i} '{24}];' ])
                eval ([     'x{10} = [' c{i} '{25};' c{i} '{26}];' ])
                eval ([     'x{11} = [' c{i} '{27};' c{i} '{28}];' ])
                eval ([     'x{12} = [' c{i} '{29};' c{i} '{30}];' ])
                eval ([     'x{13} = [' c{i} '{31};' c{i} '{32}];' ])
                eval ([     'x{14} = [' c{i} '{33};' c{i} '{34}];' ])
                eval ([     'x{15} = [' c{i} '{35};' c{i} '{36};' c{i} '{37};' c{i} '{38}];'])
                eval ([     'x{16} = [' c{i} '{39};' c{i} '{40};' c{i} '{41};' c{i} '{42}];'])
                
%                 eval ([     'y{1} = [' d{i} '{1};' d{i} '{2};' d{i} '{3};' d{i} '{4};' d{i} '{5}];' ])
%                 eval ([     'y{2} = [' d{i} '{6};' d{i} '{7};' d{i} '{8};' d{i} '{9}];'])
%                 eval ([     'y{3} = [' d{i} '{10};' d{i} '{11}];' ])
%                 eval ([     'y{4} = [' d{i} '{12};' d{i} '{13}];' ])
%                 eval ([     'y{5} = [' d{i} '{14};' d{i} '{15}];' ])
%                 eval ([     'y{6} = [' d{i} '{16};' d{i} '{17};' d{i} '{18}];'  ])
%                 eval ([     'y{7} = [' d{i} '{19};' d{i} '{20}];' ])
%                 eval ([     'y{8} = [' d{i} '{21};' d{i} '{22}];' ])
%                 eval ([     'y{9} = [' d{i} '{23};' d{i} '{24}];' ])
%                 eval ([     'y{10} = [' d{i} '{25};' d{i} '{26}];' ])
%                 eval ([     'y{11} = [' d{i} '{27};' d{i} '{28}];' ])
%                 eval ([     'y{12} = [' d{i} '{29};' d{i} '{30}];' ])
%                 eval ([     'y{13} = [' d{i} '{31};' d{i} '{32}];' ])
%                 eval ([     'y{14} = [' d{i} '{33};' d{i} '{34}];' ])
%                 eval ([     'y{15} = [' d{i} '{35};' d{i} '{36};' d{i} '{37};' d{i} '{38}];'])
%                 eval ([     'y{16} = [' d{i} '{39};' d{i} '{40};' d{i} '{41};' d{i} '{42}];'])
                
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
                
                eval ([     'y{1}  = '  d{i}  '{1};' ])
                eval ([     'y{2}  = [' d{i} '{2};'  d{i} '{3}];' ])
                eval ([     'y{3}  = [' d{i} '{4};'  d{i} '{5}];' ])
                eval ([     'y{4}  = [' d{i} '{6};'  d{i} '{7}];' ])
                eval ([     'y{5}  = [' d{i} '{8};'  d{i} '{9}];' ])
                eval ([     'y{6}  = [' d{i} '{10};' d{i} '{11};' d{i} '{12};' d{i} '{13};' d{i} '{14}];' ])
                eval ([     'y{7}  = [' d{i} '{15};' d{i} '{16}];'])
                eval ([     'y{8}  = [' d{i} '{17};' d{i} '{18}];'])
                eval ([     'y{9}  = [' d{i} '{19};' d{i} '{20}];' ])
                eval ([     'y{10} = [' d{i} '{21};' d{i} '{22}];' ])
                eval ([     'y{11} = [' d{i} '{23};' d{i} '{24}];' ])
                eval ([     'y{12} = [' d{i} '{25};' d{i} '{26};' d{i} '{27};' d{i} '{28}];'])
                eval ([     'y{13} = [' d{i} '{29};' d{i} '{30};'  d{i} '{31}];'])
                eval ([     'y{14} = [' d{i} '{32};' d{i} '{33};' d{i} '{34};' d{i} '{35}];'])
                eval ([     'y{15} = [' d{i} '{36};' d{i} '{37};' d{i} '{38};' d{i} '{39}];'])
                eval ([     'y{16} = [' d{i} '{40};' d{i} '{41};' d{i} '{42};' d{i} '{43}];'])

        end


        eval ([ 'clear ' c{i} ]);
        eval ([ 'c{i} = x;']);
        
        eval ([ 'clear ' d{i} ]);
        %eval ([ 'd{i} = y;']);

    end


    for i = 1:length(contrData) 

        % % average sessions for each subject
        out_c{i} = c{i};
        %out_id{i} = d{i};

        % % % no average
        %out_c{i} = contrData{i};

    end

else


    for i = 1:length(contrData) 

     out_c{i} = contrData{i};
     %out_id{i} = idData{i};

    end



end



%%














