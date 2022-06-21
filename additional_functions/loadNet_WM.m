f2t = strsplit(f2sav, '_'); 
if strcmp(f2t{1}, 'Alex')
    if ~exist('act_CH')  | size(act_CH, 1) ~= 8
        ver = 'imageNet'; % 'imageNet' or 'ecoset'
        [act_CH act_FR] = load_alex_activ(ver);
    end
elseif strcmp(f2t{1}, 'RNN')
    if ~exist('act_CH')  | size(act_CH, 1) ~= 56
        [act_CH act_FR] = load_rnn_activ;
    end
else
    act_CH = []; 
    act_FR = [];
end

disp('network loaded')