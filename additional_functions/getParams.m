function [cfg] = getParams(f2sav)
    
    f2t = strsplit(f2sav, '_'); 

    cfg.brainROI    = f2t{1}; 
    lays = strsplit((strrep(f2t{2}, '-', ':')), ':')'; 
    cfg.lays2load   = double(string(lays));
    freqs = strsplit(f2t{3}, '-');
    cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    cfg.fR          = double(string((f2t{4}))); %frequency resolved
    cfg.fitMode     = double(string((f2t{5}))); %1= trials; 0 = no trials; 
    cfg.timeRes     = double(string((f2t{6})));
    if isnan(cfg.timeRes) cfg.timeRes = 'all'; end
    cfg.win_width   = double(string((f2t{7})));
    cfg.mf          = double(string((f2t{8})));

end