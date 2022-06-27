function [cfg] = getParams(f2sav)
    
    f2t = strsplit(f2sav, '_'); 
    cfg.net2load     = f2t{1}; 
    cfg.brainROI    = f2t{2}; 
    %lays = strsplit((strrep(f2t{2}, '-', ':')), ':')'; 
    lays = eval(strrep(f2t{3}, '-', ':'))
    cfg.lays2load   = double(string(lays));
    freqs = strsplit(f2t{4}, '-');
    cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    cfg.avRep       = double(string((f2t{5}))); %average repetitions
    cfg.avTFV       = double(string((f2t{6}))); %average repetitions
    cfg.fR          = double(string((f2t{7}))); %frequency resolved
    cfg.fitMode     = double(string((f2t{8}))); %1= trials; 0 = no trials; 
    cfg.timeRes     = double(string((f2t{9})));
    if isnan(cfg.timeRes) cfg.timeRes = 'all'; end
    cfg.win_width   = double(string((f2t{10})));
    mf = strsplit(f2t{11}, '.');
    cfg.mf          = double(string((mf{1})));

end