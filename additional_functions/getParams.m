function [cfg] = getParams(f2sav)
    
    f2t = strsplit(f2sav, '_'); 
    cfg.net2load     = f2t{1}; 
    cfg.brainROI    = f2t{2}; 
    cfg.period  = f2t{3};
    %lays = strsplit((strrep(f2t{2}, '-', ':')), ':')'; 
    lays = eval(strrep(f2t{4}, '-', ':'))
    cfg.lays2load   = double(string(lays));
    freqs = strsplit(f2t{5}, '-');
    cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    cfg.avRep       = double(string((f2t{6}))); %average repetitions
    cfg.avTFV       = double(string((f2t{7}))); %average repetitions
    cfg.fR          = double(string((f2t{8}))); %frequency resolved
    cfg.fitMode     = double(string((f2t{9}))); %1= trials; 0 = no trials; 
    cfg.timeRes     = double(string((f2t{10})));
    if isnan(cfg.timeRes) cfg.timeRes = 'all'; end
    cfg.win_width   = double(string((f2t{11})));
    mf = strsplit(f2t{12}, '.');
    cfg.mf          = double(string((mf{1})));

end