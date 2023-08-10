function [x]=log2mat(file)

[trial, event_type,code,  touch_x,touch_y,time] = textread(file, '%*s%d%s%s%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s', 'headerlines', 5,'delimiter','\t','whitespace','');%---Liest das Logfile ein;
x{:,1}=trial;
x{:,2}=event_type;
x{:,3}=code;
x{:,4}=time; 
x{:,5}=[touch_x,touch_y];
