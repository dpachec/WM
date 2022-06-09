function[files] = getFilesWM(path)


files = dir([path]); cbytes = [files.bytes]; 
files = files(cbytes > 2000000);files = sort({files.name})'; %gets rid of console_log

end