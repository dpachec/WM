function[files] = getFiles(path)

files = dir([path]); cbytes = [files.bytes]; 
files = files(cbytes > 10000);files = sort({files.name})';

end