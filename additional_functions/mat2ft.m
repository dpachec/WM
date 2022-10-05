function ft_data = mat2ft(data, fs, labels, offset)
%MAT2FT transform a matrix into fieldtrip data structure 
%   ft_data = MAT2FT(X, fs) transforms the data in X (samples x channels x
%   trials) into a fieldtrip compatible data structure
%   
%   ft_data = MAT2FT(X, fs, varargin)
%       accepts optional arguments:
%       'labels' - 1xChannels cell array of channels labels
%       'offset' - pre-stimulus offset, in seconds
%
%   shenin@gc.cuny.edu
%   simon.henin@nyumc.org
%   (c) September, 25th, 2017
%   Sept. 25th 2007 - version 1.0
%
if ~exist('labels', 'var')|isempty(labels),
    labels = cellstr(num2str([1:size(data,1)]'))';
end
if ~exist('offset')|isempty(offset),
    offset = 0;
end

ft_data = struct('trial', [], 'time', [], 'label', [], 'fsample', []);

t = (0:size(data,2)-1)./fs-offset;
for i=1:size(data,3),
   ft_data.trial{end+1} = data(:, :,i);
   ft_data.time{end+1} = t;
end
ft_data.label = {labels{:}};
ft_data.fsample = fs;
ft_data.sampleinfo= [];

