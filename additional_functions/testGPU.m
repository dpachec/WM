%%
fprintf('Benchmarking GPU-accelerated Cross-Correlation.\n');
if ~(parallel.gpu.GPUDevice.isAvailable)
    fprintf(['\n\t**GPU not available. Stopping.**\n']);
    return;
else
    dev = gpuDevice;
    fprintf(...
    'GPU detected (%s, %d multiprocessors, Compute Capability %s)',...
    dev.Name, dev.MultiprocessorCount, dev.ComputeCapability);
end

%%
fprintf('\n\n *** Benchmarking vector-vector cross-correlation*** \n\n');
fprintf('Benchmarking function :\n');
%type('benchXcorrVec');
fprintf('\n\n');

sizes = [2000 1e4 1e5 5e5 1e6];
tc = zeros(1,numel(sizes));
tg = zeros(1,numel(sizes));
numruns = 10;

for s=1:numel(sizes);
    fprintf('Running xcorr of %d elements...\n', sizes(s));
    delchar = repmat('\b', 1,numruns);

    a = rand(sizes(s),1);
    b = rand(sizes(s),1);
    tc(s) = benchXcorrVec(a, b, numruns);
    fprintf([delchar '\t\tCPU  time : %.2f ms\n'], 1000*tc(s));
    tg(s) = benchXcorrVec(gpuArray(a), gpuArray(b), numruns);
    fprintf([delchar '\t\tGPU time :  %.2f ms\n'], 1000*tg(s));
end

%Plot the results
fig = figure;
ax = axes('parent', fig);
semilogx(ax, sizes, tc./tg, 'r*-');
ylabel(ax, 'Speedup');
xlabel(ax, 'Vector size');
title(ax, 'GPU Acceleration of XCORR');
drawnow;



%%
function t = benchXcorrVec(u,v, numruns)
%Used to benchmark xcorr with vector inputs on the CPU and GPU.
    
%   Copyright 2012 The MathWorks, Inc.

    timevec = zeros(1,numruns);
    gdev = gpuDevice;
    for ii=1:numruns
        ts = tic;
        %o = xcorr(u,v); %#ok<NASGU>
        o = corr(u,v, 'Type', 's'); %#ok<NASGU>
        wait(gdev)
        timevec(ii) = toc(ts);
        fprintf('.');
    end
    t = min(timevec);
end