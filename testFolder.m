clear all;

addpath('baseline', genpath('workspace'), genpath('PEAQ'));

folder = 'audio/Speech/';

ODG = [];

files = dir([folder '*.wav']);
for i = 1:size(files,1)
    filename = files(i).name(1:end-4);
    if (filename(end-3:end) == '_48k')
        continue;
    end
    ODG = [ODG PEAQTest(filename, folder)];
end

disp(mean(ODG));