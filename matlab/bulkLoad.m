clear
close all
clc

%% load data
parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
trainfolders = {folders(strncmpi({folders.name}, 'train', 5)).name};
freefolders = {folders(strncmpi({folders.name}, 'free', 4)).name};
testfolders = {folders(strncmpi({folders.name}, 'test', 4)).name};
allfolders = [trainfolders testfolders freefolders];

%%
for i = 1:length(allfolders)
    file = char(allfolders(i));
    respath = strcat(parent,'\results\');
    resfile = strcat(respath, file, '.mat');
    if exist(resfile, 'file') == 0
        disp(['calculating ' file ' results...'])
        [vars, daq, pvd, results] = loadFile(strcat(parent,'\',file));
        %% save results
        if exist(respath,'dir') == 0
            mkdir(respath);
        end
        
        save(resfile, 'vars', 'daq', 'pvd', 'results');
    end
end