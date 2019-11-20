clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
ABfolders = {folders(strncmpi({folders.name}, 'AB', 2)).name};

%%
data = struct;
for i = 1:length(ABfolders)
    folders = dir([parent '/' ABfolders{i}]);
    noise_condition = {folders(strncmpi({folders.name}, 'AB', 2)).name};
    for ii = 1:length(noise_condition)
        respath = [parent '/' ABfolders{i} '/' noise_condition{ii} '/results/'];
        all_results = dir(fullfile(respath,'*.mat'));
        for iii = 1:length(all_results)
            disp(all_results(iii).name)
            load([respath all_results(iii).name]);
            data(iii).noise = vars.noise;
            data(iii).ctrl = vars.ctrl;
            data(iii).daq = daq;
            data(iii).pvd = pvd;
            data(iii).vars = vars;
            data(iii).results = calcMetrics(vars);
        end
        save([parent '/all data/' noise_condition{ii} ' data'], 'data');
        disp(['saved ' ABfolders{i}])
        clear data
    end
end