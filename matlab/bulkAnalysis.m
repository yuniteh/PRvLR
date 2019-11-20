clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
testfolders = {folders(strncmpi({folders.name}, 'test', 4)).name};

%%
total.noise = zeros(length(testfolders),1);
total.complete_rate = total.noise;
total.complete_time = total.noise;
total.path_eff = total.noise;

for i = 1:length(testfolders)
    file = char(testfolders(i));
    respath = strcat(parent,'\results\');
    resfile = strcat(respath, file, '.mat');
    if exist(resfile, 'file') == 0
        disp(['getting ' file ' results...'])
        [vars, daq, pvd, results] = loadFile(strcat(parent,'\',file));
        %% save results
        if exist(respath,'dir') == 0
            mkdir(respath);
        end
        
        save(resfile, 'vars', 'daq', 'pvd', 'results');
    else
        disp(['loading ' file ' results...'])
        load(resfile);
    end
    %%
    total.noise(i) = vars.noise;
    total.complete_rate(i) = results.complete_rate;
    total.complete_time(i) = mean(results.complete_time);
    total.path_eff(i) = mean(results.path_eff);
end

%% all plots
c = colormap(lines);
total.noise = total.noise/100;

ind = strfind(parent,'\');
bulk = parent(ind(end) + 1:end);
[total.noise, noise_ind] = sort(total.noise);

% completion rate
subplot(311)
plot(total.noise, total.complete_rate(noise_ind), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
ylabel('Completion Rate')
axis tight
ylim([0 1])
title([bulk ' Results'])

% completion time
subplot(312)
plot(total.noise, total.complete_time(noise_ind), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Completion Time (s)')
axis tight
ylim([0 10])

% path efficiency
subplot(313)
plot(total.noise, total.path_eff(noise_ind), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(3,:))
ylabel('Path Efficiency')
xlabel('Noise (V)')
axis tight
ylim([0 1])
