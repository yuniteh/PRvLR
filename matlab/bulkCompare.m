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
pr_ct = length(find(strncmpi({folders.name}, 'test_1', 6)));
lr_ct = length(find(strncmpi({folders.name}, 'test_2', 6)));

%%
group1 = input('Label 1 (default = LDA):','s');
group2 = input('Label 2 (default = LR):','s');
if isempty(group1)
    group1 = 'LDA';
end
if isempty(group2)
    group2 = 'LR';
end

%%
pr.noise = zeros(pr_ct,1);
pr.complete_rate = pr.noise;
pr.complete_time = pr.noise;
pr.path_eff = pr.noise;
pr.numTrials = pr.noise;
pr.complete_ind = cell(pr_ct,1);
pr.error = pr.complete_ind;
prTemp.complete_ind = pr.complete_ind;
prTemp.error = pr.complete_ind;
pr_i = 1;

lr.noise = zeros(lr_ct,1);
lr.complete_rate = lr.noise;
lr.complete_time = lr.noise;
lr.path_eff = lr.noise;
lr.numTrials = lr.noise;
lr.complete_ind = cell(lr_ct,1);
lr.error = lr.complete_ind;
lrTemp.complete_ind = lr.complete_ind;
lrTemp.error = lr.complete_ind;
lr_i = 1;

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
    if vars.ctrl == 1
        pr.noise(pr_i) = vars.noise;
        pr.complete_rate(pr_i) = results.complete_rate;
        pr.complete_time(pr_i) = mean(results.complete_time);
        pr.path_eff(pr_i) = mean(results.path_eff(results.complete_ind));
        pr.numTrials(pr_i) = results.numTrials;
        prTemp.complete_ind{pr_i} = results.complete_ind;
        prTemp.error{pr_i} = results.error;
        pr.in_time(pr_i) = mean(results.in_time);
        pr_i = pr_i + 1;
    elseif vars.ctrl == 2
        lr.noise(lr_i) = vars.noise;
        lr.complete_rate(lr_i) = results.complete_rate;
        lr.complete_time(lr_i) = mean(results.complete_time);
        lr.path_eff(lr_i) = mean(results.path_eff(results.complete_ind));
        lr.numTrials(lr_i) = results.numTrials;
        lrTemp.complete_ind{lr_i} = results.complete_ind;
        lrTemp.error{lr_i} = results.error;
        lr.in_time(lr_i) = mean(results.in_time);
        lr_i = lr_i + 1;
    end
    clear test daq pvd results
end

%% all plots
c = colormap(lines);
pr.noise = pr.noise/100;
lr.noise = lr.noise/100;

ind = strfind(parent,'\');
name = parent(ind(end) + 1:end);
[pr.noise, noise_pr] = sort(pr.noise);
[lr.noise, noise_lr] = sort(lr.noise);

pr.complete_rate = pr.complete_rate(noise_pr);
pr.complete_time = pr.complete_time(noise_pr);
pr.path_eff = pr.path_eff(noise_pr);
pr.numTrials = pr.numTrials(noise_pr);
pr.in_time = pr.in_time(noise_pr);

lr.complete_rate = lr.complete_rate(noise_lr);
lr.complete_time = lr.complete_time(noise_lr);
lr.path_eff = lr.path_eff(noise_lr);
lr.numTrials = lr.numTrials(noise_lr);
lr.in_time = lr.in_time(noise_lr);

for i = 1:pr_ct
    pr.complete_ind{i} = prTemp.complete_ind{noise_pr(i)};
    pr.error{i} = prTemp.error{noise_pr(i)};
    lr.complete_ind{i} = lrTemp.complete_ind{noise_lr(i)};
    lr.error{i} = lrTemp.error{noise_lr(i)};
end

% completion rate
subplot(411)
hold all
plot(pr.noise, pr.complete_rate, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.complete_rate, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Completion Rate')
axis tight
ylim([0 1])
title([name ' Results'])
legend(group1, group2)

% completion time
subplot(412)
hold all
plot(pr.noise, pr.complete_time, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.complete_time, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Completion Time (s)')
axis tight
ylim([0 10])

% path efficiency
subplot(413)
hold all
plot(pr.noise, pr.path_eff, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.path_eff, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Path Efficiency')
axis tight
ylim([0 1])

% time in target
subplot(414)
hold all
plot(pr.noise, pr.in_time, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.in_time, '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Time in Target (s)')
xlabel('Noise (V)')
axis tight
saveas(gcf, [respath name '_summary.png'])

%% euclidean distance for completed trials (LDA)
pr.c = linspecer(pr_ct + 2,'blue');
lr.c = linspecer(lr_ct + 2,'red');
pr.c = pr.c(3:end,:);
lr.c = lr.c(3:end,:);
t = 0:.1:10;
t = t(1:end - 1);
figure
for i = 1:pr_ct
    subplot(pr_ct,1,i)
    hold all
    error = pr.error{i}(pr.complete_ind{i},:);
    if ~isempty(error)
        plot(t, error, 'Color', pr.c(i,:), 'LineWidth', 1.5)
        plot(t, mean(error,1), 'k', 'LineWidth',1.5)
    end
    xlim([0 10])
    ylim([0 3])
    ylabel([num2str(pr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Completed Trials (' group1 ')'])
    end
end
xlabel('Time (s)')
saveas(gcf, [respath name '_ldaComplete.png'])
clear error
%% euclidean distance for failed trials LDA
pr.mean_failed = zeros(pr_ct,100);
figure
for i = 1:pr_ct
    subplot(pr_ct,1,i)
    hold all
    fail_ind = setdiff(1:pr.numTrials(i),pr.complete_ind{i});
    error = pr.error{i}(fail_ind,:);
    if ~isempty(error)
        plot(t, error, 'Color', pr.c(i,:), 'LineWidth', 1)
        plot(t, mean(error,1), 'k', 'LineWidth',1.5)
        pr.mean_failed(i,:) = mean(error,1);
    end
    xlim([0 10])
    ylim([0 6])
    ylabel([num2str(pr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Failed Trials (' group1 ')'])
    end
end
xlabel('Time (s)')
saveas(gcf, [respath name '_ldaFailed.png'])
clear error
%% euclidean distance for completed trials LR
figure
for i = 1:lr_ct
    subplot(lr_ct,1,i)
    hold all
    error = lr.error{i}(lr.complete_ind{i},:);
    if ~isempty(error)
        plot(t, error, 'Color', lr.c(i,:), 'LineWidth', 1.5)
        plot(t, mean(error,1), 'k', 'LineWidth',1.5)
    end
    xlim([0 10])
    ylim([0 3])
    ylabel([num2str(lr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Completed Trials (' group2 ')'])
    end
end
xlabel('Time (s)')
saveas(gcf, [respath name '_lrComplete.png'])
clear error
%%
lr.mean_failed = zeros(lr_ct,100);
figure
for i = 1:lr_ct
    subplot(lr_ct,1,i)
    hold all
    fail_ind = setdiff(1:lr.numTrials(i), lr.complete_ind{i});
    error = lr.error{i}(fail_ind,:);
    if ~isempty(error)
        plot(t, error, 'Color', lr.c(i,:), 'LineWidth', 1)
        plot(t, mean(error,1), 'k', 'LineWidth',1.5)
        lr.mean_failed(i,:) = mean(error,1);
    end
    ylabel([num2str(lr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Failed Trials (' group2 ')'])
    end
    xlim([0 10])
    ylim([0 6])
end
xlabel('Time (s)')
saveas(gcf, [respath name '_lrFailed.png'])
clear error
%% mean error for failed trials
figure
subplot(121)
hold all
for i = 1:pr_ct
    if sum(pr.mean_failed(i,:) > 0)
        plot(t, pr.mean_failed(i,:),'Color', pr.c(i,:), 'LineWidth',1.5)
    else
        plot(0,0,'Color', pr.c(i,:))
    end
end
legend(num2str(pr.noise))
title([name ' Ave Euclidean Distance for Failed Trials (' group1 ')'])
ylim([0 6])

subplot(122)
hold all
for i = 1:lr_ct
    if sum(lr.mean_failed(i,:) > 0)
        plot(t,lr.mean_failed(i,:),'Color', lr.c(i,:), 'LineWidth',1.5)
    else
        plot(0,0,'Color', lr.c(i,:))
    end
end
legend(num2str(lr.noise))
title([name ' Ave Euclidean Distance for Failed Trials (' group2 ')'])
xlabel('Time (s)')
ylim([0 6])
saveas(gcf, [respath name '_aveFailed.png'])

%% mean error for all trials
figure
subplot(121)
hold all
for i = 1:pr_ct
    plot(t, mean(pr.error{i}, 1), 'Color', pr.c(i,:), 'LineWidth', 1.5)
end
legend(num2str(pr.noise))
title([name ' Ave Euclidean Distance (' group1 ')'])
ylim([0 6])

subplot(122)
hold all
for i = 1:lr_ct
    plot(t, mean(lr.error{i}, 1), 'Color', lr.c(i,:), 'LineWidth', 1.5)
end
legend(num2str(lr.noise))
title([name ' Ave Euclidean Distance (' group2 ')'])
xlabel('Time (s)')
ylim([0 6])
saveas(gcf, [respath name '_aveAll.png'])
