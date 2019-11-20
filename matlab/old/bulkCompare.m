clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
testfolders = {};
pr_ct = 0;
lr_ct = 0;

for i = 1:length(folders)
    if ~isempty(strfind(folders(i).name,'test'))
        testfolders = [testfolders; folders(i).name];
        if ~isempty(strfind(folders(i).name,'test_1'))
            pr_ct = pr_ct + 1;
        else
            lr_ct = lr_ct + 1;
        end
        
    end
end

%%
pr.noise = zeros(pr_ct,1);
pr.complete_rate = pr.noise;
pr.complete_time = pr.noise;
pr.path_eff = pr.noise;
pr.numTrials = pr.noise;
pr.complete_ind = cell(pr_ct,1);
pr.error = pr.complete_ind;
pr_i = 1;

lr.noise = zeros(lr_ct,1);
lr.complete_rate = lr.noise;
lr.complete_time = lr.noise;
lr.path_eff = lr.noise;
lr.numTrials = lr.noise;
lr.complete_ind = cell(lr_ct,1);
lr.error = lr.complete_ind;
lr_i = 1;

for i = 1:length(testfolders)
    file = char(testfolders(i));
    respath = strcat(parent,'\results\');
    resfile = strcat(respath, file, '.mat');
    if exist(resfile, 'file') == 0
        disp(['getting ' file ' results...'])
        [test, daq, pvd] = loadTest(strcat(parent,'\',file));
        %% calculate metrics
        results = calcMetrics(test);
        
        %% save results
        if exist(respath,'dir') == 0
            mkdir(respath);
        end
        
        save(resfile, 'test', 'daq', 'pvd', 'results');
    else
        disp(['loading ' file ' results...'])
        load(resfile);
    end
    %%
    if test.ctrl == 1
        pr.noise(pr_i) = test.noise;
        pr.complete_rate(pr_i) = results.complete_rate;
        pr.complete_time(pr_i) = mean(results.complete_time);
        pr.path_eff(pr_i) = mean(results.path_eff(results.complete_ind));
        pr.numTrials(pr_i) = results.numTrials;
        pr.complete_ind{pr_i} = results.complete_ind;
        pr.error{pr_i} = results.error;
        pr_i = pr_i + 1;
    elseif test.ctrl == 2
        lr.noise(lr_i) = test.noise;
        lr.complete_rate(lr_i) = results.complete_rate;
        lr.complete_time(lr_i) = mean(results.complete_time);
        lr.path_eff(lr_i) = mean(results.path_eff(results.complete_ind));
        lr.numTrials(lr_i) = results.numTrials;
        lr.complete_ind{lr_i} = results.complete_ind;
        lr.error{lr_i} = results.error;
        lr_i = lr_i + 1;
    end
end

%% all plots
c = colormap(lines);
pr.noise = pr.noise/100;
lr.noise = lr.noise/100;

ind = strfind(parent,'\');
name = parent(ind(end) + 1:end);
[pr.noise, noise_pr] = sort(pr.noise);
[lr.noise, noise_lr] = sort(lr.noise);

% completion rate
subplot(311)
hold all
plot(pr.noise, pr.complete_rate(noise_pr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.complete_rate(noise_lr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Completion Rate')
axis tight
ylim([0 1])
title([name ' Results'])
legend('LDA','LR')

% completion time
subplot(312)
hold all
plot(pr.noise, pr.complete_time(noise_pr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.complete_time(noise_lr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Completion Time (s)')
axis tight
ylim([0 10])

% path efficiency
subplot(313)
hold all
plot(pr.noise, pr.path_eff(noise_pr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(1,:))
plot(lr.noise, lr.path_eff(noise_lr), '.-', 'MarkerSize', 25, 'LineWidth', 1.5, 'Color', c(2,:))
ylabel('Path Efficiency')
xlabel('Noise (V)')
axis tight
ylim([0 1])
%% euclidean distance for completed trials (LDA)
pr.c = linspecer(pr_ct + 2,'blue');
lr.c = linspecer(lr_ct + 2,'red');
pr.c = pr.c(3:end,:);
lr.c = lr.c(3:end,:);
t = 0:.1:12;
figure
for i = 1:pr_ct
    subplot(pr_ct,1,i)
    hold all
    for j = 1:length(pr.complete_ind{i})
        error = pr.error{i}{pr.complete_ind{i}(j)};
        if ~isempty(error)
            plot(t(1:length(error)),error, 'Color', pr.c(i,:), 'LineWidth', 1.5)
        end
    end
    xlim([0 10])
    ylim([0 3])
    ylabel([num2str(pr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Completed Trials (LDA)'])
    end
end
xlabel('Time (s)')
saveas(gcf, [respath name '_ldaComplete.png'])
%% euclidean distance for failed trials LDA
pr_error = zeros(pr_ct,100);
figure
for i = 1:pr_ct
    subplot(pr_ct,1,i)
    hold all
    fail_ind = setdiff(1:pr.numTrials(i),pr.complete_ind{i});
    error = zeros(length(fail_ind), 100);
    if ~isempty(fail_ind)
        for j = 1:length(fail_ind)
            error(j,:) = pr.error{i}{fail_ind(j)}(1:100);
            plot(t(1:length(error(j,:))), error(j,:), 'Color', pr.c(i,:), 'LineWidth',1)
        end
        pr_error(i,:) = mean(error);
        plot(t(1:length(pr_error(i,:))),mean(error), 'k', 'LineWidth',1.5)
    end
    ylabel([num2str(pr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Failed Trials (LDA)'])
    end
    xlim([0 10])
    ylim([0 6])
end
xlabel('Time (s)')
saveas(gcf, [respath name '_ldaFailed.png'])
%% euclidean distance for completed trials LR
figure
for i = 1:lr_ct
    subplot(lr_ct,1,i)
    hold all
    for j = 1:length(lr.complete_ind{i})
        error = lr.error{i}{lr.complete_ind{i}(j)};
        if ~isempty(error)
            plot(t(1:length(error)), error, 'Color', lr.c(i,:), 'LineWidth', 1.5)
        end
    end
    xlim([0 10])
    ylim([0 3])
    ylabel([num2str(lr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Completed Trials (LR)'])
    end
end
xlabel('Time (s)')
saveas(gcf, [respath name '_lrComplete.png'])

lr_error = zeros(lr_ct,100);
figure
for i = 1:lr_ct
    subplot(lr_ct,1,i)
    hold all
    fail_ind = setdiff(1:lr.numTrials(i),lr.complete_ind{i});
    error = zeros(length(fail_ind), 100);
    if ~isempty(fail_ind)
        for j = 1:length(fail_ind)
            error(j,:) = lr.error{i}{fail_ind(j)}(1:100);
            plot(t(1:length(error(j,:))), error(j,:), 'Color', lr.c(i,:), 'LineWidth',1)
        end
        lr_error(i,:) = mean(error);
        plot(t(1:length(lr_error(i,:))), mean(error), 'k', 'LineWidth',1.5)
        
    end
    ylabel([num2str(lr.noise(i)) ' V'])
    if i == 1
        title([name ' Euclidean Distance for Failed Trials (LR)'])
    end
    xlim([0 10])
    ylim([0 6])
end
xlabel('Time (s)')
saveas(gcf, [respath name '_lrFailed.png'])
%%
figure
subplot(121)
hold all
for i = 1:pr_ct
    if sum(pr_error(i,:) > 0)
        plot(t(1:length(pr_error(i,:))), pr_error(i,:),'Color', pr.c(i,:), 'LineWidth',1.5)
    else
        plot(0,0,'Color', pr.c(i,:))
    end
end
legend(num2str(pr.noise))
title([name ' Mean Euclidean Distance (LDA)'])
ylim([0 6])


% c_mean = colormap(winter(lr_ct));
subplot(122)
hold all
for i = 1:lr_ct
    if sum(lr_error(i,:) > 0)
        plot(t(1:length(lr_error(i,:))),lr_error(i,:),'Color', lr.c(i,:), 'LineWidth',1.5)
    else
        plot(0,0,'Color', lr.c(i,:))
    end
end
legend(num2str(lr.noise))
title([name ' Mean Euclidean Distance (LR)'])
xlabel('Time (s)')
ylim([0 6])
saveas(gcf, [respath name '_meanFailed.png'])