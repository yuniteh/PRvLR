clear
close all
clc

%%
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',12)

%% load data
parent = uigetdir();
files = dir(fullfile(parent,'*data.mat'));
%%
ax = [3 7 5 1 4 6 2 8];
for i = 1:length(files)
    load([parent '/' files(i).name])
    noise = unique([data.noise]);
    for ii = 1:length(data)
        for iii = 1:8
            targs = data(ii).results.targets(data(ii).results.complete_ind);
            path_eff = data(ii).results.path_eff(data(ii).results.complete_ind);
            if data(ii).ctrl == 1
                pr.t_rate(i,noise == data(ii).noise, ax == iii) = length(find(targs == iii))/3;
                pr.path_eff(i,noise == data(ii).noise, ax == iii) = nanmean(path_eff(targs == iii));
                pr.in_time(i,noise == data(ii).noise, ax == iii) = nanmean(data(ii).results.in_time(data(ii).results.targets == iii));
                pr.complete_time(i,noise == data(ii).noise, ax == iii) = nanmean(data(ii).results.complete_time(targs == iii));
            else
                lr.t_rate(i,noise == data(ii).noise, ax == iii) = length(find(targs == iii))/3;
                lr.path_eff(i,noise == data(ii).noise, ax == iii) = nanmean(path_eff(targs == iii));
                lr.in_time(i,noise == data(ii).noise, ax == iii) = nanmean(data(ii).results.in_time(data(ii).results.targets == iii));
                lr.complete_time(i,noise == data(ii).noise, ax == iii) = nanmean(data(ii).results.complete_time(targs == iii));
            end
        end
    end
end
lr.t_rate = lr.t_rate*100;
lr.path_eff = lr.path_eff * 100;
pr.t_rate = pr.t_rate*100;
pr.path_eff = pr.path_eff * 100;
if mean(noise) > 10
    noise = noise/100;
end

%%
pr_c = linspecer(length(noise) + 2,'blue');
lr_c = linspecer(length(noise) + 2,'red');
pr_c = pr_c(3:end,:);
lr_c = lr_c(3:end,:);
%%
for i = 1:length(noise)
    leg{i} = num2str(noise(i));
end
% for i = 1:length(pr.t_rate)
%     figure(1)
subplot(1,2,1)
hold all
meanrate = nanmean(pr.t_rate,1);
meanrate = reshape(meanrate,[5,8]);
meanrate(isnan(meanrate)) = 0;
[f, ca, o1,t] = spider(meanrate','LDA',ones(8,1),[],leg',1);

subplot(1,2,2)
hold all
meanrate = nanmean(lr.t_rate,1);
meanrate = reshape(meanrate,[5,8]);
meanrate(isnan(meanrate)) = 0;
[f, ca, o2,t1] = spider(meanrate','LR',ones(8,1),[],leg',1);
for ii = 1:length(noise)
    set(o1(ii),'color',pr_c(ii,:),'linewidth',1.5);
    set(o2(ii),'color',lr_c(ii,:),'linewidth',1.5);
end
% end

%%
leg2{1} = 'LDA';
leg2{2} = 'LR';
c = colormap(lines);
c(2,:) = [c(2,1)+.15 c(2,2)+.15 c(2,3)];
fields = fieldnames(pr);
labels = {'Completion Rate (%)', 'Path Efficiency (%)', 'Time in Target (s)', 'Completion Time (s)'};
r = [100; 100; 10; 10];

for i = 1:length(fields)
    prfield = pr.(fields{i});
    lrfield = lr.(fields{i});
    
    meanpr = reshape(nanmean(prfield,1),[5,8]);
    meanpr(isnan(meanpr)) = 0;
    meanlr = reshape(nanmean(lrfield,1),[5,8]);
    meanlr(isnan(meanlr)) = 0;
    
    meanall = [meanpr(1,:); meanlr(1,:)];
    figure(i)
    [f, ca, o1,t] = spider(meanall',labels{i},r(i)*ones(8,1),[],leg2',i);
    set(o1(1),'color',c(1,:),'linewidth',1.5);
    set(o1(2),'color',c(2,:),'linewidth',1.5);
    
end
