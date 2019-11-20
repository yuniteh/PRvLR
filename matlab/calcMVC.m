clear
close
clc
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',12)

%%
ground_noise = .1:.1:.4;
electrode_noise = .25:.25:1;
colors = colormap(lines);
figure(1)

parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
ABfolders = {folders(strncmpi({folders.name}, 'AB', 2)).name};

ctrl_all = zeros(length(ABfolders)*2,1);
uq_all = ctrl_all;

for i = 1:length(ABfolders)
    path = [parent '/' ABfolders{i} '/mvc/results/'];
    files = dir(fullfile(path,'*.mat'));
    
    for ii = 1:length(files)
        load([path files(ii).name]);
        nEMG = size(daq.DAQ_DATA,2);
        win = 25;
        tDown = daq.t(1:win:end);
        tDown = tDown(1:size(pvd.CTRL,1));
        disp(files(ii).name)
        
        if ii == 1
            ctrl_bool = input('LDA (0) or LR(1)? ');
        else
            ctrl_bool = ~ctrl_bool;
        end
        ctrl = ctrl_bool + 1;
        
        ind = (pvd.MODE > 0) & (pvd.COLLECT == 1); % get data collection indices
        mode = pvd.MODE(ind);
        ave = pvd.CUR_VAL(ind);
        chan_ave = pvd.CHAN_MAV(ind,1:2:end);
        uq = mean(ave) + std(ave);
        
        ctrl_all((i - 1)*2 + ii) = ctrl;
        uq_all((i - 1)*2 + ii) = uq;
        
        figure(1)
        a(i) = subplot(2,length(ABfolders),i);
        title(ABfolders{i})
        hold all
        bar(ctrl, uq, 1, 'FaceColor', colors(ctrl,:))
        
        b(i) = subplot(2,length(ABfolders),i + length(ABfolders));
        hold all
        plot(ctrl*ones(size(ground_noise)), uq./ground_noise, '.-', 'MarkerSize', 15, 'Color', colors(ctrl,:), 'LineWidth', 1.5)
        
    end
end
ylabel(a(1),'MVC (V)');
ylabel(b(1),'SNR');
set([a(2:end) b(2:end)],'yticklabel',[])
set([a b],'xtick',[])
b(1).YLim = [0,7];
a(1).YLim = [0,.6];
a(1).XLim = [0,3];
linkaxes(a,'y')
linkaxes(b,'y')
linkaxes([a,b],'x')

%%
figure
ax = tight_subplot(2,1,.04,.1,.25);
uq_pr = uq_all(ctrl_all == 1);
uq_lr = uq_all(ctrl_all == 2);
mean_pr = mean(uq_pr);
mean_lr = mean(uq_lr);
ctrl = [1,2];
axes(ax(1));
hold on
d(1) = bar(1,mean_pr,1);
d(2) = bar(2, mean_lr,1);
d(1).FaceColor = colors(1,:);
d(2).FaceColor = colors(2,:);

pause(.01);

prx = 1 + d(1).XOffset;
lrx = 2 + d(2).XOffset;

errorbar(prx,mean_pr,nanstd(uq_pr)/sqrt(length(uq_pr)),'.','Color','k');
errorbar(lrx,mean_lr,nanstd(uq_lr)/sqrt(length(uq_lr)),'.','Color','k');

p = ranksum(uq_pr,uq_lr);

axes(ax(2))
hold all
plot(ones(size(ground_noise)), mean_pr./ground_noise, '.-', 'MarkerSize', 15, 'Color', colors(1,:), 'LineWidth', 1.5)
plot(2*ones(size(ground_noise)), mean_lr./ground_noise, '.-', 'MarkerSize', 15, 'Color', colors(2,:), 'LineWidth', 1.5)

linkaxes(ax, 'x');

set(ax(1),'xlim',[0 3])
set(ax,'xtick',[1 2])
set(ax(1),'xtick',[])
set(ax(2),'xticklabel',[{'LDA'},{'LR'}])
ylabel(ax(1),'MVC (V)');
ylabel(ax(2),'SNR');

%%
figure
ax = tight_subplot(2,1,.04,.1,.25);
uq_pr = uq_all(ctrl_all == 1);
uq_lr = uq_all(ctrl_all == 2);
mean_pr = mean(uq_pr);
mean_lr = mean(uq_lr);

sem_pr = zeros(size(ground_noise));
sem_lr = sem_pr;

for i = 1:length(ground_noise)
    sem_pr(i) = std(uq_pr./ground_noise(i))/sqrt(6);
    sem_lr(i) = std(uq_lr./ground_noise(i))/sqrt(6);
end

ctrl = [1,2];
axes(ax(1));
hold on
d(1) = bar(1,mean_pr,1);
d(2) = bar(2, mean_lr,1);
d(1).FaceColor = colors(1,:);
d(2).FaceColor = colors(2,:);

pause(.01);

prx = 1 + d(1).XOffset;
lrx = 2 + d(2).XOffset;

errorbar(prx,mean_pr,nanstd(uq_pr)/sqrt(length(uq_pr)),'.','Color','k');
errorbar(lrx,mean_lr,nanstd(uq_lr)/sqrt(length(uq_lr)),'.','Color','k');

p = ranksum(uq_pr,uq_lr);

axes(ax(2))
hold all
plot(ones(size(ground_noise)), mean_pr./ground_noise, '.', 'MarkerSize', 10, 'Color', colors(1,:), 'LineWidth', 1.5)
for i = 1:length(ground_noise)
    plot([1 1], [mean_pr./ground_noise(i) - sem_pr(i) mean_pr./ground_noise(i) + sem_pr(i)],'-','Color',colors(1,:),'Linewidth',1.5)
    plot([2 2], [mean_lr./ground_noise(i) - sem_lr(i) mean_lr./ground_noise(i) + sem_lr(i)],'-','Color',colors(2,:),'Linewidth',1.5)
    plot([1 2], [mean_pr./ground_noise(i) - sem_pr(i) mean_lr./ground_noise(i) - sem_lr(i)],'k--','Linewidth',1)
    plot([1 2], [mean_pr./ground_noise(i) + sem_pr(i) mean_lr./ground_noise(i) + sem_lr(i)],'k--','Linewidth',1)
end
plot(2*ones(size(ground_noise)), mean_lr./ground_noise, '.', 'MarkerSize', 10, 'Color', colors(2,:), 'LineWidth', 1.5)


linkaxes(ax, 'x');

set(ax(1),'xlim',[0 3])
set(ax,'xtick',[1 2])
set(ax(1),'xtick',[])
set(ax(2),'xticklabel',[{'LDA'},{'LR'}])
ylabel(ax(1),'MVC (V)');
ylabel(ax(2),'SNR');

%%
figure
ax2 = tight_subplot(2,1,.04,.1,.25);
uq_pr = uq_all(ctrl_all == 1);
uq_lr = uq_all(ctrl_all == 2);
mean_pr = mean(uq_pr);
mean_lr = mean(uq_lr);

sem_pr = zeros(size(electrode_noise));
sem_lr = sem_pr;

for i = 1:length(electrode_noise)
    sem_pr(i) = std(uq_pr./electrode_noise(i))/sqrt(6);
    sem_lr(i) = std(uq_lr./electrode_noise(i))/sqrt(6);
end

ctrl = [1,2];
axes(ax2(1));
hold on
d(1) = bar(1,mean_pr,1);
d(2) = bar(2, mean_lr,1);
d(1).FaceColor = colors(1,:);
d(2).FaceColor = colors(2,:);

pause(.01);

prx = 1 + d(1).XOffset;
lrx = 2 + d(2).XOffset;

errorbar(prx,mean_pr,nanstd(uq_pr)/sqrt(length(uq_pr)),'.','Color','k');
errorbar(lrx,mean_lr,nanstd(uq_lr)/sqrt(length(uq_lr)),'.','Color','k');

p = ranksum(uq_pr,uq_lr);

axes(ax2(2))
hold all
plot(ones(size(electrode_noise)), mean_pr./electrode_noise, '.', 'MarkerSize', 10, 'Color', colors(1,:), 'LineWidth', 1.5)
for i = 1:length(ground_noise)
    plot([1 1], [mean_pr./electrode_noise(i) - sem_pr(i) mean_pr./electrode_noise(i) + sem_pr(i)],'-','Color',colors(1,:),'Linewidth',1.5)
    plot([2 2], [mean_lr./electrode_noise(i) - sem_lr(i) mean_lr./electrode_noise(i) + sem_lr(i)],'-','Color',colors(2,:),'Linewidth',1.5)
    plot([1 2], [mean_pr./electrode_noise(i) - sem_pr(i) mean_lr./electrode_noise(i) - sem_lr(i)],'k--','Linewidth',1)
    plot([1 2], [mean_pr./electrode_noise(i) + sem_pr(i) mean_lr./electrode_noise(i) + sem_lr(i)],'k--','Linewidth',1)
end
plot(2*ones(size(electrode_noise)), mean_lr./electrode_noise, '.', 'MarkerSize', 10, 'Color', colors(2,:), 'LineWidth', 1.5)


linkaxes(ax2, 'x');

set(ax2(1),'xlim',[0 3])
set(ax2,'xtick',[1 2])
set(ax2(1),'xtick',[])
set(ax2(2),'xticklabel',[{'LDA'},{'LR'}])
ylabel(ax2(1),'MVC (V)');
ylabel(ax2(2),'SNR');
