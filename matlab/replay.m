clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Cambria')

%% load data
ax = tight_subplot(1,2);

%%
[vars, daq, pvd, results] = loadFile();

%% plot full trajectory
ctrl_name = {'LDA','LR'};
co = colormap(lines);
figure(1)
axes(ax(vars.ctrl))
hold all
circles(0,0,3.5,'facecolor','none')
theta = [pi/2, -pi/2, pi, 0, 3*pi/4, -3*pi/4, pi/4, -pi/4];
x = .25*7*cos(theta);
y = .25*7*sin(theta);
grid off
c = circles(x,y,.4,'facecolor','none','edgecolor','green','linewidth',1.5);
%p = circles(0,0,.1, 'edgecolor',co(vars.ctrl,:),'facecolor', co(vars.ctrl,:));

end_trial = find(diff(vars.completeN) == 1);
start_trial = find(diff(vars.inTrial) == 1);
for i = 1:8
    start = start_trial(find(results.targets == i,1));
    end_i = end_trial(find(results.targets == i,1));
    plot(vars.xPos(start:end_i),vars.yPos(start:end_i),'LineWidth',1.5,'Color', co(vars.ctrl,:));
end
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XColor','none')
set(gca,'YColor','none')
axis square
title(ctrl_name(vars.ctrl))
%% plot targets and trajectories
figure
hold all
circles(0,0,3.5,'facecolor','none')
theta = [pi/2, -pi/2, pi, 0, 3*pi/4, -3*pi/4, pi/4, -pi/4];
x = .25*7*cos(theta);
y = .25*7*sin(theta);
grid off
for i = 1:length(vars.xPos)
    if vars.inTrial(i) > 0
        c = circles(x(vars.testTarget(i)),y(vars.testTarget(i)),.4,'facecolor','none','edgecolor','green');
    else
        c = [];
    end
    p = circles(vars.xPos(i), vars.yPos(i), .1, 'facecolor', 'black');
    title(['t = ' num2str(floor(vars.t(i)))])
    drawnow
    pause(.01)
    delete(p)
    delete(c)
end