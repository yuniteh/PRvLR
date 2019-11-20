%%
clear
close
clc

%%
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Calibri')
set(0,'DefaultAxesFontSize',12)

%%
[file, path] = uigetfile();
load([path file])

%% INITIALIZE VARIABLES
labels = {'Completion Rate (%)', 'Path Efficiency (%)', 'Dwelling Time (s)', 'Completion Time (s)'};
lims = [100;100;10;10];
c = colormap(lines);
c(1,:) = [0 .298 0.6];
c(2,:) = [1 0.475 0.098];
% c(2,:) = [c(2,1)+.15 c(2,2)+.15 c(2,3)];
noise = unique(results.noise)/100;

%% LOOP THROUGH METRICS

figure(1)
ax = tight_subplot(2,2,.09,.15,.13);

figure(2)
ax2 = tight_subplot(1,4,.05,.15,.13);

for i = 1:size(results,2) - 3
    prfield = table2array(results(results.ctrl == 1, i + 3));
    lrfield = table2array(results(results.ctrl == 2, i + 3));
    prmean = zeros(length(noise), 2);
    lrmean = prmean;
    prsem = zeros(length(noise), 1);
    lrsem = prsem;
    
    if i < 3
        scale = 100;
    else
        scale = 1;
    end
    
    % CALC MEANS
    for ii = 1:length(noise)
        prmean(ii,1) = nanmean(prfield(ii:length(noise):end));
        lrmean(ii,1) = nanmean(lrfield(ii:length(noise):end));
        prmean(ii,2) = nanstd(prfield(ii:length(noise):end));
        lrmean(ii,2) = nanstd(lrfield(ii:length(noise):end));
        prsem(ii) = prmean(ii,2)/sqrt(max(results.sub));
        lrsem(ii) = lrmean(ii,2)/sqrt(max(results.sub));
    end
    
    prmean = prmean*scale;
    lrmean = lrmean*scale;
    prsem = prsem*scale;
    lrsem = lrsem*scale;
    
    % SUBPLOT FOR LDA VS LR BARS ONLY
    if i < 3
        ind = 1;
    else
        ind = 2;
    end
    axes(ax(i));
    hold all
    stdshade(prmean(:,1),prsem,.3,c(1,:),noise);
    stdshade(lrmean(:,1),lrsem,.3,c(2,:),noise);
    ylabel(labels{i})
    if i > 2
        xlabel('Noise (V)')
        %set(ax,'XTickLabel',[]);
        set(ax(i),'XTickLabel',noise);
        set(ax,'XTick',noise);
        ylim([0 10])
    else
        ylim([0 100])
    end
    
end

%% SET AXES PROPERTIES
set(ax,'XTickLabel',[]);
set(ax(end-1:end),'XTickLabel',noise);
set(ax,'XTick',noise);
linkaxes(ax,'x')
axis tight

%% BASELINE PLOTS
set(ax2,'xtick',[])
linkaxes(ax2,'x')
set(ax2(1),'xlim',[0 3])