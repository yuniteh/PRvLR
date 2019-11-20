%
clear
%close
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',12)

%
[file, path] = uigetfile();
load([path file])

% INITIALIZE VARIABLES
labels = {'Completion Rate (%)', 'Path Efficiency (%)', 'Dwelling Time (s)', 'Completion Time (s)'};
c = colormap(lines);
c(2,:) = [c(2,1)+.15 c(2,2)+.15 c(2,3)];
noise = unique(results.noise)/100;

% LOOP THROUGH METRICS

figure(1)
ax = tight_subplot(4,1,.04,.09,.13);

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
    
%     CALC MEANS
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
    
%     SUBPLOT FOR LDA VS LR BARS ONLY
    axes(ax(i));
    hold all
    b = bar(noise,[prmean(:,1),lrmean(:,1)],1);
    b(1).FaceColor = c(1,:);
    b(2).FaceColor = c(2,:);
    pause(.01);
    
    prx = noise + b(1).XOffset;
    lrx = noise + b(2).XOffset;
    
    errorbar(prx,prmean(:,1),prsem,'.','Color','k');
    errorbar(lrx,lrmean(:,1),lrsem,'.','Color','k');
    
    for ii = 2:length(noise)
        if ii == 2 && i == 1
            p = .028;
        else
            p = .001;
        end
        if p < .05
            x = [noise(1) noise(ii)];
            sigstar(x,p);
        end
    end
    
    ylabel(labels{i})
    if i == 1
        leg = legend('LDA','LR');
        set(leg,'Location','northeastoutside');
        legend('boxoff')
    end
    
    axes(ax2(i));
    hold all
    b1 = bar(1,prmean(1,1),1);
    b2 = bar(2,lrmean(1,1),1);
    b1.FaceColor = c(1,:);
    b2.FaceColor = c(2,:);

    title(labels{i})
    pause(.01);
    
    prx = 1 + b1.XOffset;
    lrx = 2 + b2.XOffset;
    
    errorbar(prx,prmean(1,1),prsem(1),'.','Color','k');
    errorbar(lrx,lrmean(1,1),lrsem(1),'.','Color','k');
end

% SET AXES PROPERTIES
set(get(ax(end),'XLabel'),'String','Noise (V)');
set(ax,'XTickLabel',[]);
set(ax(end),'XTickLabel',noise);
set(ax,'XTick',noise);
linkaxes(ax,'x')
axis tight

% BASELINE PLOTS
set(ax2,'xtick',[])
linkaxes(ax2,'x')
set(ax2(1),'xlim',[0 3])