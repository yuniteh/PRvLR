clear
close
clc
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontName','Cambria')

%%
noise = .25:.25:1;
c = colormap(lines);
figure(1)

parent = uigetdir();
folders = dir(parent);
flags = [folders.isdir];
folders = folders(flags);
ABfolders = {folders(strncmpi({folders.name}, 'AB', 2)).name};

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
        mode(ctrl,:) = pvd.MODE(ind);
        ave = pvd.CHAN_MAV(ind,1:2:end);
        uq(ctrl,:) = mean(ave) + std(ave);
        
    end
    ch = 1:3;
    a(i) = subplot(2,length(ABfolders),i);
    title(ABfolders{i})
    hold all
    aa = bar(ch',[uq(1,:)', uq(2,:)']);
    aa(1).FaceColor = c(1,:);
    aa(2).FaceColor = c(2,:);
    pause(.01)
    prx = ch' + aa(1).XOffset;
    lrx = ch' + aa(2).XOffset;
    
    b(i) = subplot(2,length(ABfolders),i + length(ABfolders));
    hold all
    for iii = 1:size(uq,2)
        plot(prx(iii)*ones(size(noise)), uq(1,iii)./noise, '.-', 'MarkerSize', 15, 'Color', c(1,:), 'LineWidth', 1.5)
        plot(lrx(iii)*ones(size(noise)), uq(2,iii)./noise, '.-', 'MarkerSize', 15, 'Color', c(2,:), 'LineWidth', 1.5)
    end
end
%%
ylabel(a(1),'MVC (V)');
ylabel(b(1),'SNR');
xlabel(b(1),'Electrode');
set([a b(2:end)],'XTickLabel',[])
set(b(1),'xtick',[1 2 3])
set([a(2:end) b(2:end)],'yticklabel',[])
a(1).YLim = [0, 1.3];
b(1).YLim = [0 6];
linkaxes(a,'y')
linkaxes(b,'y')
linkaxes([a,b],'x')

