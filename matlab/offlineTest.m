clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
[vars, daq, pvd, results] = loadFile();

%% plot raw data
figure
hold all

for i = 1:vars.nEMG
    a(i) = subplot(vars.nEMG,1,i);
    plot(daq.t,daq.DAQ_DATA(:,i))
    if i == 1
        title('Raw EMG')
    end
end
xlabel('Time (s)')
linkaxes(a, 'x')

%% plot targets and trajectories
plotTargets();
plot(vars.xPos, vars.yPos)
title('Cursor Trajectory')

%% load matlab lda classifier
[ldafile, ldapath] = uigetfile('*lda.mat');
if ldafile ~= 0
    load(strcat(ldapath,ldafile));
    offPred = predict(lda,pvd.FEAT_RAW);
    
    % compare
    bool = offPred == pvd.CLASS_EST;
    matpy_rat = sum(bool)/size(bool,1);
end
%% calculate metrics
results = calcMetrics(vars);
plotError(results);
%% save metrics
%save(strcat(vars.datafolder,'\results\',vars.dataname,'.mat'), 'vars', 'daq', 'pvd', 'results');

%% classify noisy data
% noise_a = input('add noise (1/0)? ');
%
% if noise_a == 1
%     classOut = cell(10,1);
%     for i = 1:10
%         % extract features
%         feat_noisy = addNoise(data.daq.DAQ_DATA, max(mvc(4:5)), i/10);
%
%         % classify
%         pceP = [feat_noisy ones(size(feat_noisy,1),1)] * [wg; cg];
%         [temp, pceOut] = max(pceP,[],2);
%         pceOut = pceOut - 1;
%         classOut{i} = pceOut;
%     end
% end


