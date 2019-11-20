clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% convert & load data
datafolder = uigetdir();
if exist(strcat(datafolder,'/DATA/MAT/data.mat'),'file') == 0
    convertDAQtoMAT(datafolder);
end
load(strcat(datafolder,'/DATA/MAT/data.mat'));

%% general variables
nEMG = size(data.daq.DAQ_DATA,2);
win = 25;
tDown = data.daq.t(1:win:end);
tDown = tDown(1:size(data.pvd.CTRL,1));

%% plot raw data
figure
hold all
for i = 1:nEMG
    a(i) = subplot(nEMG,1,i);
    plot(data.daq.t,data.daq.DAQ_DATA(:,i))
end
xlabel('Time (s)')
linkaxes(a, 'x')

%% align flags
ind = (data.pvd.MODE ~= -1) & (data.pvd.Y_LABEL ~= -1) & (data.pvd.COLLECT == 1); % get data collection indices
ylabel = -ones(size(data.pvd.MODE));
mode = ylabel;
collect = ylabel;
ylabel(ind) = data.pvd.Y_LABEL(ind);
mode(ind) = data.pvd.MODE(ind);
collect(ind) = data.pvd.COLLECT(ind);
%% plot flags over data
figure
hold all
subplot(nEMG+1,1,1)
hold all
plot(tDown,ylabel)
plot(tDown,mode)
for i = 1:nEMG
    a(i) = subplot(nEMG+1,1,i+1);
    plot(data.daq.t,data.daq.DAQ_DATA(:,i))
end
xlabel('Time (s)')
linkaxes(a,'x')

%% get training features
nModes = 4;
feat = cell(nModes,1);
for i = 1:nModes
    y = ylabel(mode == i);
    temp = data.pvd.FEAT_RAW(mode == i,:);
    zeroMat = data.pvd.FEAT_RAW(mode ~= i & mode > 0,:);
    zeroMat = [zeroMat zeros(size(zeroMat,1),1)];
    feat{i} = [zeroMat; temp y];
end
%% pca & LR
for i = 1:nModes
    coeff{i} = pca(feat{i}(:,1:end-1));
    test = feat{i}(:,1:end-1)*coeff{i};
    plot3(test(:,1),test(:,2),test(:,3),'o')
    model{i} = fitlm(feat{i}(:,1:end-1),feat{i}(:,end));
    figure
    plot(model{i});
end
