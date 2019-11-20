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
ind = (data.pvd.MODE ~= -1) & (data.pvd.COLLECT == 1); % get data collection indices
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
plot(tDown,mode)
for i = 1:nEMG
    a(i) = subplot(nEMG+1,1,i+1);
    plot(data.daq.t,data.daq.DAQ_DATA(:,i))
end
xlabel('Time (s)')
linkaxes(a,'x')

%% get training features
feat = data.pvd.FEAT_RAW(mode ~= -1,:);
class = mode(mode ~= -1);

%% build lda classifier
lda = fitcdiscr(feat,class);
matOut = predict(lda, feat);
disp('matlab lda: ') 
disp(sum(matOut == class)/size(class,1))

save(strcat(datafolder,'\lda.mat'), 'lda');

%% get weights from pce
wg = csvread(strcat(datafolder,'\wg.csv'));
cg = csvread(strcat(datafolder,'\cg.csv'));

pceP = [feat ones(size(feat,1),1)] * [wg; cg];
[temp, pceOut] = max(pceP,[],2);
pceOut = pceOut - 1;
disp('python lda: ')
disp(sum(pceOut == class)/size(class,1))
%% get helper function weights
ldaW = LDA(feat, class)';
ldaP = [ones(size(feat,1),1) feat] * ldaW;
[temp, ldaOut] = max(ldaP,[],2);
ldaOut = ldaOut - 1;
%% compare classifiers
figure

subplot(311)
hold all
plot(class)
plot(pceOut)
subplot(312)
hold all
plot(class)
plot(matOut)
subplot(313)
hold all
plot(class)
plot(ldaOut)
