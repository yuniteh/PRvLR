clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
[test, daq, pvd] = loadTest();

%% classify noisy data
classOut = cell(5,1);
for i = 1:5
    % extract features
    feat_noisy = addNoise(daq.daqUINT16, i/5);
    
    % classify
    pceP = [feat_noisy ones(size(feat_noisy,1),1)] * [test.wg; test.cg];
    [temp, pceOut] = max(pceP,[],2);
    pceOut = pceOut - 1;
    classOut{i} = pceOut;
end

%% get trial targets



