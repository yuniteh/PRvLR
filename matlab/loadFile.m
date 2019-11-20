function [vars, daq, pvd, results] = loadFile(varargin)

%% convert & load data
if isempty(varargin)
    datafolder = uigetdir();
else
    datafolder = varargin{1};
end

if ~exist(strcat(datafolder,'\DATA\MAT\data.mat'),'file')
    convertDAQtoMAT(datafolder);
end
load(strcat(datafolder,'\DATA\MAT\data.mat'));

ind = strfind(datafolder, '\');
dataname = datafolder(ind(end) + 1:end);
if strncmpi(dataname, 'train', 5)
    type = 0;
elseif strncmpi(dataname, 'test', 4)
    type = 1;
elseif strncmpi(dataname, 'free', 4)
    type = 2;
end

%% general variables
vars.nEMG = size(data.daq.DAQ_DATA,2);
win = 25;
vars.tDown = data.daq.t(1:win:end);

%% read csv file
delim = strfind(dataname, '_');
if dataname(delim(1) + 1) == '1'
    vars.ctrl = 1;
    vars.wg = csvread(strcat(datafolder,'\wg.csv'));
    vars.cg = csvread(strcat(datafolder,'\cg.csv'));
elseif dataname(delim(1) + 1) == '2'
    vars.ctrl = 2;
    vars.w = csvread(strcat(datafolder,'\w.csv'));
end
vars.mvc = csvread(strcat(datafolder,'\mvc.csv'));

%% read text file
% time noise testchannel intrial target intarget completeN x y
if type > 0
    testdata = dlmread(strcat(datafolder,'\',dataname,'.txt'));
    
    vars.datafolder = datafolder;
    vars.dataname = dataname;
    
    %% correct time
    vars.t = testdata(:,1) - testdata(1,1);
    skip = find(diff(vars.t) > 4000) + 1;
    for i = 1:length(skip)
        vars.t(skip(i):end) = vars.t(skip(i):end) - 4000;
    end
    skip = find(diff(vars.t) > 40) + 1;
    for i = 1:length(skip)
        vars.t(skip(i):end) = vars.t(skip(i):end) - 40;
    end
    
    %% get data
    vars.noise = testdata(1,2);
    vars.testCh = testdata(:,3);
    vars.inTrial = testdata(:,4);
    vars.testTarget = testdata(:,5);
    vars.inTarget = testdata(:,6);
    vars.completeN = testdata(:,7);
    vars.xPos = testdata(:,8);
    vars.yPos = testdata(:,9);
    
    % old format
    %     test.inTrial = testdata(:,3);
    %     test.testTarget = testdata(:,4);
    %     test.inTarget = testdata(:,5);
    %     test.completeN = testdata(:,6);
    %     test.xPos = testdata(:,7);
    %     test.yPos = testdata(:,8);
end

%% set outputs
daq = data.daq;
pvd = data.pvd;

if type == 1
    results = calcMetrics(vars);
else
    results = {};
end
end