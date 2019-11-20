function [test, daq, pvd, results] = loadFile(datafolder)
% type 0 = train; type 1 = test; type 2 = free
% varargin = folders

%% convert & load data
if isempty(datafolder)
    datafolder = uigetdir();
end

if exist(strcat(datafolder,'\DATA\MAT\data.mat'),'file') == 0
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

%% read text file
% time noise testchannel intrial target intarget completeN x y
if type > 0
    if ~isempty(strfind(datafolder,'test'))
        testdata = dlmread(strcat(datafolder,'\',dataname,'.txt'));
        
        test.datafolder = datafolder;
        test.dataname = dataname;
        
        %% correct time
        test.t = testdata(:,1) - testdata(1,1);
        skip = find(diff(test.t) > 4000) + 1;
        for i = 1:length(skip)
            test.t(skip(i):end) = test.t(skip(i):end) - 4000;
        end
        skip = find(diff(test.t) > 40) + 1;
        for i = 1:length(skip)
            test.t(skip(i):end) = test.t(skip(i):end) - 40;
        end
        
        %% get data
        test.noise = testdata(1,2);
        test.testCh = testdata(:,3);
        test.inTrial = testdata(:,4);
        test.testTarget = testdata(:,5);
        test.inTarget = testdata(:,6);
        test.completeN = testdata(:,7);
        test.xPos = testdata(:,8);
        test.yPos = testdata(:,9);
        
        % old format
        %     test.inTrial = testdata(:,3);
        %     test.testTarget = testdata(:,4);
        %     test.inTarget = testdata(:,5);
        %     test.completeN = testdata(:,6);
        %     test.xPos = testdata(:,7);
        %     test.yPos = testdata(:,8);
    end
end

%% read csv file
if dataname(6) == '1'
    test.ctrl = 1;
    test.wg = csvread(strcat(datafolder,'\wg.csv'));
    test.cg = csvread(strcat(datafolder,'\cg.csv'));
elseif dataname(6) == '2'
    test.ctrl = 2;
    test.w = csvread(strcat(datafolder,'\w.csv'));
end
test.mvc = csvread(strcat(datafolder,'\mvc.csv'));

%% general variables
test.nEMG = size(data.daq.DAQ_DATA,2);
win = 25;
test.tDown = data.daq.t(1:win:end);

%% set outputs
daq = data.daq;
pvd = data.pvd;

if type == 1
    results = calcMetrics(test);
else
    results = {};
end
end