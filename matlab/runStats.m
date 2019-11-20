clear
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% load data
[file, path] = uigetfile();
load([path file]);

%% pr vs. lr
h = zeros(5,8);
p = h;

for i = 1:5
    for j = 3:2:9
        [h(i,j),p(i,j)] = ttest(results{i:5:end,j},results{i:5:end,j+1});
    end
end

%% control vs. noise 
h1 = zeros(4,8);
p1 = h1;

for i = 1:4
    for j = 3:10
        [h1(i,j),p1(i,j)] = ttest(results{1:5:end,j},results{i+1:5:end,j});
    end
end

