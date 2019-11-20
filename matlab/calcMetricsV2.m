%% grab variables
t = test.t;
inTrial = test.inTrial;
testTarget = test.testTarget;
inTarget = test.inTarget;
completeN = test.completeN;
xPos = test.xPos;
yPos = test.yPos;

%% calculate metrics
% completion rate
numTrials = sum(diff(inTrial)==1);
complete_total = completeN(end);
complete_rate = complete_total/numTrials;

% trial time stamps
end_trial = find(diff(inTrial) == -1);
start_trial = find(diff(inTrial) == 1);
complete_end = find(diff(completeN) == 1);
if complete_rate ~= 1
    complete_ind = zeros(1,length(complete_end));
    for i = 1:length(complete_end)
        complete_ind(i) = find(end_trial == complete_end(i));
    end
    complete_start = start_trial(complete_ind);
else
    complete_ind = 1:numTrials;
    complete_start = start_trial;
end

% completion time
complete_time = t(complete_end) - t(complete_start);

%% calculate path efficiency
dist = zeros(numTrials, 1);
path_eff = zeros(numTrials, 1);
error = cell(numTrials,1);
theta = [pi/2, -pi/2, pi, 0, 3*pi/4, -3*pi/4, pi/4, -pi/4];

for i = 1:numTrials
    x_trial = xPos(start_trial(i):end_trial(i));
    y_trial = yPos(start_trial(i):end_trial(i));
    target = testTarget(end_trial(i));
    x_target = ones(length(x_trial),1)*.25*7*cos(theta(target));
    y_target = ones(length(x_trial),1)*.25*7*sin(theta(target));
    error{i} = sqrt((x_trial-x_target).^2 + (y_trial-y_target).^2);
    dist(i) = sum(sqrt(diff(x_trial).^2 + diff(y_trial).^2));
    path_eff(i) = 1.45/dist(i);
end

%%
in_time = zeros(numTrials, 1);
dt = diff(test.t);

for i = 1:numTrials
    stat = inTarget(start_trial(i):end_trial(i));
    in_time(i) = sum(dt(stat == 1));
end


%% create results structure
results.complete_rate = complete_rate;
results.complete_time = complete_time;
results.path_eff = path_eff;
results.error = error;
results.complete_ind = complete_ind;
results.numTrials = numTrials;