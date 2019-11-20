function results = calcMetrics(test)
%% parameters
trial_t = 10; % length of trial
fs = 10; % sampling rate

%% grab variables
t = test.t;
inTrial = test.inTrial;
testTarget = test.testTarget;
testCh = test.testCh;
inTarget = test.inTarget;
completeN = test.completeN;
xPos = test.xPos;
yPos = test.yPos;
testTarget = test.testTarget;

%% calculate metrics
% trial time stamps
end_trial = find(diff(inTrial) == -1);
start_trial = find(diff(inTrial) == 1);
start_trial = start_trial(1:length(end_trial));
complete_end = find(diff(completeN) == 1);
targs = testTarget(end_trial);

%%
% completion rate
numTrials = length(start_trial);
complete_total = completeN(end);
complete_rate = complete_total/numTrials;

if complete_rate ~= 1
    complete_ind = zeros(1,length(complete_end));
    for i = 1:length(complete_end)
        complete_ind(i) = find(complete_end(i)<=end_trial & complete_end(i)>start_trial);
    end
    complete_start = start_trial(complete_ind);
else
    complete_ind = 1:numTrials;
    complete_start = start_trial;
end

fail_ind = setdiff(1:numTrials, complete_ind);

% target completion rate
chan = testCh(complete_end);

% completion time
complete_time = zeros(numTrials,1);
complete_time(complete_ind) = t(complete_end) - t(complete_start);
complete_time(fail_ind) = t(end_trial(fail_ind)) - t(start_trial(fail_ind));

% path efficiency & error
dist = zeros(numTrials, 1);
temp = cell(numTrials,1);
error = zeros(numTrials, trial_t*fs);
path_eff = zeros(numTrials, 1);
theta = [pi/2, -pi/2, pi, 0, 3*pi/4, -3*pi/4, pi/4, -pi/4];

for i = 1:numTrials
    x_trial = xPos(start_trial(i):end_trial(i));
    y_trial = yPos(start_trial(i):end_trial(i));
    target = testTarget(end_trial(i));
    x_target = ones(length(x_trial),1)*.25*7*cos(theta(target));
    y_target = ones(length(x_trial),1)*.25*7*sin(theta(target));
    temp{i} = sqrt((x_trial-x_target).^2 + (y_trial-y_target).^2);
    error(i,:) = temp{i}(1:trial_t*fs)';
    
    if any(complete_ind==i)
        x_trial = xPos(start_trial(i):complete_end(complete_ind == i));
        y_trial = yPos(start_trial(i):complete_end(complete_ind == i));
    end
    dist(i) = sum(sqrt(diff(x_trial).^2 + diff(y_trial).^2));
    if dist(i) < 1.5
        path_eff(i) = 1;
    else
        path_eff(i) = 1.5/dist(i);
    end
    
    if any(fail_ind==i)
        path_eff(i) = 0;
    end
        
end

% PATH EFFICIENCY FOR COMPLETED TRIALS ONLY
% path_eff = zeros(length(complete_end), 1);
% for i = 1:length(complete_end)
%     x_trial = xPos(complete_start(i):complete_end(i));
%     y_trial = yPos(complete_start(i):complete_end(i));
%     dist(i) = sum(sqrt(diff(x_trial).^2 + diff(y_trial).^2));
%     if dist(i) < 1.5
%         path_eff(i) = 1;
%     else
%         path_eff(i) = 1.5/dist(i);
%     end
% end

% time in target
in_time = zeros(numTrials, 1);
for i = 1:numTrials
    in_time(i) = sum(inTarget(start_trial(i):end_trial(i)) == 1) * .1;
end

results.in_time = in_time;
results.complete_rate = complete_rate;
results.complete_time = complete_time;
results.path_eff = path_eff;
results.error = error;
results.complete_ind = complete_ind;
results.numTrials = numTrials;
results.targets = targs;
results.channel = chan;

end