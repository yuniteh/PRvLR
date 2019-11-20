function plotError(results)
%% get variables
complete_ind = results.complete_ind;
numTrials = results.numTrials;
fail_ind = setdiff(1:numTrials, complete_ind);
error = results.error;

%% plot error over time
% complete_lines = linspecer(length(complete_ind), 'green');
% fail_lines = linspecer(length(fail_ind), 'red');
complete_lines = colormap(winter(length(complete_ind)));
fail_lines = colormap(autumn(length(fail_ind)));
t_plot = 0:.1:10;

figure
if results.complete_rate ~= 1
    % plot completed trials
    b(2) = subplot(2,1,1);
    title('Completed Trials')
    hold all
    for i = 1:length(complete_ind)
        y_plot = error{complete_ind(i)};
        plot(t_plot(1:length(y_plot)),y_plot,'Color', complete_lines(i,:), 'LineWidth', 1.5)
    end
    set(gca,'xticklabel',[])
    
    % plot failed trials
    b(1) = subplot(2,1,2);
    title('Failed Trials')
    hold all
    for i = 1:numTrials - length(complete_ind)
        plot(t_plot,error{fail_ind(i)},'Color', fail_lines(i,:), 'LineWidth', 1.5)
    end
    axis tight
    linkaxes(b,'x')
    xlabel('Time (s)')
else
    hold all
    for i = 1:numTrials
        plot(t_plot(1:length(error{i})),error{i}, 'Color', complete_lines(i,:), 'LineWidth', 1.5)
    end
    title('Completed Trials')
    xlabel('Time (s)')
    ylabel('Euclidean distance')
end
end