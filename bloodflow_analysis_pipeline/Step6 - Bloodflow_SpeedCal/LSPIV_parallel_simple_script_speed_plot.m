% Create combined velocity plot
velocity_fig = figure('Name', 'All Velocity Profiles', 'Position', [100, 100, 800, 400]);
hold on;
box on;
set(gca, 'TickDir', 'out');

% Plot all velocity curves
y_min = inf;
y_max = -inf;
x_min = inf;
x_max = -inf;
mean_velocities = zeros(filenum,1);

plotNo = [1, 2, 3];
%plotNo = [1, 2, 6, 8, 9, 11, 14, 19, 21, 23, 24];
%plotNo = [2, 3, 4, 6, 7];
%plotNo = [1, 2, 3, 4, 5, 6, 7, 8, 9];
%plotNo = [1, 3, 5, 6, 7, 8, 9, 10];

legend_entries = {};  % 初始化
plot_handles = [];    % 用于 legend

for k = 1:length(plotNo)
    i = plotNo(k);  % 实际数据索引
    if ~isempty(all_velocities{i})
        h = plot(all_times{i}, all_velocities{i}, 'Color', colors(i,:), 'LineWidth', 1.5);
        plot_handles(end+1) = h;  % 保存句柄用于 legend
        mean_velocities(i) = mean(all_velocities{i}, 'omitnan');
        legend_entries{end+1} = sprintf('Vessel %d (%.2f mm/s)', i, mean_velocities(i));

        y_min = min(y_min, min(all_velocities{i}));
        y_max = max(y_max, max(all_velocities{i}));
        x_min = min(x_min, min(all_times{i}));
        x_max = max(x_max, max(all_times{i}));
    end
end

% Set axis limits and labels
xlim([x_min x_max]);
y_range = y_max - y_min;
ylim([y_min-0.1*y_range y_max+0.1*y_range]);
xlabel('Time (ms)');
ylabel('Velocity (mm/s)');
title('Blood Flow Velocity Profiles');

% Add legend
legend(plot_handles, legend_entries, 'Location', 'best');
