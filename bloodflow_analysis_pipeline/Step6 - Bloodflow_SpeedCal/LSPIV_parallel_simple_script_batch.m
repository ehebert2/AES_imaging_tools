% Script to run LSPIV_parallel_simple.m
clear;close all;
clc;

% Base path and filename
base_pathname = 'D:\Code\Matlab\20250609 3DBloodVesselDataAnalysis\Step5 - Nonlinear_Transform\Kymograph_Transformed\';
base_filename_1 = 'mouse1_AES_dualPort_50pt4Hz_684-732um_3pt4mW_0pt53mW_00002_MIP_chunk_13of40_kymograph_from_file_segment_';
base_filename_2 = '_frames_50-250_width5_transformed'; 

% Define colors for each curve [R G B]
% colors = [
%     192 000 000;  % 
%     237 125 049;  % 
%     091 155 213;  % 
%     084 130 053;  % 
%     255 255 000;  % 
%     146 208 080   % 
% ]/255;

colors = [
    091 155 213;  %
    084 130 053;  % 
    192 000 000;  % 
    237 125 049;  % 
    255 255 000;  % 
    146 208 080;  %
    255 192 000;
    143 170 220;
    255 000 000
]/255;

% colors = [
%     091 155 213;  %
%     084 130 053;  % 
%     192 000 000;  % 
%     237 125 049;  % 
%     255 255 000;  % 
%     146 208 080;  %
%     255 192 000;
%     143 170 220;
%     255 000 000;
%     191 144 000
% ]/255;
 
% colors = [
%     091 155 213;  %
%     084 130 053;  % 
%     192 000 000;  % 
%     237 125 049;  % 
%     255 255 000;  % 
%     146 208 080;  %
%     255 192 000;
%     143 170 220;
%     255 000 000;
%     191 144 000;
%     091 155 213;  %
%     084 130 053;  % 
%     192 000 000;  % 
%     237 125 049;  % 
%     255 255 000;  % 
%     146 208 080;  %
%     255 192 000;
%     143 170 220;
%     255 000 000;
%     191 144 000;
%     091 155 213;  %
%     084 130 053;  % 
%     192 000 000;  % 
%     237 125 049;  % 
%     255 255 000;  % 
%     146 208 080;  %
%     255 192 000
% ]/255;

% colors = [
%     084 130 053 %
% ]/255;

% colors = [
%     237 125 049 %
% ]/255;

filenum = 3;

% Initialize cell arrays to store velocity data
all_velocities = cell(filenum,1);
all_times = cell(filenum,1);

% Process files with suffixes 1 to 6
for i = 1:filenum
    % Construct full filename
    filename = [base_filename_1 num2str(i) base_filename_2 '.tif'];
    filepath = fullfile(base_pathname, filename);
    
    % Check if file exists
    if ~exist(filepath, 'file')
        warning('File not found: %s', filepath);
        continue;
    end
    
    % Read the kymograph
    disp(['Processing file: ' filename]);
    kymograph = imread(filepath);
    
    % Set parameters
    params = struct();
    params.frmRate = 100.8;         % Frame rate in Hz
    params.pxlSize = 0.8702;       % Pixel size in um
    params.numavgs = 10;          % Number of lines to average
    params.shiftamt = 1;          % Shift amount for correlation
    params.skipamt = 2;           % Skip amount for processing
    params.numWorkers = 12;       % Number of parallel workers
    params.kymo_path = filepath;  % Save the kymograph path for saving results
    
    % Run LSPIV analysis
    [velocity_mm, time, shifted_kymograph] = LSPIV_parallel_simple(kymograph, params);
    
    % Store velocity data
    all_velocities{i} = velocity_mm;
    all_times{i} = time;
    
    % Clear some variables to free memory
    clear kymograph shifted_kymograph;
end

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

for i = 1:filenum %plotNo%
    if ~isempty(all_velocities{i})
        plot(all_times{i}, all_velocities{i}, 'Color', colors(i,:), 'LineWidth', 1.5);
        y_min = min(y_min, min(all_velocities{i}));
        y_max = max(y_max, max(all_velocities{i}));
        x_min = min(x_min, min(all_times{i}));
        x_max = max(x_max, max(all_times{i}));
        mean_velocities(i) = mean(all_velocities{i}, 'omitnan');
    end
end

% Set axis limits and labels
xlim([x_min x_max]);
y_range = y_max - y_min;
ylim([y_min-0.1*y_range y_max+0.1*y_range]);
xlabel('Time (ms)');
ylabel('Velocity (mm/s)');
title('Blood Flow Velocity Profiles');

% Add legend with mean velocities
legend_entries = cell(filenum,1);
for i = 1:filenum %plotNo%
    if ~isempty(all_velocities{i})
        legend_entries{i} = sprintf('Vessel %d (%.2f mm/s)', i, mean_velocities(i));
    else
        legend_entries{i} = sprintf('Vessel %d (N/A)', i);
    end
end
legend(legend_entries, 'Location', 'best');

% Add mean velocity lines and text
% for i = 1:6
%     if ~isempty(all_velocities{i})
%         % Draw mean velocity line
%         line([x_min x_max], [mean_velocities(i) mean_velocities(i)], ...
%             'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1);
%     end
% end

% Save the combined velocity plot
[~, base_name] = fileparts(base_filename_1);
saveas(velocity_fig, fullfile(base_pathname, [base_name 'combined_velocity.png']));

% Save velocity data including mean velocities
save(fullfile(base_pathname, [base_name 'velocity_data.mat']), ...
    'all_velocities', 'all_times', 'colors', 'legend_entries', 'mean_velocities');

disp('Processing complete!');