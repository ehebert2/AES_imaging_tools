% Script to run LSPIV_parallel_simple.m
clear;close all;
clc;

% Base path and filename
base_pathname = 'E:\Volumetric blood vessel imaging data\kymo_mouseA2_video4_1-3200vols_Copy\Kymograph_Transformed\';
base_filename_1 = 'MAX_mouseA2_subY_aM0pt8_50pt4Hz_890-940um_flip_video4_vol1-9600_kymo_segment_';
base_filename_2 = '_frames_0-3200_width5_transformed'; 

%filenum = [1 2 3 4 5 6 7 8 9 10];
filenum = [27];

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
 
colors = [
    091 155 213;  %
    084 130 053;  % 
    192 000 000;  % 
    237 125 049;  % 
    255 255 000;  % 
    146 208 080;  %
    255 192 000;
    143 170 220;
    255 000 000;
    191 144 000;
    091 155 213;  %
    084 130 053;  % 
    192 000 000;  % 
    237 125 049;  % 
    255 255 000;  % 
    146 208 080;  %
    255 192 000;
    143 170 220;
    255 000 000;
    191 144 000;
    091 155 213;  %
    084 130 053;  % 
    192 000 000;  % 
    237 125 049;  % 
    255 255 000;  % 
    146 208 080;  %
    255 192 000
]/255;

colors = parula(44);



% Initialize cell arrays to store velocity data
all_velocities = cell(size(filenum,2),1);
all_times = cell(size(filenum,2),1);

% Process files with suffixes 1 to 6
for i = filenum
    % Construct full filename
    filename = [base_filename_1 num2str(i) base_filename_2 '.tif'];
    filepath = fullfile(base_pathname, filename);
    % filepath = 'E:\Volumetric blood vessel imaging data\testFolder5\testKymo_long.tif';
    
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
    params.frmRate = 50.4*2;         % Frame rate in Hz
    params.pxlSize = 0.8702;       % Pixel size in um
    params.numavgs = 30;          % Number of lines to average
    params.shiftamt = 3;          % Shift amount for correlation
    params.skipamt = 20;           % Skip amount for processing
    params.numWorkers = 12;       % Number of parallel workers
    params.kymo_path = filepath;  % Save the kymograph path for saving results
    
    % Run LSPIV analysis
    [velocity_mm, time, shifted_kymograph] = LSPIV_parallel_simple(kymograph, params);
    
    % Store velocity data
    all_velocities{find(filenum==i)} = velocity_mm;
    all_times{find(filenum==i)} = time;
    
    % Clear some variables to free memory
    clear kymograph shifted_kymograph;
end

% Create combined velocity plot
velocity_fig = figure('Name', 'All Velocity Profiles', 'Position', [200, 200, 350, 200]);
hold on;
box on;
set(gca, 'TickDir', 'out');

% Plot all velocity curves
y_min = inf;
y_max = -inf;
x_min = inf;
x_max = -inf;
mean_velocities = zeros(size(filenum,2),1);

for i = filenum %plotNo%
    if ~isempty(all_velocities{find(filenum==i)})
        plot(all_times{find(filenum==i)}, all_velocities{find(filenum==i)}, 'Color', colors(i,:), 'LineWidth', 1.5);
        y_min = min(y_min, min(all_velocities{find(filenum==i)}));
        y_max = max(y_max, max(all_velocities{find(filenum==i)}));
        x_min = min(x_min, min(all_times{find(filenum==i)}));
        x_max = max(x_max, max(all_times{find(filenum==i)}));
        mean_velocities(find(filenum==i)) = mean(all_velocities{find(filenum==i)}, 'omitnan');
    end
end

% Set axis limits and labels
xlim([x_min x_max]);
y_range = y_max - y_min;
ylim([y_min-0.1*y_range y_max+0.1*y_range]);
xlabel('Time (ms)');
ylabel('Velocity (mm/s)');
%title('Blood Flow Velocity Profiles');

% Add legend with mean velocities
legend_entries = cell(size(filenum,2),1);
for i = filenum %plotNo%
    if ~isempty(all_velocities{find(filenum==i)})
        legend_entries{find(filenum==i)} = sprintf('Vessel %d (%.2f mm/s)', i, mean_velocities(find(filenum==i)));
    else
        legend_entries{find(filenum==i)} = sprintf('Vessel %d (N/A)', i);
    end
end
legend(legend_entries, 'Location', 'best');

% Adjust font
ft_size = 12;
set(gca, 'FontSize', ft_size)
set(get(gca, 'title'), 'FontSize', ft_size, 'FontName', 'AvantGarde','FontWeight','bold');
set(get(gca, 'xlabel'),'FontSize', ft_size, 'FontName', 'AvantGarde','FontWeight','bold');
set(get(gca, 'ylabel'),'FontSize', ft_size, 'FontName', 'AvantGarde','FontWeight','bold');
set(get(gca, 'legend'),'FontSize', ft_size, 'FontName', 'AvantGarde','FontWeight','bold');

% Add mean velocity lines and text
% for i = 1:6
%     if ~isempty(all_velocities{i})
%         % Draw mean velocity line
%         line([x_min x_max], [mean_velocities(i) mean_velocities(i)], ...
%             'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1);
%     end
% end

% % Save the combined velocity plot
% [~, base_name] = fileparts(base_filename);
% saveas(velocity_fig, fullfile(base_pathname, [base_filename_1 num2str(filenum) '_velocity.png']));

% Save velocity data including mean velocities
save(fullfile(base_pathname, [base_name 'velocity_data.mat']), ...
    'all_velocities', 'all_times', 'colors', 'legend_entries', 'mean_velocities');

disp('Processing complete!');