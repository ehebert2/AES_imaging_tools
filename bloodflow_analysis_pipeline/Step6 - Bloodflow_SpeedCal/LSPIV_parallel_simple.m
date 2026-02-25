function [velocity_mm, time, shifted_kymograph] = LSPIV_parallel_simple(kymograph, params)
% LSPIV_parallel_simple - Simplified version for blood flow velocity calculation from kymograph
% Input:
%   kymograph - 2D matrix containing the kymograph data
%   params - structure containing parameters:
%       .frmRate - frame rate in Hz (default: 1000)
%       .pxlSize - pixel size in um (default: 0.476)
%       .numavgs - number of lines to average (default: 50)
%       .shiftamt - shift amount for correlation (default: 2)
%       .skipamt - skip amount for processing (default: 2)
%       .numWorkers - number of parallel workers (default: 12)
% Output:
%   velocity_mm - velocity in mm/s
%   time - time points in ms
%   shifted_kymograph - kymograph with each line shifted based on calculated velocity

% Set default parameters if not provided
if nargin < 2
    params = struct();
end
if ~isfield(params, 'frmRate')
    params.frmRate = 1000;
end
if ~isfield(params, 'pxlSize')
    params.pxlSize = 0.476;
end
if ~isfield(params, 'numavgs')
    params.numavgs = 50;
end
if ~isfield(params, 'shiftamt')
    params.shiftamt = 2;
end
if ~isfield(params, 'skipamt')
    params.skipamt = 2;
end
if ~isfield(params, 'numWorkers')
    params.numWorkers = 12;
end

% Limit to first 500 frames (why???)
max_frames = min(5*1e5, size(kymograph, 1));
kymograph = kymograph(1:max_frames, :);

% Start parallel pool
try
    parpool('local', params.numWorkers)
catch
    disp('Parallel pool already exists');
end

% Convert kymograph to double type if it's not already
kymograph = double(kymograph);

% Remove DC offset
disp('DC correction')
DCoffset = mean(kymograph, 1);
imageLinesDC = kymograph - repmat(DCoffset, size(kymograph,1), 1);

% Perform LSPIV correlation
disp('LSPIV calculation...')
scene_fft = fft(imageLinesDC(1:end-params.shiftamt,:), [], 2);
test_img = zeros(size(scene_fft));
test_img(:,:) = imageLinesDC(params.shiftamt+1:end, :);
test_fft = fft(test_img, [], 2);
W = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only

LSPIVresultFFT = scene_fft .* conj(test_fft) .* W;
LSPIVresult = ifft(LSPIVresultFFT, [], 2);

% Calculate velocities
disp('Finding peaks...')
maxpxlshift = round(size(kymograph,2)/2)-1;
index_vals = params.skipamt:params.skipamt:(size(LSPIVresult,1) - params.numavgs);
numpixels = size(LSPIVresult,2);
velocity = nan(size(index_vals));

% Iterate through segments
parfor index = 1:length(index_vals)
    if mod(index_vals(index),100) == 0
        fprintf('Processing line: %d\n', index_vals(index))
    end
    
    % Average correlation results
    LSPIVresult_AVG = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+params.numavgs,:),1)) ...
        / max(sum(LSPIVresult(index_vals(index):index_vals(index)+params.numavgs,:),1));
    
    % Find peak
    c = zeros(1, numpixels);
    c(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift) = ...
        LSPIVresult_AVG(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift);
    [~, maxindex] = max(c);
    
    % Fit gaussian for subpixel accuracy
    options = fitoptions('gauss1');
    options.Lower = [0 numpixels/2-maxpxlshift 0 0];
    options.Upper = [1e9 numpixels/2+maxpxlshift 20 1];
    options.StartPoint = [1 maxindex 10 .1];
    [q,~] = fit((1:length(LSPIVresult_AVG))', LSPIVresult_AVG', 'a1*exp(-((x-b1)/c1)^2) + d1', options);
    
    % Calculate velocity
    velocity(index) = (q.b1 - size(LSPIVresult,2)/2 - 1)/params.shiftamt;
end

% Correct velocity direction
%medV = median(velocity);
%velocity = medV/abs(medV)*velocity;

% Convert to mm/s
time = (index_vals + params.numavgs/2)/params.frmRate*1000; % time in ms
velocity_mm = params.pxlSize/1000*params.frmRate*velocity; % velocity in mm/s

% Apply median filter to smooth results
velocity_mm = medfilt1(velocity_mm, 3);

% Create shifted kymograph
disp('Creating shifted kymograph...')

% Interpolate velocity for all lines
start_line = min(index_vals) + floor(params.numavgs/2);
all_lines = start_line:max_frames;
velocity_valid = interp1(index_vals, velocity, all_lines, 'linear', 'extrap');

% Calculate cumulative displacement in pixels
pixel_shifts = velocity_valid;  % Compensate flow direction  * params.shiftamt/2
cumulative_shift = cumsum(pixel_shifts);

% Calculate maximum shift to determine extended image size
max_shift = max(abs(cumulative_shift));
extended_width = size(kymograph,2) + ceil(max_shift)*2;

% Create extended kymograph
shifted_kymograph = zeros(max_frames, extended_width);

% Copy first part (before velocity calculation starts)
base_pos = ceil(max_shift);  % Center the image
for i = 1:start_line-1
    shifted_kymograph(i, base_pos:base_pos+size(kymograph,2)-1) = kymograph(i,:);
end

% Apply shifts to each line
for i = 1:length(all_lines)
    current_line = all_lines(i);
    shift = round(cumulative_shift(i));
    pos = base_pos + shift;
    
    % Ensure we don't exceed array bounds
    if pos < 1
        pos = 1;
    elseif pos + size(kymograph,2) - 1 > extended_width
        pos = extended_width - size(kymograph,2) + 1;
    end
    
    shifted_kymograph(current_line, pos:pos+size(kymograph,2)-1) = kymograph(current_line,:);
end

% Display results
mean_velocity = mean(velocity_mm, 'omitnan');
hf = figure('Name', 'Velocity Profile', 'Position', [100, 100, 800, 300]);
plot(time, velocity_mm, 'black-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Velocity (mm/s)');
title(sprintf('Blood Flow Velocity (Mean: %.2f mm/s)', mean_velocity));
box on;  % Add box around plot
ax = gca;
ax.XLim = [min(time) max(time)];  % Set x-axis limits to data range
ax.YLim = [min(velocity_mm)-0.1*abs(mean_velocity) max(velocity_mm)+0.1*abs(mean_velocity)];  % Add 10% padding to y-axis
set(gca, 'TickDir', 'out');  % Set tick marks to point outward

% Save velocity plot in the same directory as kymograph
if isfield(params, 'kymo_path')
    [kymo_dir, kymo_name] = fileparts(params.kymo_path);
    % Save velocity plot
    save_path = fullfile(kymo_dir, [kymo_name '_velocity.png']);
    saveas(hf, save_path);
    
    % Save shifted kymograph
    shifted_fig = figure('Name', 'Shifted Kymograph', 'Position', [100, 500, 800, 300]);
    imagesc(shifted_kymograph);
    colormap('gray');
    title(sprintf('Shifted Kymograph (Mean velocity: %.2f mm/s)', mean_velocity));
    xlabel('Position (pixels)');
    ylabel('Time (frames)');
    hold on;
    % Add a line indicating where velocity calculation starts
    line([1 extended_width], [start_line start_line], 'Color', 'r', 'LineStyle', '--');
    text(10, start_line-5, 'Velocity calculation starts', 'Color', 'r');
    colorbar;
    
    % Save shifted kymograph as both figure and data
    Img = uint16(shifted_kymograph);  
    saveas(shifted_fig, fullfile(kymo_dir, [kymo_name '_shifted.png']));
    if numel(Img) >= 2^31-1
        % if the tiff is larger than 4GB, substitute with the following codes
        filename = fullfile(kymo_dir, [kymo_name '_shifted_big.tif']);
        t = Tiff(filename, 'w8');   % BigTIFF

        tagstruct.ImageLength = size(Img,1);
        tagstruct.ImageWidth  = size(Img,2);
        tagstruct.SamplesPerPixel = 1;                % grayscale
        tagstruct.BitsPerSample   = 16;               % uint16
        tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Compression     = Tiff.Compression.None; % or Tiff.Compression.LZW if desired
        tagstruct.RowsPerStrip    = min(size(Img,1), 1024); % reasonable strip size

        t.setTag(tagstruct);
        t.write(Img);
        t.close();
    else
        imwrite(Img, fullfile(kymo_dir, [kymo_name '_shifted.tif']));
    end

end

end 