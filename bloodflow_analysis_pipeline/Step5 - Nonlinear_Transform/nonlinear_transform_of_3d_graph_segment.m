% Written by Jiang Lan Fan
% Visualize blood vessel skeletons and non-linearly transform kymographs
% based on a vessel segment's 3D length.

% Create the output directory if it doesn't exist
output_dir_1 = 'E:\Volumetric blood vessel imaging data\testFolder11';
output_dir = [output_dir_1,'\Kymograph_Transformed'];
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% load skeleton data
skel = readtiff('E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\volumeImagingAnalysisCodes\Step1 - Vessel_skeletonization\mouseL4_structure_465-385um_00002_medFilt_skelInterp.tif');
% skel = readtiff('E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\volumeImagingAnalysisCodes\Step1 - Vessel_skeletonization\sampledata\mouseA2_structure_955-875um_00002_medFilt_bgSubt_skelInterp.tif');

% Load segment coordinates from Segments.txt
% segment_coords = readmatrix('Segments.txt');
segment_coords = readmatrix('Segments - mouseL4_400um.txt');

% Get a list of original kymograph files
kymo_files = dir([output_dir_1,'\*segment_*.tif']);

% convert to list of positions as nodes
[x,y,z] = ind2sub(size(skel),find(skel == 1));
Position = [x y z];

% get rid of edges that connect nodes more than 1 pixel apart
d = pdist(Position);
d(d>1.8)=0;

% create graph
dd = squareform(d);
G = graph(dd);
G.Nodes = array2table(Position,'VariableNames',{'X','Y','Z'});

% Loop through each segment
for i = 1:size(segment_coords, 1)
    % Get segment start and end points
    a = segment_coords(i, 1:3);
    b = segment_coords(i, 4:6);

    % Find the corresponding kymograph file
    % This assumes files are named like '..._segment_1_...', '..._segment_2_...', etc.
    current_kymo_file = '';
    for f = 1:length(kymo_files)
        if contains(kymo_files(f).name, ['segment_' num2str(i) '_'])
            current_kymo_file = fullfile(kymo_files(f).folder, kymo_files(f).name);
            break;
        end
    end

    if isempty(current_kymo_file)
        warning('Could not find kymograph for segment %d. Skipping.', i);
        continue;
    end
    
    %load and get number of pixels along length of kymograph (from same vessel
    %segment)
    kymo = readtiff(current_kymo_file);
    n_kymo = size(kymo,2);

    % find the nodes corresponding to segment of interest and calculate shortest path
    loc1 = Position(:, 1) == (a(2)+1) & Position(:, 2) == (a(1)+1) & Position(:, 3) == (a(3)+1);
    loc2 = Position(:, 1) == (b(2)+1) & Position(:, 2) == (b(1)+1) & Position(:, 3) == (b(3)+1);
    ind1 = find(loc1==1);
    ind2 = find(loc2==1);
    [P,L] = shortestpath(G,ind1,ind2);
    
    % If you want to visualize each segment, uncomment the following block
    %{
    figure
    h = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z);
    h.NodeColor = 'none';
    h.EdgeColor = 'k';
    h.LineWidth = 0.5;
    daspect([1 1 1])
    % Highlight the result in the graph plot
    highlight(h,P,'EdgeColor','r','LineWidth',1)
    alpha(1)
    % Display the path length
    disp(L)
    % remove axis ticks
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    %}

    % coords is x,y,z coordinates of shortest path
    coords = Position(P,:);

    % Smooth the path to eliminate artifacts/burrs
    % Using a Gaussian weighted moving average with a window of 5 points.
    % You can adjust the window size for more or less smoothing.
    coords = smoothdata(coords, 1, 'gaussian', 5);

    % Optional: visualize the original vs smoothed path
    
    figure;
    hold on;
    plot3(Position(P,1), Position(P,2), Position(P,3), '.-', 'DisplayName', 'Original Path');
    plot3(coords(:,1), coords(:,2), coords(:,3), '.-', 'DisplayName', 'Smoothed Path');
    hold off;
    axis equal;
    legend;
    title(['Segment ' num2str(i) ' Path Smoothing']);
    %}

    % 3d projected edge lengths
    l_3d = zeros(1, length(coords)-1);
    for j=2:length(coords)
        l_3d(j-1) = sqrt((coords(j,1) - coords(j-1,1))^2 + (coords(j,2) - coords(j-1,2))^2 + (coords(j,3) - coords(j-1,3))^2);
    end
    L_3d = sum(l_3d);
    l_3d_cs = cumsum(l_3d);

    % map every micron in 3D space to a point in 2D space
    real_coords = zeros(floor(L_3d)+1, 3);
    real_coords(1,:) = coords(1,:);
    % starting at 1st micron (0 um is edge case, corresponds to coords(1,:))
    for j=1:floor(L_3d)
        ind = sum(j >= l_3d_cs);
        if ind == 0
            remainder = j;
        else
            remainder = j - l_3d_cs(ind);
        end
        fraction = remainder/l_3d(ind+1);
        real_coords(j+1,1) = coords(ind+1,1) + fraction*(coords(ind+2,1)-coords(ind+1,1));
        real_coords(j+1,2) = coords(ind+1,2) + fraction*(coords(ind+2,2)-coords(ind+1,2));
        real_coords(j+1,3) = coords(ind+1,3) + fraction*(coords(ind+2,3)-coords(ind+1,3));
    end

    % 2d projected edge lengths from real_coords
    l_2d = zeros(1, length(real_coords)-1);
    for j=2:length(real_coords)
        l_2d(j-1) = sqrt((real_coords(j,1) - real_coords(j-1,1))^2 + (real_coords(j,2) - real_coords(j-1,2))^2);
    end
    L_2d = sum(l_2d);
    l_2d_cs = cumsum(l_2d);

    % convert to index of kymo pixel, per micron in 3D
    kymo_ind = round(l_2d_cs * (n_kymo-1) / L_2d) + 1;

    % make new, transformed kymo
    kymo_transformed = zeros(size(kymo,1),length(kymo_ind));
    for k=1:length(kymo_ind)
        kymo_transformed(:,k) = kymo(:,kymo_ind(k));
    end
    
    [~, name, ext] = fileparts(current_kymo_file);
    output_filename = fullfile(output_dir, [name '_transformed' ext]);
    writetiff(uint16(kymo_transformed), output_filename);
    
    fprintf('Processed segment %d and saved to %s\n', i, output_filename);
end
