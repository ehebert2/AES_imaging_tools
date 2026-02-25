% Written by Jiang Lan Fan
% Visualize blood vessel skeletons and non-linearly transform kymographs
% based on a vessel segment's 3D length.

% load skeleton data
skel = readtiff('E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\volumeImagingAnalysisCodes\Step1 - Vessel_skeletonization\mouseL4_structure_465-385um_00002_medFilt_skelInterp.tif');
% x y z from Fiji dictating start and stop position of a vessel segment

% load segments coordinates from text file
% Format: Start x, Start y, Start z, End x, End y, End z
segments_coords = readmatrix('Segments.txt');

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

% plot graph
figure
h = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z);
h.NodeColor = 'none';
h.EdgeColor = 'k';
h.LineWidth = 0.5;
daspect([1 1 1])
hold on; % Hold on to plot all segments

% remove axis ticks
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

all_segments_data = struct('path_nodes', {}, 'path_coords', {}, 'real_coords', {}, 'path_length_voxels', {}, 'path_length_um', {}, 'sampled_points', {});
colors = lines(size(segments_coords, 1)); % Generate different colors for each segment

for i = 1:size(segments_coords, 1)
    a = segments_coords(i, 1:3);
    b = segments_coords(i, 4:6);

    % find the nodes corresponding to segment of interest and calculate shortest path
    loc1 = Position(:, 1) == (a(2)+1) & Position(:, 2) == (a(1)+1) & Position(:, 3) == (a(3)+1);
    loc2 = Position(:, 1) == (b(2)+1) & Position(:, 2) == (b(1)+1) & Position(:, 3) == (b(3)+1);
    ind1 = find(loc1, 1);
    ind2 = find(loc2, 1);

    if isempty(ind1) || isempty(ind2)
        warning('Segment %d: Start or end point not found in skeleton. Skipping.', i);
        continue;
    end

    [P,L] = shortestpath(G,ind1,ind2);

    if isempty(P)
        warning('Segment %d: No path found between start and end points. Skipping.', i);
        continue;
    end

    % Highlight the result in the graph plot
    highlight(h,P,'EdgeColor',colors(i,:),'LineWidth',1.5)
    alpha(1)
    % Display the path length
    fprintf('Segment %d Path Length: %f\n', i, L);

    % coords is x,y,z coordinates of shortest path
    coords = Position(P,:);
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
    if floor(L_3d) > 0
        for j=1:floor(L_3d)
            % Find which segment the j-th micron falls into
            % A segment 'k' is the line between coords(k) and coords(k+1)
            segment_idx = find(j <= l_3d_cs, 1, 'first');
        
            % Distance travelled before reaching the current segment
            if segment_idx == 1
                dist_before = 0;
            else
                dist_before = l_3d_cs(segment_idx - 1);
            end
            
            % Distance into the current segment
            remainder = j - dist_before;
            
            % Length of the current segment
            segment_length = l_3d(segment_idx);
            
            % Avoid division by zero
            if segment_length == 0
                fraction = 0;
            else
                fraction = remainder / segment_length;
            end
        
            % Start and end points of the current segment
            start_node = coords(segment_idx, :);
            end_node = coords(segment_idx + 1, :);
            
            % Interpolate to find the coordinate
            real_coords(j+1, :) = start_node + fraction * (end_node - start_node);
        end
    end
    
    % Uniformly sample 10 points along the segment
    sampled_points_xy = [];
    num_samples = 10;
    total_points = size(real_coords, 1);
    if total_points > 1
        sample_indices = round(linspace(1, total_points, num_samples));
        sampled_points_xy = real_coords(sample_indices, 1:2); % Get only X and Y
    end

    % Store data for this segment
    all_segments_data(i).path_nodes = P;
    all_segments_data(i).path_coords = coords;
    all_segments_data(i).real_coords = real_coords;
    all_segments_data(i).path_length_voxels = L;
    all_segments_data(i).path_length_um = L_3d;
    all_segments_data(i).sampled_points = sampled_points_xy;

end

hold off;

% Save all segments data to a .mat file
save('SegmentsData.mat', 'all_segments_data');

% Save the sampled XY points to a text file
fileID = fopen('SegmentsPrompt.txt', 'w');
if fileID == -1
    error('Cannot open SegmentsPrompt.txt for writing.');
end

for i = 1:numel(all_segments_data)
    fprintf(fileID, '# Segment %d\n', i);
    points = all_segments_data(i).sampled_points;
    if ~isempty(points)
        % Write to file, rounding to one decimal place
        for k = 1:size(points, 1)
            fprintf(fileID, '%.1f,%.1f\n', points(k,1), points(k,2));
        end
    end
end
fclose(fileID);
disp('Sampled XY points saved to SegmentsPrompt.txt');

disp('All segments processed and data saved to SegmentsData.mat');

% The following part is commented out as it was for single segment analysis
%{
figure
plot3(coords(:,1), coords(:,2), coords(:,3))
plot3(real_coords(:,1), real_coords(:,2), real_coords(:,3), '.')
axis equal

% 2d projected edge lengths from real_coords
l_2d = zeros(1, length(real_coords)-1);
for i=2:length(real_coords)
    l_2d(i-1) = sqrt((real_coords(i,1) - real_coords(i-1,1))^2 + (real_coords(i,2) - real_coords(i-1,2))^2);
end
L_2d = sum(l_2d);
l_2d_cs = cumsum(l_2d);
%}

