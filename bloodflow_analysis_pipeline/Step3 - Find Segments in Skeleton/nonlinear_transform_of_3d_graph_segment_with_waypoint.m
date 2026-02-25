% Written by Jiang Lan Fan (Modified)
% Visualize blood vessel skeletons and non-linearly transform kymographs
% based on a vessel segment's 3D length with waypoint support.

% load skeleton data
skel = readtiff('sampledata\kf31_0z_roi2_gauss_1umz_1p3238x_skelInterp.tif');
skel = readtiff('sampledata\mouse1_fineStack_Hz_740-675um_6pt5mW_00001_med2_skelInterp.tif');
% x y z from Fiji dictating start and stop position of a vessel segment

a = [215 270 13];b = [237 231 14];
a = [83 265 28];b = [74 180 50];
%a = [88 265 29];b = [69 199 46];

% New: Specify intermediate point c (waypoint)
% If no waypoint is needed, set to empty array []
c = [105 192 45]; % Example: c = [95 220 35]; or c = [];

%a = [217 64 16];b = [183 83 29];
%load and get number of pixels along length of kymograph (from same vessel
%segment)
kymo_file = 'sampledata\mouse1_AES_dualPort_2.tif';
kymo = readtiff(kymo_file);
n_kymo = size(kymo,2);

% convert to list of positions as nodes
[x,y,z] = ind2sub(size(skel),find(skel == 1));
Position = [x y z];

% get rid of edges that connect nodes more than 1.8 pixel apart
d = pdist(Position);
d(d>1.8)=0;

% create graph
dd = squareform(d);
G = graph(dd);
G.Nodes = array2table(Position,'VariableNames',{'X','Y','Z'});

% plot graph
figure
h = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z);
h.NodeColor = [0.5 0.5 0.5]; % Lighter gray to reduce visual interference
h.Marker = 'o'; % Set node marker
h.MarkerSize = 0.5; % Reduce node size
h.EdgeColor = [0.8 0.8 0.8]; % Light gray edges
h.LineWidth = 0.3; % Reduce line width
daspect([1 1 1])

% Set figure properties to improve visual appearance
set(gca, 'Color', 'white'); % White background
grid on; % Add grid
view(3); % 3D view
lighting gouraud; % Add lighting effect

% find the nodes corresponding to segment of interest
loc1 = Position(:, 1) == (a(2)+1) & Position(:, 2) == (a(1)+1) & Position(:, 3) == (a(3)+1);
loc2 = Position(:, 1) == (b(2)+1) & Position(:, 2) == (b(1)+1) & Position(:, 3) == (b(3)+1);
ind1 = find(loc1==1);
ind2 = find(loc2==1);

% New: Handle intermediate point c
if ~isempty(c)
    loc3 = Position(:, 1) == (c(2)+1) & Position(:, 2) == (c(1)+1) & Position(:, 3) == (c(3)+1);
    ind3 = find(loc3==1);
    
    if isempty(ind3)
        error('Waypoint c not found in skeleton, please check coordinates');
    end
    
    % Calculate shortest path from a to c
    [P1,L1] = shortestpath(G,ind1,ind3);
    if isempty(P1)
        error('Cannot find path from start point a to waypoint c');
    end
    
    % Calculate shortest path from c to b
    [P2,L2] = shortestpath(G,ind3,ind2);
    if isempty(P2)
        error('Cannot find path from waypoint c to end point b');
    end
    
    % Concatenate paths (remove duplicate waypoint)
    P = [P1(1:end-1), P2];
    L = L1 + L2;
    
    fprintf('Path information:\n');
    fprintf('- Path length from a to c: %.2f\n', L1);
    fprintf('- Path length from c to b: %.2f\n', L2);
    fprintf('- Total path length: %.2f\n', L);
    
    % Highlight path in graph
    highlight(h,P1,'EdgeColor','r','LineWidth',3); % Increase line width
    highlight(h,P2,'EdgeColor','b','LineWidth',3); % Increase line width
    highlight(h,ind3,'NodeColor','g','MarkerSize',5); % Enlarge waypoint
    highlight(h,ind1,'NodeColor','k','MarkerSize',5); % Enlarge start point
    highlight(h,ind2,'NodeColor','r','MarkerSize',5); % Enlarge end point
    
else
    % Original shortest path calculation
    [P,L] = shortestpath(G,ind1,ind2);
    if isempty(P)
        error('Cannot find path from start point a to end point b');
    end
    
    % Highlight the result in the graph plot
    highlight(h,P,'EdgeColor','r','LineWidth',4); % Increase line width
    highlight(h,ind1,'NodeColor','k','MarkerSize',15); % Enlarge start point
    highlight(h,ind2,'NodeColor','r','MarkerSize',15); % Enlarge end point
    
    fprintf('Shortest path length: %.2f\n', L);
end

alpha(1)
% Display the path length
disp(L)

% remove axis ticks
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

% coords is x,y,z coordinates of the path
coords = Position(P,:);

% 3d projected edge lengths
l_3d = zeros(1, length(coords)-1);
for i=2:length(coords)
    l_3d(i-1) = sqrt((coords(i,1) - coords(i-1,1))^2 + (coords(i,2) - coords(i-1,2))^2 + (coords(i,3) - coords(i-1,3))^2);
end
L_3d = sum(l_3d);
l_3d_cs = cumsum(l_3d);

% map every micron in 3D space to a point in 2D space
real_coords = zeros(floor(L)+1, 3);
real_coords(1,:) = coords(1,:);
% starting at 1st micron (0 um is edge case, corresponds to coords(1,:))
for j=1:floor(L)
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

figure
plot3(coords(:,1), coords(:,2), coords(:,3), 'b-', 'LineWidth',2)
hold on
plot3(real_coords(:,1), real_coords(:,2), real_coords(:,3), 'r.', 'MarkerSize',3)
if ~isempty(c)
    % Mark waypoint - use coordinates from Position array
    plot3(Position(ind3,1), Position(ind3,2), Position(ind3,3), 'go', 'MarkerSize',8, 'MarkerFaceColor','g')
end
plot3(coords(1,1), coords(1,2), coords(1,3), 'ks', 'MarkerSize',8, 'MarkerFaceColor','k') % Start point
plot3(coords(end,1), coords(end,2), coords(end,3), 'rs', 'MarkerSize',8, 'MarkerFaceColor','r') % End point
legend('RawTrace', 'Resample point', 'Way point', 'Start', 'End')
axis equal
grid on

% 2d projected edge lengths from real_coords
l_2d = zeros(1, length(real_coords)-1);
for i=2:length(real_coords)
    l_2d(i-1) = sqrt((real_coords(i,1) - real_coords(i-1,1))^2 + (real_coords(i,2) - real_coords(i-1,2))^2);
end
L_2d = sum(l_2d);
l_2d_cs = cumsum(l_2d);

% convert to index of kymo pixel, per micron in 3D
kymo_ind = round(l_2d_cs * (n_kymo-1) / L_2d) + 1;

% make new, transformed kymo
kymo_transformed = zeros(length(kymo),length(kymo_ind));
for k=1:length(kymo_ind)
    kymo_transformed(:,k) = kymo(:,kymo_ind(k));
end

% Generate output filename
if ~isempty(c)
    output_filename = [kymo_file(1:end-4) '_transformed_with_waypoint.tif'];
else
    output_filename = [kymo_file(1:end-4) '_transformed.tif'];
end

writetiff(uint16(kymo_transformed), output_filename)

fprintf('Transformation complete! Output file: %s\n', output_filename); 