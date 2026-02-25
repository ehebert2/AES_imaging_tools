function skel_connected = connect_broken_vessels(skel, spacing, max_gap_distance)
% 连接断开的血管段
% skel: 输入骨架（二值图像）
% spacing: 体素尺寸 [y,x,z]
% max_gap_distance: 最大连接距离（微米）

if nargin < 3
    max_gap_distance = 20; % 默认20微米
end

% 转换为像素距离
max_gap_pixels = max_gap_distance ./ spacing;

% 找到所有端点（只有1个邻居的点）
[endpoints, ~] = find_endpoints(skel);

if size(endpoints, 1) < 2
    skel_connected = skel;
    return;
end

% 计算端点间的距离矩阵
distance_matrix = calculate_endpoint_distances(endpoints, spacing);

% 找到需要连接的端点对
connections = find_potential_connections(endpoints, distance_matrix, max_gap_pixels, skel);

% 执行连接
skel_connected = skel;
for i = 1:size(connections, 1)
    p1 = endpoints(connections(i,1), :);
    p2 = endpoints(connections(i,2), :);
    
    % 使用3D直线连接
    connecting_line = connect_3d_line(p1, p2);
    skel_connected(connecting_line) = 1;
end

end

function [endpoints, endpoint_coords] = find_endpoints(skel)
% 找到骨架的端点
se = strel('sphere', 1);
dilated = imdilate(skel, se);
eroded = imerode(skel, se);

% 端点是膨胀后连通但腐蚀后断开的点
endpoints_mask = skel & ~eroded;

% 获取端点坐标
[y, x, z] = ind2sub(size(skel), find(endpoints_mask));
endpoint_coords = [y, x, z];

% 进一步筛选真正的端点（只有1个26连通邻居）
true_endpoints = [];
for i = 1:size(endpoint_coords, 1)
    y = endpoint_coords(i, 1);
    x = endpoint_coords(i, 2);
    z = endpoint_coords(i, 3);
    
    % 检查26连通邻域
    neighbors = 0;
    for dy = -1:1
        for dx = -1:1
            for dz = -1:1
                if dy == 0 && dx == 0 && dz == 0
                    continue;
                end
                ny = y + dy; nx = x + dx; nz = z + dz;
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3)
                    if skel(ny, nx, nz)
                        neighbors = neighbors + 1;
                    end
                end
            end
        end
    end
    
    if neighbors <= 1 % 真正的端点
        true_endpoints = [true_endpoints; y, x, z];
    end
end

endpoints = true_endpoints;
endpoint_coords = endpoints;
end

function distance_matrix = calculate_endpoint_distances(endpoints, spacing)
% 计算端点间的真实距离（考虑体素尺寸）
n = size(endpoints, 1);
distance_matrix = zeros(n, n);

for i = 1:n
    for j = i+1:n
        p1 = endpoints(i, :);
        p2 = endpoints(j, :);
        
        % 计算真实距离（微米）
        real_distance = sqrt(sum(((p1 - p2) .* spacing').^2));
        distance_matrix(i, j) = real_distance;
        distance_matrix(j, i) = real_distance;
    end
end
end

function connections = find_potential_connections(endpoints, distance_matrix, max_gap_pixels, skel)
% 找到应该连接的端点对
n = size(endpoints, 1);
connections = [];

for i = 1:n
    for j = i+1:n
        distance = distance_matrix(i, j);
        
        % 距离筛选
        if distance > max(max_gap_pixels)
            continue;
        end
        
        p1 = endpoints(i, :);
        p2 = endpoints(j, :);
        
        % 方向一致性检查
        if check_direction_consistency(p1, p2, skel)
            connections = [connections; i, j];
        end
    end
end
end

function is_consistent = check_direction_consistency(p1, p2, skel)
% 检查两个端点的方向是否一致（指向对方）
is_consistent = false;

% 获取端点的局部方向
dir1 = get_endpoint_direction(p1, skel);
dir2 = get_endpoint_direction(p2, skel);

if isempty(dir1) || isempty(dir2)
    is_consistent = true; % 如果无法确定方向，默认允许连接
    return;
end

% 计算端点间的向量
connection_vector = p2 - p1;
connection_vector = connection_vector / norm(connection_vector);

% 检查方向一致性
dot1 = dot(dir1, connection_vector);
dot2 = dot(dir2, -connection_vector);

% 如果两个端点都指向对方，则认为方向一致
if dot1 > 0.3 && dot2 > 0.3
    is_consistent = true;
end
end

function direction = get_endpoint_direction(point, skel)
% 获取端点的局部方向
y = point(1); x = point(2); z = point(3);
direction = [];

% 在端点附近找到连接的骨架点
for r = 1:3
    neighbors = [];
    for dy = -r:r
        for dx = -r:r
            for dz = -r:r
                if dy == 0 && dx == 0 && dz == 0
                    continue;
                end
                ny = y + dy; nx = x + dx; nz = z + dz;
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3)
                    if skel(ny, nx, nz)
                        neighbors = [neighbors; ny, nx, nz];
                    end
                end
            end
        end
    end
    
    if ~isempty(neighbors)
        % 计算主方向
        centered = neighbors - point;
        if size(centered, 1) > 1
            [~, ~, V] = svd(centered);
            direction = V(:, 1)';
        else
            direction = centered / norm(centered);
        end
        break;
    end
end
end

function line_indices = connect_3d_line(p1, p2)
% 在3D空间中连接两点的直线
line_indices = [];

% Bresenham 3D线算法的简化版本
steps = max(abs(p2 - p1));
if steps == 0
    return;
end

for i = 0:steps
    t = i / steps;
    point = round(p1 + t * (p2 - p1));
    
    % 转换为线性索引（这里需要知道图像尺寸，简化处理）
    % 实际使用时需要传入图像尺寸
    line_indices = [line_indices; point];
end

% 去除重复点
line_indices = unique(line_indices, 'rows');
end 