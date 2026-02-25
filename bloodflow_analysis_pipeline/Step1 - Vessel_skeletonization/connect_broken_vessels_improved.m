function skel_connected = connect_broken_vessels_improved(skel, spacing, max_gap_distance, search_radius)
% 改进版血管连接函数 - 支持扩展邻域搜索
% skel: 输入骨架（二值图像）
% spacing: 体素尺寸 [y,x,z]
% max_gap_distance: 最大连接距离（微米）
% search_radius: 邻域搜索半径（体素），默认为1（3x3x3），可设为2（5x5x5）等

if nargin < 3
    max_gap_distance = 20; % 默认20微米
end

if nargin < 4
    search_radius = 2; % 默认使用5x5x5邻域搜索
end

fprintf('开始连接断开的血管 (搜索半径: %d, 邻域: %dx%dx%d)...\n', ...
    search_radius, 2*search_radius+1, 2*search_radius+1, 2*search_radius+1);

% 找到所有端点 - 使用扩展邻域搜索
endpoints = find_endpoints_extended(skel, search_radius);
fprintf('发现 %d 个端点\n', size(endpoints, 1));

if size(endpoints, 1) < 2
    fprintf('端点数量不足，无需连接\n');
    skel_connected = skel;
    return;
end

% 转换为像素距离
max_gap_pixels = max_gap_distance / min(spacing);

% 计算端点间的距离
distance_matrix = calculate_distances(endpoints, spacing);

% 找到需要连接的端点对 - 使用扩展搜索
connections = find_connections_extended(endpoints, distance_matrix, max_gap_pixels, skel, search_radius);
fprintf('找到 %d 个潜在连接\n', size(connections, 1));

% 执行连接
skel_connected = skel;
img_size = size(skel);

if ~isempty(connections)
    for i = 1:size(connections, 1)
        p1 = endpoints(connections(i,1), :);
        p2 = endpoints(connections(i,2), :);
        
        % 使用3D直线连接
        line_points = bresenham_3d(p1, p2);
        
        % 将点添加到骨架中
        for j = 1:size(line_points, 1)
            y = line_points(j, 1);
            x = line_points(j, 2);
            z = line_points(j, 3);
            
            if y >= 1 && y <= img_size(1) && ...
               x >= 1 && x <= img_size(2) && ...
               z >= 1 && z <= img_size(3)
                skel_connected(y, x, z) = 1;
            end
        end
    end
else
    fprintf('未找到合适的连接点\n');
end

fprintf('连接完成！\n');
fprintf('- 发现端点: %d 个\n', size(endpoints, 1));
fprintf('- 建立连接: %d 条\n', size(connections, 1));

% 只有在有距离数据时才显示距离范围
if any(distance_matrix(:) > 0)
    valid_distances = distance_matrix(distance_matrix > 0);
    fprintf('- 连接距离范围: %.1f - %.1f 微米\n', min(valid_distances), max(valid_distances));
end

% 计算连接效果统计
original_components = bwconncomp(skel, 26);
connected_components = bwconncomp(skel_connected, 26);
fprintf('- 连通组件数: %d -> %d (减少了 %d)\n', ...
    original_components.NumObjects, connected_components.NumObjects, ...
    original_components.NumObjects - connected_components.NumObjects);

end

function endpoints = find_endpoints_extended(skel, search_radius)
% 扩展的端点查找函数 - 支持更大的邻域搜索
endpoints = [];
[y_coords, x_coords, z_coords] = ind2sub(size(skel), find(skel));

fprintf('使用扩展邻域搜索 (半径=%d) 查找端点...\n', search_radius);

for i = 1:length(y_coords)
    y = y_coords(i);
    x = x_coords(i);
    z = z_coords(i);
    
    % 计算扩展邻域中的骨架点数量
    neighbor_count = 0;
    neighbor_positions = [];
    
    for dy = -search_radius:search_radius
        for dx = -search_radius:search_radius
            for dz = -search_radius:search_radius
                if dy == 0 && dx == 0 && dz == 0
                    continue;
                end
                
                ny = y + dy;
                nx = x + dx;
                nz = z + dz;
                
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3) && ...
                   skel(ny, nx, nz)
                    neighbor_count = neighbor_count + 1;
                    neighbor_positions = [neighbor_positions; ny, nx, nz];
                end
            end
        end
    end
    
    % 端点定义：在扩展邻域内只有少量邻居
    % 对于5x5x5邻域，如果邻居数≤3，认为是端点
    max_neighbors_for_endpoint = max(1, floor(search_radius * 1.5));
    
    if neighbor_count <= max_neighbors_for_endpoint
        endpoints = [endpoints; y, x, z];
    end
end

fprintf('在扩展邻域中找到 %d 个潜在端点\n', size(endpoints, 1));
end

function distance_matrix = calculate_distances(endpoints, spacing)
% 计算端点间的欧氏距离（考虑体素尺寸）
n = size(endpoints, 1);
distance_matrix = zeros(n, n);

for i = 1:n
    for j = i+1:n
        p1 = endpoints(i, :);
        p2 = endpoints(j, :);
        
        % 计算真实距离（微米）
        diff = (p1 - p2) .* spacing';
        real_distance = sqrt(sum(diff.^2));
        
        distance_matrix(i, j) = real_distance;
        distance_matrix(j, i) = real_distance;
    end
end
end

function connections = find_connections_extended(endpoints, distance_matrix, max_gap_pixels, skel, search_radius)
% 找到应该连接的端点对 - 使用扩展搜索
n = size(endpoints, 1);
connections = [];

% 对每个端点，找到最近的几个候选端点
for i = 1:n
    % 获取当前端点到所有其他端点的距离
    distances = distance_matrix(i, :);
    distances(i) = inf; % 排除自己
    
    % 找到距离最近的端点
    [sorted_dist, sorted_idx] = sort(distances);
    
    % 考虑更多候选点，因为搜索范围扩大了
    max_candidates = min(5, length(sorted_idx)); % 增加到5个候选
    
    for k = 1:max_candidates
        j = sorted_idx(k);
        distance = sorted_dist(k);
        
        % 距离筛选 - 更宽松的距离控制
        if distance > max_gap_pixels * 1.5  % 扩大搜索范围
            break; % 后面的距离只会更大
        end
        
        p1 = endpoints(i, :);
        p2 = endpoints(j, :);
        
        % 检查是否已经有这个连接（避免重复）
        if ~isempty(connections) && any(connections(:,1) == j & connections(:,2) == i)
            continue;
        end
        
        % 检查路径是否合理 - 使用扩展搜索
        if is_reasonable_connection_extended(p1, p2, skel, distance, search_radius)
            connections = [connections; i, j];
            break; % 每个端点只连接一个最佳匹配
        end
    end
end
end

function is_reasonable = is_reasonable_connection_extended(p1, p2, skel, distance, search_radius)
% 改进的连接合理性检查 - 使用扩展搜索
is_reasonable = false;

% 1. 基本距离检查 - 更宽松
max_distance = 40 + search_radius * 5; % 基于搜索半径调整最大距离
if distance > max_distance
    return;
end

% 2. 计算连接路径
line_points = bresenham_3d(p1, p2);

% 3. 检查路径上现有骨架点的密度
existing_points = 0;
total_points = size(line_points, 1);

for i = 1:total_points
    y = line_points(i, 1);
    x = line_points(i, 2);
    z = line_points(i, 3);
    
    if y >= 1 && y <= size(skel,1) && ...
       x >= 1 && x <= size(skel,2) && ...
       z >= 1 && z <= size(skel,3) && ...
       skel(y, x, z)
        existing_points = existing_points + 1;
    end
end

% 4. 路径密度检查 - 更宽松的标准
existing_ratio = existing_points / total_points;
max_existing_ratio = 0.3 + search_radius * 0.1; % 基于搜索半径调整
if existing_ratio > max_existing_ratio
    return;  % 路径上已有太多点，不连接
end

% 5. 扩展的端点方向一致性检查
dir1 = get_endpoint_direction_extended(p1, skel, search_radius);
dir2 = get_endpoint_direction_extended(p2, skel, search_radius);

if ~isempty(dir1) && ~isempty(dir2)
    connection_vector = (p2 - p1) / norm(p2 - p1);
    
    % 检查两个端点是否大致指向对方
    dot1 = dot(dir1, connection_vector);
    dot2 = dot(dir2, -connection_vector);
    
    % 更宽松的角度判断
    angle_threshold = 0.1 - search_radius * 0.02; % 搜索范围越大，角度要求越宽松
    if dot1 > angle_threshold || dot2 > angle_threshold
        is_reasonable = true;
    end
else
    % 如果无法确定方向，基于距离和路径密度决定 - 更宽松
    distance_threshold = 20 + search_radius * 5;
    density_threshold = 0.15 + search_radius * 0.05;
    if distance < distance_threshold && existing_ratio < density_threshold
        is_reasonable = true;
    end
end
end

function direction = get_endpoint_direction_extended(point, skel, search_radius)
% 扩展的端点方向计算
direction = [];
y = point(1); x = point(2); z = point(3);

% 在端点的扩展邻域内寻找连接的骨架点
connected_points = [];

% 使用扩展的搜索半径
for r = 1:search_radius+1  % 比端点检测稍大的搜索范围
    for dy = -r:r
        for dx = -r:r
            for dz = -r:r
                if abs(dy) + abs(dx) + abs(dz) > r || (dy == 0 && dx == 0 && dz == 0)
                    continue;
                end
                
                ny = y + dy; nx = x + dx; nz = z + dz;
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3) && ...
                   skel(ny, nx, nz)
                    connected_points = [connected_points; ny, nx, nz];
                end
            end
        end
    end
    
    % 如果找到足够的点，计算主方向
    if size(connected_points, 1) >= 3  % 需要更多点来可靠计算方向
        break;
    end
end

if size(connected_points, 1) >= 3
    % 使用PCA计算主方向
    centered = connected_points - point;
    [~, ~, V] = svd(centered, 'econ');
    direction = V(:, 1)';
elseif size(connected_points, 1) >= 1
    % 如果点不够多，使用质心方向
    centroid = mean(connected_points, 1);
    direction = (centroid - point) / norm(centroid - point);
end
end

function line_points = bresenham_3d(p1, p2)
% 3D Bresenham算法实现
line_points = [];

% 计算方向和步长
dx = abs(p2(1) - p1(1));
dy = abs(p2(2) - p1(2));
dz = abs(p2(3) - p1(3));

% 确定步长
steps = max([dx, dy, dz]);

if steps == 0
    line_points = p1;
    return;
end

% 计算增量
x_inc = (p2(1) - p1(1)) / steps;
y_inc = (p2(2) - p1(2)) / steps;
z_inc = (p2(3) - p1(3)) / steps;

% 生成路径点
current_point = p1;
for i = 0:steps
    rounded_point = round(current_point);
    line_points = [line_points; rounded_point];
    current_point = current_point + [x_inc, y_inc, z_inc];
end

% 去除重复点
line_points = unique(line_points, 'rows');
end 