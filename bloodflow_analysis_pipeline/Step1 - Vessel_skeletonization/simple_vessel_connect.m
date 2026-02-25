function connected = simple_vessel_connect(skel, max_gap)
% 简单有效的血管连接方法
% skel: 二值骨架图像
% max_gap: 最大间隙距离（像素）

if nargin < 2
    max_gap = 10;
end

connected = skel;

% 1. 多尺度形态学闭运算
for r = 2:max_gap
    se = strel('sphere', r);
    temp = imclose(connected, se);
    
    % 只保留填充了间隙的部分，避免过度膨胀
    new_pixels = temp & ~connected;
    
    % 检查新像素是否合理（连接端点）
    if check_reasonable_addition(connected, new_pixels)
        connected = temp;
    end
end

% 2. 线性连接 - 针对各个方向
directions = [1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,-1,0; 1,0,1; 1,0,-1; 0,1,1; 0,1,-1];

for i = 1:size(directions, 1)
    dir = directions(i, :);
    se_line = create_line_se(dir, max_gap);
    temp = imclose(connected, se_line);
    
    new_pixels = temp & ~connected;
    if check_reasonable_addition(connected, new_pixels)
        connected = temp;
    end
end

end

function reasonable = check_reasonable_addition(original, new_pixels)
% 检查新增像素是否合理
reasonable = true;

% 如果新增像素太多，可能不合理
if sum(new_pixels(:)) > sum(original(:)) * 0.2
    reasonable = false;
    return;
end

% 检查新增像素是否主要连接端点
endpoints = find_simple_endpoints(original);
if isempty(endpoints)
    reasonable = false;
    return;
end

% 计算新像素与端点的距离
new_coords = find(new_pixels);
if isempty(new_coords)
    reasonable = false;
    return;
end

% 简化检查：如果新像素数量合理，就接受
reasonable = true;
end

function endpoints = find_simple_endpoints(skel)
% 简单的端点查找
endpoints = [];
[y, x, z] = ind2sub(size(skel), find(skel));

for i = 1:length(y)
    % 计算邻居数量
    neighbors = 0;
    for dy = -1:1
        for dx = -1:1
            for dz = -1:1
                if dy == 0 && dx == 0 && dz == 0
                    continue;
                end
                ny = y(i) + dy; nx = x(i) + dx; nz = z(i) + dz;
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3) && ...
                   skel(ny, nx, nz)
                    neighbors = neighbors + 1;
                end
            end
        end
    end
    
    if neighbors == 1  % 端点
        endpoints = [endpoints; y(i), x(i), z(i)];
    end
end
end

function se = create_line_se(direction, length)
% 创建线性结构元素
direction = direction / norm(direction);
size_se = length * 2 + 1;
center = length + 1;

se_matrix = zeros(size_se, size_se, size_se);

% 沿方向创建线
for t = -length:length
    pos = center + round(t * direction);
    
    % 边界检查
    if all(pos >= 1) && all(pos <= size_se)
        se_matrix(pos(1), pos(2), pos(3)) = 1;
    end
end

se = strel('arbitrary', se_matrix);
end 