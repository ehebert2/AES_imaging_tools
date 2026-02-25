function result = morphological_bridge(binary_img, spacing, bridge_distance)
% 使用形态学操作连接断开的血管段
% binary_img: 二值图像
% spacing: 体素尺寸
% bridge_distance: 桥接距离（微米）

if nargin < 3
    bridge_distance = 15;
end

% 转换为像素尺寸
bridge_pixels = round(bridge_distance ./ spacing);

% 多尺度连接策略
result = binary_img;

% 1. 小尺度连接 - 连接小间隙
small_se = strel('sphere', min(bridge_pixels));
temp1 = imclose(result, small_se);

% 2. 中等尺度连接 - 各向异性结构元素
if length(unique(spacing)) > 1  % 各向异性体素
    % 为每个方向创建不同的结构元素
    for dim = 1:3
        se_size = round(bridge_pixels(dim));
        if se_size > 0
            se_dim = zeros(1, 1, 1);
            if dim == 1
                se_dim = strel('line', se_size, 90);  % Y方向
            elseif dim == 2
                se_dim = strel('line', se_size, 0);   % X方向
            else
                se_dim = strel('line', se_size, 0);   % Z方向近似
            end
            temp1 = imclose(temp1, se_dim);
        end
    end
end

% 3. 方向性连接 - 基于主方向的连接
temp2 = directional_closing(temp1, spacing, bridge_pixels);

% 4. 组合结果，保留原始细节
result = temp1 | temp2;

% 5. 可选：保持原始连通性
if sum(result(:)) > sum(binary_img(:)) * 1.5
    % 如果增加太多，使用更保守的方法
    result = binary_img | imclose(binary_img, strel('sphere', 2));
end

end

function result = directional_closing(img, spacing, bridge_pixels)
% 基于主方向的方向性闭运算
result = img;

% 分析图像的主方向
[main_directions, strengths] = analyze_main_directions(img);

% 为每个主方向创建结构元素
for i = 1:size(main_directions, 1)
    if strengths(i) > 0.1  % 只处理显著方向
        direction = main_directions(i, :);
        
        % 创建线性结构元素
        se = create_linear_se(direction, bridge_pixels, spacing);
        
        % 应用闭运算
        temp = imclose(img, se);
        result = result | temp;
    end
end
end

function [directions, strengths] = analyze_main_directions(img)
% 分析图像的主要方向
directions = [1, 0, 0; 0, 1, 0; 0, 0, 1;  % 主轴方向
              1, 1, 0; 1, -1, 0; 1, 0, 1;  % 对角方向
              1, 0, -1; 0, 1, 1; 0, 1, -1];
directions = directions ./ sqrt(sum(directions.^2, 2));

strengths = zeros(size(directions, 1), 1);

% 计算每个方向的强度（简化版本）
for i = 1:size(directions, 1)
    strengths(i) = calculate_direction_strength(img, directions(i, :));
end

% 归一化强度
strengths = strengths / max(strengths);
end

function strength = calculate_direction_strength(img, direction)
% 计算特定方向的强度
strength = 0;

% 使用方向性滤波器估计
kernel_size = 7;
kernel = create_directional_kernel(direction, kernel_size);

% 应用滤波器
filtered = imfilter(double(img), kernel, 'same');
strength = sum(filtered(:).^2) / sum(img(:));
end

function kernel = create_directional_kernel(direction, size)
% 创建方向性滤波器核
center = ceil(size / 2);
kernel = zeros(size, size, size);

% 简化的线性核
for i = 1:size
    for j = 1:size
        for k = 1:size
            pos = [i, j, k] - center;
            proj = abs(dot(pos, direction));
            if proj < 1
                kernel(i, j, k) = 1;
            end
        end
    end
end

kernel = kernel / sum(kernel(:));
end

function se = create_linear_se(direction, bridge_pixels, spacing)
% 创建线性结构元素
max_extent = max(bridge_pixels);
se_size = round(max_extent * 2 + 1);

% 确保奇数尺寸
if mod(se_size, 2) == 0
    se_size = se_size + 1;
end

center = ceil(se_size / 2);
se_matrix = zeros(se_size, se_size, se_size);

% 沿着方向创建线
for t = -max_extent:max_extent
    pos = center + round(t * direction ./ spacing);
    
    % 检查边界
    if all(pos >= 1) && all(pos <= se_size)
        se_matrix(pos(1), pos(2), pos(3)) = 1;
    end
end

se = strel('arbitrary', se_matrix);
end 