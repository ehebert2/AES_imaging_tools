function analyze_resolution_effects(spacing, skel)
% 分析低Z分辨率对血管骨架分析的影响
% spacing: [y, x, z] 体素尺寸 (微米)
% skel: 骨架图像

fprintf('\n=== 分辨率影响分析 ===\n');

% 1. 分析体素各向异性
xy_resolution = mean(spacing(1:2));
z_resolution = spacing(3);
anisotropy_ratio = z_resolution / xy_resolution;

fprintf('体素尺寸分析:\n');
fprintf('- XY平面分辨率: %.3f μm\n', xy_resolution);
fprintf('- Z轴分辨率: %.3f μm\n', z_resolution);
fprintf('- 各向异性比率: %.2f (1.0为理想)\n', anisotropy_ratio);

if anisotropy_ratio > 1.5
    fprintf('⚠️  警告: Z分辨率明显低于XY分辨率 (>1.5倍)\n');
    fprintf('   影响: 可能导致Z方向连接断裂和角度估算不准确\n');
elseif anisotropy_ratio > 1.2
    fprintf('⚠️  注意: Z分辨率略低于XY分辨率\n');
    fprintf('   影响: 轻微影响3D骨架质量\n');
else
    fprintf('✅ 体素分辨率相对均匀\n');
end

% 2. 分析骨架的方向分布
if nargin > 1 && ~isempty(skel)
    fprintf('\n骨架方向分析:\n');
    analyze_skeleton_directions(skel, spacing);
end

% 3. 提供改进建议
fprintf('\n=== 改进建议 ===\n');

if anisotropy_ratio > 1.5
    fprintf('1. 体素插值重采样:\n');
    fprintf('   - 建议将Z方向插值到与XY相同分辨率\n');
    fprintf('   - 使用三次样条或线性插值\n');
    
    fprintf('2. 调整算法参数:\n');
    fprintf('   - 增加Z方向搜索范围\n');
    fprintf('   - 使用各向异性的形态学核\n');
    fprintf('   - 调整vesselness滤波器的sigma值\n');
    
    fprintf('3. 连接策略调整:\n');
    fprintf('   - Z方向使用更大的连接距离阈值\n');
    fprintf('   - 优先连接Z方向的断裂\n');
end

% 4. 计算建议的插值参数
target_resolution = xy_resolution;
z_scale_factor = z_resolution / target_resolution;

fprintf('\n插值参数建议:\n');
fprintf('- 目标等方体素尺寸: %.3f μm\n', target_resolution);
fprintf('- Z方向缩放因子: %.2f\n', z_scale_factor);
fprintf('- 新的图像尺寸倍数: [1, 1, %.2f]\n', z_scale_factor);

end

function analyze_skeleton_directions(skel, spacing)
% 分析骨架的主要方向分布

% 找到所有骨架点
[y, x, z] = ind2sub(size(skel), find(skel));

if length(y) < 100
    fprintf('骨架点太少，无法进行方向分析\n');
    return;
end

% 随机采样部分点进行方向分析
if length(y) > 1000
    sample_idx = randperm(length(y), 1000);
    y = y(sample_idx);
    x = x(sample_idx);
    z = z(sample_idx);
end

directions = [];

% 对每个点，找到其邻域内的方向
for i = 1:length(y)
    center = [y(i), x(i), z(i)];
    neighbors = [];
    
    % 在3x3x3邻域内查找相邻点
    for dy = -1:1
        for dx = -1:1
            for dz = -1:1
                if dy == 0 && dx == 0 && dz == 0
                    continue;
                end
                
                ny = y(i) + dy;
                nx = x(i) + dx;
                nz = z(i) + dz;
                
                if ny >= 1 && ny <= size(skel,1) && ...
                   nx >= 1 && nx <= size(skel,2) && ...
                   nz >= 1 && nz <= size(skel,3) && ...
                   skel(ny, nx, nz)
                    neighbors = [neighbors; ny, nx, nz];
                end
            end
        end
    end
    
    % 如果有足够的邻居，计算主方向
    if size(neighbors, 1) >= 2
        % 转换为真实坐标（考虑spacing）
        center_real = center .* spacing';
        neighbors_real = neighbors .* spacing';
        
        % 计算相对位置向量
        rel_positions = neighbors_real - center_real;
        
        % 使用PCA找主方向
        if size(rel_positions, 1) >= 2
            [~, ~, V] = svd(rel_positions, 'econ');
            main_direction = V(:, 1)';
            directions = [directions; abs(main_direction)]; % 取绝对值避免方向歧义
        end
    end
end

if ~isempty(directions)
    % 分析方向分布
    mean_dir = mean(directions, 1);
    mean_dir = mean_dir / norm(mean_dir); % 归一化
    
    fprintf('主方向分布 (归一化):\n');
    fprintf('- Y方向 (行): %.3f\n', mean_dir(1));
    fprintf('- X方向 (列): %.3f\n', mean_dir(2));
    fprintf('- Z方向 (层): %.3f\n', mean_dir(3));
    
    % 检查是否有明显的偏向
    [max_component, max_idx] = max(mean_dir);
    direction_names = {'Y', 'X', 'Z'};
    
    if max_component > 0.6
        fprintf('⚠️  血管主要沿%s方向分布 (%.1f%%)\n', ...
            direction_names{max_idx}, max_component * 100);
        if max_idx == 3
            fprintf('   注意: 主要沿Z方向可能加剧低分辨率影响\n');
        end
    else
        fprintf('✅ 血管方向分布相对均匀\n');
    end
    
    % 分析Z方向偏好
    z_preference = mean_dir(3);
    if z_preference > 0.4
        fprintf('⚠️  血管有较强的Z方向成分 (%.1f%%)\n', z_preference * 100);
        fprintf('   建议: 特别注意Z方向的连接质量\n');
    end
end

end 