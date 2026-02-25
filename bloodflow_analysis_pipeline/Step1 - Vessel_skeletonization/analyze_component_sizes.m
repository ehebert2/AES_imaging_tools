function [stats, largest_component] = analyze_component_sizes(CC)
% 分析连通组件的大小统计
% CC: bwconncomp的输出结果
% stats: 统计信息结构体
% largest_component: 最大连通组件的信息

if CC.NumObjects == 0
    stats = struct('max_size', 0, 'max_percentage', 0, 'mean_size', 0, ...
                   'small_threshold', 50, 'num_small', 0, 'sizes', []);
    largest_component = struct('idx', [], 'size', 0);
    return;
end

% 计算每个连通组件的大小
component_sizes = cellfun(@numel, CC.PixelIdxList);

% 基本统计
max_size = max(component_sizes);
mean_size = mean(component_sizes);
total_voxels = sum(component_sizes);

% 百分比
max_percentage = (max_size / total_voxels) * 100;

% 小组件统计（小于50个体素认为是小组件）
small_threshold = 50;
num_small = sum(component_sizes < small_threshold);

% 找到最大组件
[~, max_idx] = max(component_sizes);
largest_component = struct('idx', CC.PixelIdxList{max_idx}, 'size', max_size);

% 组装统计结果
stats = struct(...
    'max_size', max_size, ...
    'max_percentage', max_percentage, ...
    'mean_size', mean_size, ...
    'small_threshold', small_threshold, ...
    'num_small', num_small, ...
    'sizes', component_sizes);

% 显示大小分布
fprintf('连通组件大小分布:\n');
fprintf('- 大组件 (>500体素): %d 个\n', sum(component_sizes > 500));
fprintf('- 中等组件 (50-500体素): %d 个\n', sum(component_sizes >= 50 & component_sizes <= 500));
fprintf('- 小组件 (<50体素): %d 个\n', sum(component_sizes < 50));

end 