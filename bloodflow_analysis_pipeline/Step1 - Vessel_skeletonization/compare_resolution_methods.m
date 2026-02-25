function compare_resolution_methods(fnrt, spacing)
% 比较不同分辨率处理方法的效果
% fnrt: 输入图像文件路径
% spacing: 体素尺寸 [y, x, z]

if nargin < 2
    spacing = [0.8702; 0.8702; 1]; % 默认体素尺寸
end

fprintf('=== 分辨率处理方法比较 ===\n');

% 方法1: 原始分辨率处理
fprintf('\n方法1: 原始分辨率处理\n');
fprintf('------------------------------\n');
tic;
skel_original = GU_extractVesselness_Skel(fnrt, 'VoxelSizes', spacing);
time1 = toc;
fprintf('处理时间: %.2f 秒\n', time1);

% 方法2: 自适应处理(不重采样)
fprintf('\n方法2: 自适应处理(不重采样)\n');
fprintf('------------------------------\n');
tic;
skel_adaptive = GU_extractVesselness_Skel_adaptive(fnrt, 'VoxelSizes', spacing, 'autoResample', false);
time2 = toc;
fprintf('处理时间: %.2f 秒\n', time2);

% 方法3: 自适应处理(自动重采样)
fprintf('\n方法3: 自适应处理(自动重采样)\n');
fprintf('------------------------------\n');
tic;
skel_resampled = GU_extractVesselness_Skel_adaptive(fnrt, 'VoxelSizes', spacing, 'autoResample', true);
time3 = toc;
fprintf('处理时间: %.2f 秒\n', time3);

% 比较结果
fprintf('\n=== 结果比较 ===\n');

% 读取生成的骨架文件进行比较
[~, basename] = fileparts(fnrt);

% 原始方法的骨架
skel1 = readtiff([basename '_skel.tif']);
CC1 = bwconncomp(skel1, 26);

% 自适应方法(不重采样)的骨架  
skel2 = readtiff([basename '_skel.tif']); % 最新生成的

% 重采样方法的骨架
if exist([basename '_skel_resampled.tif'], 'file')
    skel3 = readtiff([basename '_skel_resampled.tif']);
    CC3 = bwconncomp(skel3, 26);
else
    skel3 = [];
    CC3 = [];
end

CC2 = bwconncomp(skel2, 26);

fprintf('\n连通性比较:\n');
fprintf('方法1 (原始): %d 个连通组件, %d 个体素\n', CC1.NumObjects, sum(skel1(:)));
fprintf('方法2 (自适应): %d 个连通组件, %d 个体素\n', CC2.NumObjects, sum(skel2(:)));
if ~isempty(CC3)
    fprintf('方法3 (重采样): %d 个连通组件, %d 个体素\n', CC3.NumObjects, sum(skel3(:)));
end

% 分析连通组件大小分布
analyze_size_distribution(CC1, '原始方法');
analyze_size_distribution(CC2, '自适应方法');
if ~isempty(CC3)
    analyze_size_distribution(CC3, '重采样方法');
end

% 性能比较
fprintf('\n性能比较:\n');
fprintf('原始方法: %.2f 秒\n', time1);
fprintf('自适应方法: %.2f 秒 (%.1fx)\n', time2, time2/time1);
if exist('time3', 'var')
    fprintf('重采样方法: %.2f 秒 (%.1fx)\n', time3, time3/time1);
end

% 质量评估建议
fprintf('\n=== 质量评估建议 ===\n');
xy_resolution = mean(spacing(1:2));
z_resolution = spacing(3);
anisotropy_ratio = z_resolution / xy_resolution;

if anisotropy_ratio > 1.5
    fprintf('检测到高各向异性 (比率: %.2f)\n', anisotropy_ratio);
    fprintf('建议:\n');
    fprintf('1. 优先使用重采样方法以获得最佳质量\n');
    fprintf('2. 如果内存限制，使用自适应方法\n');
    fprintf('3. 避免使用原始方法\n');
elseif anisotropy_ratio > 1.2
    fprintf('检测到中等各向异性 (比率: %.2f)\n', anisotropy_ratio);
    fprintf('建议:\n');
    fprintf('1. 自适应方法通常足够\n');
    fprintf('2. 对于高质量要求可考虑重采样\n');
else
    fprintf('体素分辨率良好 (比率: %.2f)\n', anisotropy_ratio);
    fprintf('建议: 任何方法都可以，原始方法最快\n');
end

% 生成可视化比较(如果有合适的工具)
try
    create_comparison_visualization(skel1, skel2, skel3, basename);
catch ME
    fprintf('无法生成可视化比较: %s\n', ME.message);
end

end

function analyze_size_distribution(CC, method_name)
% 分析连通组件大小分布

if isempty(CC) || CC.NumObjects == 0
    fprintf('%s: 无连通组件\n', method_name);
    return;
end

sizes = cellfun(@numel, CC.PixelIdxList);
fprintf('\n%s 组件大小分析:\n', method_name);
fprintf('- 最大组件: %d 体素\n', max(sizes));
fprintf('- 平均组件: %.1f 体素\n', mean(sizes));
fprintf('- 中位数: %.1f 体素\n', median(sizes));
fprintf('- 大组件 (>500): %d 个\n', sum(sizes > 500));
fprintf('- 中组件 (50-500): %d 个\n', sum(sizes >= 50 & sizes <= 500));
fprintf('- 小组件 (<50): %d 个\n', sum(sizes < 50));

end

function create_comparison_visualization(skel1, skel2, skel3, basename)
% 创建比较可视化图像

fprintf('\n创建比较可视化...\n');

% 创建RGB合成图像进行比较
if ~isempty(skel3)
    % 三方法比较
    [ny, nx, nz] = size(skel1);
    mid_slice = round(nz/2);
    
    % 提取中间切片
    slice1 = skel1(:, :, mid_slice);
    slice2 = skel2(:, :, mid_slice);
    slice3 = skel3(:, :, mid_slice);
    
    % 创建RGB合成图像 (R=原始, G=自适应, B=重采样)
    rgb_image = cat(3, slice1, slice2, slice3);
    rgb_image = uint8(rgb_image * 255);
    
    % 保存比较图像
    imwrite(rgb_image, [basename '_comparison_RGB.png']);
    fprintf('已保存RGB比较图像: %s_comparison_RGB.png\n', basename);
else
    % 双方法比较
    [ny, nx, nz] = size(skel1);
    mid_slice = round(nz/2);
    
    slice1 = skel1(:, :, mid_slice);
    slice2 = skel2(:, :, mid_slice);
    
    % 创建对比图像
    comparison = zeros(ny, nx, 3);
    comparison(:, :, 1) = slice1; % 红色: 原始
    comparison(:, :, 2) = slice2; % 绿色: 自适应
    
    imwrite(uint8(comparison * 255), [basename '_comparison_RG.png']);
    fprintf('已保存RG比较图像: %s_comparison_RG.png\n', basename);
end

end 