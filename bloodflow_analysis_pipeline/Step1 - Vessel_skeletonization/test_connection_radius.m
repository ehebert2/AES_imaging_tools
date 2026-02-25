function test_connection_radius(skel, spacing, test_radii)
% 测试不同搜索半径的连接效果
% skel: 输入骨架
% spacing: 体素尺寸
% test_radii: 要测试的搜索半径数组，如 [1, 2, 3]

if nargin < 3
    test_radii = [1, 2, 3]; % 默认测试3x3x3, 5x5x5, 7x7x7
end

fprintf('\n=== 测试不同搜索半径的连接效果 ===\n');

% 分析原始骨架
original_CC = bwconncomp(skel, 26);
fprintf('原始骨架连通组件数: %d\n', original_CC.NumObjects);

results = [];

for i = 1:length(test_radii)
    radius = test_radii(i);
    neighborhood_size = 2*radius + 1;
    
    fprintf('\n--- 测试搜索半径 %d (%dx%dx%d 邻域) ---\n', ...
        radius, neighborhood_size, neighborhood_size, neighborhood_size);
    
    % 使用当前半径进行连接
    tic;
    skel_connected = connect_broken_vessels_improved(skel, spacing, 15, radius);
    elapsed_time = toc;
    
    % 分析连接后的结果
    connected_CC = bwconncomp(skel_connected, 26);
    
    % 计算统计信息
    component_reduction = original_CC.NumObjects - connected_CC.NumObjects;
    reduction_percentage = component_reduction / original_CC.NumObjects * 100;
    
    % 计算添加的体素数量
    added_voxels = sum(skel_connected(:)) - sum(skel(:));
    
    fprintf('结果统计:\n');
    fprintf('- 连通组件数: %d -> %d (减少 %d, %.1f%%)\n', ...
        original_CC.NumObjects, connected_CC.NumObjects, component_reduction, reduction_percentage);
    fprintf('- 添加体素数: %d\n', added_voxels);
    fprintf('- 处理时间: %.2f 秒\n', elapsed_time);
    
    % 保存结果
    results = [results; struct(...
        'radius', radius, ...
        'neighborhood_size', neighborhood_size, ...
        'original_components', original_CC.NumObjects, ...
        'connected_components', connected_CC.NumObjects, ...
        'component_reduction', component_reduction, ...
        'reduction_percentage', reduction_percentage, ...
        'added_voxels', added_voxels, ...
        'processing_time', elapsed_time, ...
        'connected_skeleton', skel_connected)];
end

% 显示比较总结
fprintf('\n=== 搜索半径比较总结 ===\n');
fprintf('半径\t邻域大小\t组件减少\t减少率%%\t新增体素\t处理时间(秒)\n');
fprintf('----\t--------\t--------\t-------\t--------\t-----------\n');

for i = 1:length(results)
    r = results(i);
    fprintf('%d\t%dx%dx%d\t%d\t\t%.1f%%\t\t%d\t\t%.2f\n', ...
        r.radius, r.neighborhood_size, r.neighborhood_size, r.neighborhood_size, ...
        r.component_reduction, r.reduction_percentage, r.added_voxels, r.processing_time);
end

% 推荐最佳半径
[~, best_idx] = max([results.reduction_percentage]);
best_radius = results(best_idx).radius;

fprintf('\n推荐搜索半径: %d (%dx%dx%d 邻域)\n', ...
    best_radius, 2*best_radius+1, 2*best_radius+1, 2*best_radius+1);
fprintf('理由: 获得了最高的组件减少率 (%.1f%%)\n', ...
    results(best_idx).reduction_percentage);

% 可视化比较（如果有图形界面）
if exist('figure', 'file')
    try
        figure('Name', '搜索半径效果比较', 'Position', [100, 100, 1000, 600]);
        
        % 子图1: 组件减少数量
        subplot(2, 3, 1);
        bar([results.radius], [results.component_reduction]);
        xlabel('搜索半径');
        ylabel('连通组件减少数量');
        title('组件减少效果');
        grid on;
        
        % 子图2: 减少百分比
        subplot(2, 3, 2);
        bar([results.radius], [results.reduction_percentage]);
        xlabel('搜索半径');
        ylabel('减少百分比 (%)');
        title('组件减少率');
        grid on;
        
        % 子图3: 新增体素数
        subplot(2, 3, 3);
        bar([results.radius], [results.added_voxels]);
        xlabel('搜索半径');
        ylabel('新增体素数');
        title('连接代价');
        grid on;
        
        % 子图4: 处理时间
        subplot(2, 3, 4);
        bar([results.radius], [results.processing_time]);
        xlabel('搜索半径');
        ylabel('处理时间 (秒)');
        title('计算效率');
        grid on;
        
        % 子图5: 效率vs效果散点图
        subplot(2, 3, 5);
        scatter([results.reduction_percentage], [results.processing_time], 100, 'filled');
        xlabel('组件减少率 (%)');
        ylabel('处理时间 (秒)');
        title('效率vs效果');
        grid on;
        
        % 添加半径标签
        for i = 1:length(results)
            text(results(i).reduction_percentage + 0.5, results(i).processing_time + 0.01, ...
                sprintf('R=%d', results(i).radius));
        end
        
        % 子图6: 综合评分
        % 综合评分 = 减少率权重 - 时间权重 - 体素增加权重
        normalize = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);
        norm_reduction = normalize([results.reduction_percentage]);
        norm_time = normalize([results.processing_time]);
        norm_voxels = normalize([results.added_voxels]);
        
        comprehensive_score = norm_reduction - 0.3 * norm_time - 0.2 * norm_voxels;
        
        subplot(2, 3, 6);
        bar([results.radius], comprehensive_score);
        xlabel('搜索半径');
        ylabel('综合评分');
        title('综合评价 (越高越好)');
        grid on;
        
        [~, best_comprehensive_idx] = max(comprehensive_score);
        hold on;
        bar(results(best_comprehensive_idx).radius, comprehensive_score(best_comprehensive_idx), 'r');
        hold off;
        
    catch ME
        fprintf('可视化失败: %s\n', ME.message);
    end
end

% 返回结果供进一步分析
assignin('base', 'connection_test_results', results);
fprintf('\n测试结果已保存到工作空间变量 "connection_test_results"\n');

end 