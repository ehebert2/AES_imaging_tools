function [degree_stats, node_degrees] = analyze_node_degrees(A, node, link)
% 分析血管网络中各个节点的度数分布
% 
% 输入:
%   A - 邻接矩阵 (来自 Skel2Graph3D)
%   node - 节点结构体数组 (来自 Skel2Graph3D)  
%   link - 连接结构体数组 (来自 Skel2Graph3D)
%
% 输出:
%   degree_stats - 度数统计信息结构体
%   node_degrees - 每个节点的度数数组

fprintf('\n=== 节点度数分析 ===\n');

% 方法1: 使用邻接矩阵计算度数
if nargin >= 1 && ~isempty(A)
    % 对于无向图，度数 = 邻接矩阵每行的非零元素个数
    node_degrees = full(sum(A > 0, 2));
    fprintf('使用邻接矩阵计算节点度数\n');
    
% 方法2: 使用节点结构体计算度数
elseif nargin >= 2 && ~isempty(node)
    node_degrees = zeros(length(node), 1);
    for i = 1:length(node)
        % 计算有效连接数量
        valid_links = node(i).links(node(i).links > 0);
        node_degrees(i) = length(valid_links);
    end
    fprintf('使用节点结构体计算节点度数\n');
    
else
    error('需要提供邻接矩阵A或节点结构体node');
end

% 统计分析
total_nodes = length(node_degrees);
min_degree = min(node_degrees);
max_degree = max(node_degrees);
mean_degree = mean(node_degrees);
median_degree = median(node_degrees);

% 度数分布统计
degree_counts = histcounts(node_degrees, 0:max_degree+1);
degree_values = 0:max_degree;
degree_distribution = degree_counts / total_nodes * 100;

% 特殊节点类型统计
isolated_nodes = sum(node_degrees == 0);  % 孤立节点
end_nodes = sum(node_degrees == 1);       % 端点
junction_nodes = sum(node_degrees >= 3);  % 分叉点
hub_nodes = sum(node_degrees >= 5);       % 枢纽节点

fprintf('\n节点度数统计:\n');
fprintf('- 总节点数: %d\n', total_nodes);
fprintf('- 最小度数: %d\n', min_degree);
fprintf('- 最大度数: %d\n', max_degree);
fprintf('- 平均度数: %.2f\n', mean_degree);
fprintf('- 中位数度数: %.1f\n', median_degree);

fprintf('\n节点类型分布:\n');
fprintf('- 孤立节点 (度数=0): %d (%.1f%%)\n', isolated_nodes, isolated_nodes/total_nodes*100);
fprintf('- 端点 (度数=1): %d (%.1f%%)\n', end_nodes, end_nodes/total_nodes*100);
fprintf('- 普通节点 (度数=2): %d (%.1f%%)\n', sum(node_degrees == 2), sum(node_degrees == 2)/total_nodes*100);
fprintf('- 分叉点 (度数≥3): %d (%.1f%%)\n', junction_nodes, junction_nodes/total_nodes*100);
fprintf('- 枢纽节点 (度数≥5): %d (%.1f%%)\n', hub_nodes, hub_nodes/total_nodes*100);

fprintf('\n度数分布:\n');
for i = 1:min(10, length(degree_values))  % 只显示前10个度数值
    if degree_counts(i) > 0
        fprintf('- 度数 %d: %d 个节点 (%.1f%%)\n', ...
            degree_values(i), degree_counts(i), degree_distribution(i));
    end
end

% 组装统计结果
degree_stats = struct(...
    'total_nodes', total_nodes, ...
    'min_degree', min_degree, ...
    'max_degree', max_degree, ...
    'mean_degree', mean_degree, ...
    'median_degree', median_degree, ...
    'isolated_nodes', isolated_nodes, ...
    'end_nodes', end_nodes, ...
    'junction_nodes', junction_nodes, ...
    'hub_nodes', hub_nodes, ...
    'degree_values', degree_values, ...
    'degree_counts', degree_counts, ...
    'degree_distribution', degree_distribution);

% 可视化度数分布
if exist('figure', 'file')
    try
        figure('Name', '节点度数分布', 'Position', [100, 100, 800, 600]);
        
        % 子图1: 度数分布直方图
        subplot(2, 2, 1);
        bar(degree_values(degree_counts > 0), degree_counts(degree_counts > 0));
        xlabel('节点度数');
        ylabel('节点数量');
        title('节点度数分布');
        grid on;
        
        % 子图2: 度数分布百分比
        subplot(2, 2, 2);
        pie_data = [isolated_nodes, end_nodes, sum(node_degrees == 2), junction_nodes];
        pie_labels = {'孤立节点', '端点', '普通节点', '分叉点'};
        pie_data = pie_data(pie_data > 0);
        pie_labels = pie_labels(pie_data > 0);
        pie(pie_data, pie_labels);
        title('节点类型分布');
        
        % 子图3: 累积分布
        subplot(2, 2, 3);
        cumulative_dist = cumsum(degree_distribution);
        plot(degree_values, cumulative_dist, 'b-o', 'LineWidth', 2);
        xlabel('节点度数');
        ylabel('累积百分比 (%)');
        title('度数累积分布');
        grid on;
        
        % 子图4: 度数vs节点ID散点图
        subplot(2, 2, 4);
        scatter(1:total_nodes, node_degrees, 'filled');
        xlabel('节点ID');
        ylabel('节点度数');
        title('各节点度数分布');
        grid on;
        
        % 高亮高度数节点
        high_degree_threshold = max(3, round(mean_degree + 2*std(node_degrees)));
        high_degree_nodes = find(node_degrees >= high_degree_threshold);
        if ~isempty(high_degree_nodes)
            hold on;
            scatter(high_degree_nodes, node_degrees(high_degree_nodes), ...
                'red', 'filled', 'SizeData', 100);
            legend('所有节点', sprintf('高度数节点 (≥%d)', high_degree_threshold), ...
                'Location', 'best');
            hold off;
        end
        
    catch ME
        fprintf('可视化失败: %s\n', ME.message);
    end
end

% 找出度数最高的节点
if max_degree > 0
    [~, highest_degree_nodes] = sort(node_degrees, 'descend');
    fprintf('\n度数最高的节点 (前5个):\n');
    for i = 1:min(5, length(highest_degree_nodes))
        node_idx = highest_degree_nodes(i);
        fprintf('- 节点 %d: 度数 %d\n', node_idx, node_degrees(node_idx));
        if nargin >= 2 && ~isempty(node)
            fprintf('  位置: (%.1f, %.1f, %.1f)\n', ...
                node(node_idx).comx, node(node_idx).comy, node(node_idx).comz);
        end
    end
end

fprintf('\n分析完成！\n');

end 