function visualize_components(skel, max_components_to_show)
% 可视化连通组件
% skel: 二值骨架图像
% max_components_to_show: 显示的最大组件数量

if nargin < 2
    max_components_to_show = 10;
end

% 检测连通组件
CC = bwconncomp(skel, 26);
fprintf('总共发现 %d 个连通组件\n', CC.NumObjects);

if CC.NumObjects == 0
    fprintf('没有找到连通组件\n');
    return;
end

% 计算组件大小并排序
component_sizes = cellfun(@numel, CC.PixelIdxList);
[sorted_sizes, sort_idx] = sort(component_sizes, 'descend');

% 显示前几个最大的组件
num_to_show = min(max_components_to_show, CC.NumObjects);

fprintf('\n前 %d 个最大的连通组件:\n', num_to_show);
for i = 1:num_to_show
    comp_idx = sort_idx(i);
    size = sorted_sizes(i);
    percentage = (size / sum(component_sizes)) * 100;
    
    fprintf('%2d. 组件 %d: %d 个体素 (%.1f%%)\n', ...
        i, comp_idx, size, percentage);
end

% 创建彩色标记图像用于可视化
labeled_img = labelmatrix(CC);

% 为了更好的可视化，重新分配标签使大组件有连续的标签
new_labeled_img = zeros(size(labeled_img));
for i = 1:num_to_show
    comp_idx = sort_idx(i);
    new_labeled_img(CC.PixelIdxList{comp_idx}) = i;
end

% 显示统计信息
fprintf('\n连通性统计:\n');
fprintf('- 最大组件占总体积的 %.1f%%\n', (sorted_sizes(1)/sum(component_sizes))*100);
fprintf('- 前3个组件共占总体积的 %.1f%%\n', ...
    (sum(sorted_sizes(1:min(3,length(sorted_sizes))))/sum(component_sizes))*100);

small_components = sum(component_sizes < 50);
fprintf('- 小组件 (<50体素) 数量: %d (%.1f%%)\n', ...
    small_components, (small_components/CC.NumObjects)*100);

% 如果有MATLAB图像处理工具箱，可以创建3D可视化
if exist('isosurface', 'file')
    try
        figure;
        
        % 显示最大的几个组件
        colors = lines(num_to_show);
        
        for i = 1:min(3, num_to_show)  % 只显示前3个最大的
            comp_mask = false(size(skel));
            comp_mask(CC.PixelIdxList{sort_idx(i)}) = true;
            
            % 创建等值面
            [faces, verts] = isosurface(comp_mask, 0.5);
            
            if ~isempty(faces)
                patch('Faces', faces, 'Vertices', verts, ...
                      'FaceColor', colors(i,:), 'EdgeColor', 'none', ...
                      'FaceAlpha', 0.7);
            end
        end
        
        title(sprintf('前 %d 个最大连通组件的3D可视化', min(3, num_to_show)));
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal;
        view(3);
        camlight;
        lighting gouraud;
        legend(arrayfun(@(x) sprintf('组件 %d (%d体素)', x, sorted_sizes(x)), ...
               1:min(3, num_to_show), 'UniformOutput', false));
        
    catch
        fprintf('3D可视化失败，可能需要图像处理工具箱\n');
    end
end

% 返回重新标记的图像用于保存
assignin('caller', 'labeled_components_vis', new_labeled_img);

end 