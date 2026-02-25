function create_connectivity_report(fnrt, CC_vessels, CC_skel_orig, CC_skel_connected, CC_skel_final, component_stats)
% 创建详细的连通性分析报告
% fnrt: 原始文件路径
% CC_*: 各阶段的连通组件信息
% component_stats: 最终组件统计

report_filename = [fnrt(1:end-4) '_connectivity_report.txt'];

fid = fopen(report_filename, 'w');
if fid == -1
    warning('无法创建报告文件');
    return;
end

try
    % 写入报告头
    fprintf(fid, '血管骨架连通性分析报告\n');
    fprintf(fid, '========================\n');
    fprintf(fid, '分析时间: %s\n', datestr(now));
    fprintf(fid, '原始文件: %s\n\n', fnrt);
    
    % 各阶段连通组件数量
    fprintf(fid, '各处理阶段连通组件数量:\n');
    fprintf(fid, '1. 血管分割后: %d 个组件\n', CC_vessels.NumObjects);
    fprintf(fid, '2. 原始骨架: %d 个组件\n', CC_skel_orig.NumObjects);
    fprintf(fid, '3. 连接处理后: %d 个组件\n', CC_skel_connected.NumObjects);
    fprintf(fid, '4. 最终骨架: %d 个组件\n\n', CC_skel_final.NumObjects);
    
    % 连接效果评估
    improvement = CC_skel_orig.NumObjects - CC_skel_connected.NumObjects;
    final_improvement = CC_skel_orig.NumObjects - CC_skel_final.NumObjects;
    
    fprintf(fid, '连接效果评估:\n');
    fprintf(fid, '- 智能连接减少了 %d 个组件 (%.1f%%改善)\n', ...
        improvement, (improvement/CC_skel_orig.NumObjects)*100);
    fprintf(fid, '- 总体减少了 %d 个组件 (%.1f%%改善)\n', ...
        final_improvement, (final_improvement/CC_skel_orig.NumObjects)*100);
    fprintf(fid, '\n');
    
    % 最终组件详细统计
    fprintf(fid, '最终骨架组件统计:\n');
    fprintf(fid, '- 总组件数: %d\n', CC_skel_final.NumObjects);
    fprintf(fid, '- 最大组件: %d 个体素 (%.1f%%)\n', ...
        component_stats.max_size, component_stats.max_percentage);
    fprintf(fid, '- 平均组件大小: %.1f 个体素\n', component_stats.mean_size);
    fprintf(fid, '- 小组件数 (<50体素): %d\n', component_stats.num_small);
    fprintf(fid, '\n');
    
    % 组件大小分布
    sizes = component_stats.sizes;
    fprintf(fid, '组件大小分布:\n');
    fprintf(fid, '- 超大组件 (>1000体素): %d 个\n', sum(sizes > 1000));
    fprintf(fid, '- 大组件 (500-1000体素): %d 个\n', sum(sizes >= 500 & sizes <= 1000));
    fprintf(fid, '- 中等组件 (100-500体素): %d 个\n', sum(sizes >= 100 & sizes < 500));
    fprintf(fid, '- 小组件 (50-100体素): %d 个\n', sum(sizes >= 50 & sizes < 100));
    fprintf(fid, '- 微小组件 (<50体素): %d 个\n', sum(sizes < 50));
    fprintf(fid, '\n');
    
    % 详细组件列表（前10个最大的）
    [sorted_sizes, sort_idx] = sort(sizes, 'descend');
    fprintf(fid, '最大的10个连通组件:\n');
    num_to_show = min(10, length(sorted_sizes));
    for i = 1:num_to_show
        fprintf(fid, '%2d. 组件 %d: %d 个体素 (%.1f%%)\n', ...
            i, sort_idx(i), sorted_sizes(i), ...
            (sorted_sizes(i)/sum(sizes))*100);
    end
    fprintf(fid, '\n');
    
    % 质量评估
    fprintf(fid, '血管网络质量评估:\n');
    
    % 主要网络占比
    main_network_ratio = component_stats.max_percentage;
    if main_network_ratio > 80
        quality = '优秀';
    elseif main_network_ratio > 60
        quality = '良好';
    elseif main_network_ratio > 40
        quality = '一般';
    else
        quality = '需要改进';
    end
    
    fprintf(fid, '- 主要网络完整性: %s (最大组件占%.1f%%)\n', quality, main_network_ratio);
    
    % 碎片化程度
    fragmentation_ratio = component_stats.num_small / CC_skel_final.NumObjects * 100;
    if fragmentation_ratio < 20
        frag_level = '低';
    elseif fragmentation_ratio < 50
        frag_level = '中等';
    else
        frag_level = '高';
    end
    
    fprintf(fid, '- 碎片化程度: %s (%.1f%%为小组件)\n', frag_level, fragmentation_ratio);
    
    % 建议
    fprintf(fid, '\n处理建议:\n');
    if CC_skel_final.NumObjects > 20
        fprintf(fid, '- 连通组件较多，可能需要调整连接参数或检查原始图像质量\n');
    end
    if component_stats.num_small > CC_skel_final.NumObjects * 0.5
        fprintf(fid, '- 小组件较多，建议增加 minBranchLength 参数以去除噪声\n');
    end
    if main_network_ratio < 50
        fprintf(fid, '- 主要网络占比较低，建议检查血管分割质量\n');
    end
    
    fprintf(fid, '\n分析完成。\n');
    
catch ME
    warning('写入报告时出错: %s', ME.message);
end

fclose(fid);
fprintf('连通性分析报告已保存: %s\n', report_filename);

end 