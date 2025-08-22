clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))
addpath(fullfile(toolbox_path, 'Violin_plot'))
addpath(fullfile(toolbox_path, 'Spider_plot_github_repo'))
addpath(fullfile(toolbox_path, 'Silent_FCD'))

FCD_target_folder = '/path/ to/ target subjects';

rm_outlier_flag = 1;
CSF_dil_flag = 1;
ori_prop_flag = 0; % 1 for non-normalized processing; 0 for normalized processing

if ori_prop_flag == 1
    out_postfix = 'ori';
else
    out_postfix = '';
end

if CSF_dil_flag == 1
    if ori_prop_flag == 1
        out_postfix = [out_postfix '_' 'CSFdil'];
    else
        out_postfix = 'CSFdil';
    end
else
    out_postfix = out_postfix;
end

%%

% data from ROI 1 of target 1
load(fullfile(FCD_target_folder, ['Target_MRF_T1_ROI_GM_voxels_' out_postfix '.mat']))
load(fullfile(FCD_target_folder, ['Target_MRF_T2_ROI_GM_voxels_' out_postfix '.mat']))
load(fullfile(FCD_target_folder, ['Target_MRF_T1_ROI_WM_voxels_' out_postfix '.mat']))
load(fullfile(FCD_target_folder, ['Target_MRF_T2_ROI_WM_voxels_' out_postfix '.mat']))

if rm_outlier_flag == 1
    Target_1_MRF_T1_ROI_GM_voxels = rmoutliers(Target_MRF_T1_ROI_GM_voxels);
    Target_1_MRF_T2_ROI_GM_voxels = rmoutliers(Target_MRF_T2_ROI_GM_voxels);
    Target_1_MRF_T1_ROI_WM_voxels = rmoutliers(Target_MRF_T1_ROI_WM_voxels);
    Target_1_MRF_T2_ROI_WM_voxels = rmoutliers(Target_MRF_T2_ROI_WM_voxels);
else
    Target_1_MRF_T1_ROI_GM_voxels = Target_MRF_T1_ROI_GM_voxels;
    Target_1_MRF_T2_ROI_GM_voxels = Target_MRF_T2_ROI_GM_voxels;
    Target_1_MRF_T1_ROI_WM_voxels = Target_MRF_T1_ROI_WM_voxels;
    Target_1_MRF_T2_ROI_WM_voxels = Target_MRF_T2_ROI_WM_voxels;
end

% data from different ROI or group...
Target_2_MRF_T1_ROI_GM_voxels;
Target_2_MRF_T2_ROI_GM_voxels;
Target_2_MRF_T1_ROI_WM_voxels;
Target_2_MRF_T2_ROI_WM_voxels;

%% norm test 
ft_groups = {Target_1_MRF_T1_ROI_GM_voxels, ...
    Target_1_MRF_T2_ROI_GM_voxels, ...
    Target_1_MRF_T1_ROI_WM_voxels, ...
    Target_1_MRF_T2_ROI_WM_voxels, ...
    ...
    Target_2_MRF_T1_ROI_GM_voxels, ...
    Target_2_MRF_T2_ROI_GM_voxels, ...
    Target_2_MRF_T1_ROI_WM_voxels, ...
    Target_2_MRF_T2_ROI_WM_voxels};

for order = 1:length(ft_groups)
    [h, p] = lillietest(ft_groups{order});
    P_sets(order, 1) = p;
    clear h p
end

%% statistical analysis(if not normal distribution)

[MRF_T1_ROI_GM_p, MRF_T1_ROI_GM_h, MRF_T1_ROI_GM_stats] = ranksum(Target_1_MRF_T1_ROI_GM_voxels, Target_2_MRF_T1_ROI_GM_voxels);
[MRF_T2_ROI_GM_p, MRF_T2_ROI_GM_h, MRF_T2_ROI_GM_stats] = ranksum(Target_1_MRF_T2_ROI_GM_voxels, Target_2_MRF_T2_ROI_GM_voxels);
[MRF_T1_ROI_WM_p, MRF_T1_ROI_WM_h, MRF_T1_ROI_WM_stats] = ranksum(Target_1_MRF_T1_ROI_WM_voxels, Target_2_MRF_T1_ROI_WM_voxels);
[MRF_T2_ROI_WM_p, MRF_T2_ROI_WM_h, MRF_T2_ROI_WM_stats] = ranksum(Target_1_MRF_T2_ROI_WM_voxels, Target_2_MRF_T2_ROI_WM_voxels);

P_values = [MRF_T1_ROI_GM_p; MRF_T2_ROI_GM_p; MRF_T1_ROI_WM_p; MRF_T2_ROI_WM_p];
Z_values = [MRF_T1_ROI_GM_stats.zval; MRF_T2_ROI_GM_stats.zval; MRF_T1_ROI_WM_stats.zval; MRF_T2_ROI_WM_stats.zval];

%% two-sample student t-test (if normal distribution)
[MRF_T1_ROI_GM_h, MRF_T1_ROI_GM_p, MRF_T1_ROI_GM_ci, MRF_T1_ROI_GM_stats] = ttest2(Target_1_MRF_T1_ROI_GM_voxels, Target_2_MRF_T1_ROI_GM_voxels);
[MRF_T2_ROI_GM_h, MRF_T2_ROI_GM_p, MRF_T2_ROI_GM_ci, MRF_T2_ROI_GM_stats] = ttest2(Target_1_MRF_T2_ROI_GM_voxels, Target_2_MRF_T2_ROI_GM_voxels);
[MRF_T1_ROI_WM_h, MRF_T1_ROI_WM_p, MRF_T1_ROI_WM_ci, MRF_T1_ROI_WM_stats] = ttest2(Target_1_MRF_T1_ROI_WM_voxels, Target_2_MRF_T1_ROI_WM_voxels);
[MRF_T2_ROI_WM_h, MRF_T2_ROI_WM_p, MRF_T2_ROI_WM_ci, MRF_T2_ROI_WM_stats] = ttest2(Target_1_MRF_T2_ROI_WM_voxels, Target_2_MRF_T2_ROI_WM_voxels);

P_values = [MRF_T1_ROI_GM_p; MRF_T2_ROI_GM_p; MRF_T1_ROI_WM_p; MRF_T2_ROI_WM_p];
T_values = [MRF_T1_ROI_GM_stats.tstat; MRF_T2_ROI_GM_stats.tstat; MRF_T1_ROI_WM_stats.tstat; MRF_T2_ROI_WM_stats.tstat];

%% effect size
MRF_T1_ROI_GM_Effect = meanEffectSize(Target_1_MRF_T1_ROI_GM_voxels, Target_2_MRF_T1_ROI_GM_voxels);
MRF_T2_ROI_GM_Effect = meanEffectSize(Target_1_MRF_T2_ROI_GM_voxels, Target_2_MRF_T2_ROI_GM_voxels);
MRF_T1_ROI_WM_Effect = meanEffectSize(Target_1_MRF_T1_ROI_WM_voxels, Target_2_MRF_T1_ROI_WM_voxels);
MRF_T2_ROI_WM_Effect = meanEffectSize(Target_1_MRF_T2_ROI_WM_voxels, Target_2_MRF_T2_ROI_WM_voxels);

Effect_sizes = [MRF_T1_ROI_GM_Effect.Effect; MRF_T2_ROI_GM_Effect.Effect; MRF_T1_ROI_WM_Effect.Effect; MRF_T2_ROI_WM_Effect.Effect];

%% Spider plot
color_set = [0, 1; 0.447, 0; 0.741, 0]';
% color_set = [1, 1, 0.466, 1; 0, 0, 0.674, 0; 0, 0, 0.188, 0]';

Target_1_mean = [mean(Target_1_MRF_T1_ROI_GM_voxels), mean(Target_1_MRF_T2_ROI_GM_voxels), mean(Target_1_MRF_T1_ROI_WM_voxels), mean(Target_1_MRF_T2_ROI_WM_voxels)];
Target_2_mean = [mean(Target_2_MRF_T1_ROI_GM_voxels), mean(Target_2_MRF_T2_ROI_GM_voxels), mean(Target_2_MRF_T1_ROI_WM_voxels), mean(Target_2_MRF_T2_ROI_WM_voxels)];

Targets_mean = [Target_1_mean; Target_2_mean];

% Target_1_std = [std(Target_1_MRF_T1_ROI_GM_voxels), std(Target_1_MRF_T2_ROI_GM_voxels), std(Target_1_MRF_T1_ROI_WM_voxels), std(Target_1_MRF_T2_ROI_WM_voxels)];
% Target_2_std = [std(Target_2_MRF_T1_ROI_GM_voxels), std(Target_2_MRF_T2_ROI_GM_voxels), std(Target_2_MRF_T1_ROI_WM_voxels), std(Target_2_MRF_T2_ROI_WM_voxels)];
% 
% Targets_std = [ROI_Lt_std; ROI_Rt_std];
% 
% Target_1_stderr = [stderr_calc(Target_1_MRF_T1_ROI_GM_voxels), stderr_calc(Target_1_MRF_T2_ROI_GM_voxels), stderr_calc(Target_1_MRF_T1_ROI_WM_voxels), stderr_calc(Target_1_MRF_T2_ROI_WM_voxels)];
% Target_2_stderr = [stderr_calc(Target_2_MRF_T1_ROI_GM_voxels), stderr_calc(Target_2_MRF_T2_ROI_GM_voxels), stderr_calc(Target_2_MRF_T1_ROI_WM_voxels), stderr_calc(Target_2_MRF_T2_ROI_WM_voxels)];
% 
% Targets_stderr = [ROI_Lt_stderr; ROI_Rt_stderr];

AxesLimits = [-0.5, -0.5, -0.5, -0.5; 2, 2, 2, 2];

s = spider_plot(Targets_mean, 'AxesLabels', {'T1_GM', 'T2_GM', 'T1_WM', 'T2_WM'}, 'Color', color_set, 'AxesPrecision', [2, 2, 2, 2] , 'AxesLimits', AxesLimits);
set(gcf, 'Position', [680 276 782 702]);

legend('Target_1', 'Target_2', 'Location', 'southoutside','Interpreter','none');
saveas(gcf, fullfile(FCD_target_folder, 'Spider_plot.png')); 

%% violin plot
% T1 GM
grp1 = Target_1_MRF_T1_ROI_GM_voxels;
grp2 = Target_2_MRF_T1_ROI_GM_voxels;
grp1_class = repmat({'Target_1'}, 1, length(Target_1_MRF_T1_ROI_GM_voxels));
grp2_class = repmat({'Target_2'}, 1, length(Target_2_MRF_T1_ROI_GM_voxels));

grouporder={'Target_1', 'Target_2'};
figure, violinplot2([grp1; grp2], [grp1_class'; grp2_class'], 'GroupOrder', grouporder);
title('MRF_T1_GM','Interpreter','none')
saveas(gcf, fullfile(FCD_target_folder, 'Violiin_plot_T1_GM.png')); 

% T1 WM
grp1 = Target_1_MRF_T1_ROI_WM_voxels;
grp2 = Target_2_MRF_T1_ROI_WM_voxels;
grp1_class = repmat({'Target_1'}, 1, length(Target_1_MRF_T1_ROI_WM_voxels));
grp2_class = repmat({'Target_2'}, 1, length(Target_2_MRF_T1_ROI_WM_voxels));

grouporder={'Target_1', 'Target_2'};
figure, violinplot2([grp1; grp2], [grp1_class'; grp2_class'], 'GroupOrder', grouporder);
title('MRF_T1_WM','Interpreter','none')
saveas(gcf, fullfile(FCD_target_folder, 'Violiin_plot_T1_WM.png')); 

% T2 GM
grp1 = Target_1_MRF_T2_ROI_GM_voxels;
grp2 = Target_2_MRF_T2_ROI_GM_voxels;
grp1_class = repmat({'Target_1'}, 1, length(Target_1_MRF_T2_ROI_GM_voxels));
grp2_class = repmat({'Target_2'}, 1, length(Target_2_MRF_T2_ROI_GM_voxels));

grouporder={'Target_1', 'Target_2'};
figure, violinplot2([grp1; grp2], [grp1_class'; grp2_class'], 'GroupOrder', grouporder);
title('MRF_T2_GM','Interpreter','none')
saveas(gcf, fullfile(FCD_target_folder, 'Violiin_plot_T2_GM.png')); 

% T2 WM
grp1 = Target_1_MRF_T2_ROI_WM_voxels;
grp2 = Target_2_MRF_T2_ROI_WM_voxels;
grp1_class = repmat({'Target_1'}, 1, length(Target_1_MRF_T2_ROI_WM_voxels));
grp2_class = repmat({'Target_2'}, 1, length(Target_2_MRF_T2_ROI_WM_voxels));

grouporder={'Target_1', 'Target_2'};
figure, violinplot2([grp1; grp2], [grp1_class'; grp2_class'], 'GroupOrder', grouporder);
title('MRF_T2_WM','Interpreter','none')
saveas(gcf, fullfile(FCD_target_folder, 'Violiin_plot_T2_WM.png')); 
