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

Target_3_MRF_T1_ROI_GM_voxels;
Target_3_MRF_T2_ROI_GM_voxels;
Target_3_MRF_T1_ROI_WM_voxels;
Target_3_MRF_T2_ROI_WM_voxels;

...

%%
% T1 GM
ft_T1_groups = {Target_1_MRF_T1_ROI_GM_voxels, Target_2_MRF_T1_ROI_GM_voxels, Target_3_MRF_T1_ROI_GM_voxels, Target_4_MRF_T1_ROI_GM_voxels};
fts_T1_groups_name = {...
    repmat({'Target_1_T1_GM'}, length(Target_1_MRF_T1_ROI_GM_voxels), 1), ...
    repmat({'Target_1_T1_GM'}, length(Target_2_MRF_T1_ROI_GM_voxels), 1), ...
    repmat({'Target_1_T1_GM'}, length(Target_3_MRF_T1_ROI_GM_voxels), 1), ...
    repmat({'Target_1_T1_GM'}, length(Target_4_MRF_T1_ROI_GM_voxels), 1)};

%norm test
for order = 1:length(ft_T1_groups)
    [~, p] = lillietest(ft_T1_groups{order});
    P_sets(order, 1) = p;
    clear p
end

% Multiple-sample tests
for order = 1:length(ft_T1_groups)
    if order == 1
        fts_T1_all = ft_T1_groups{order};
        fts_T1_grp = fts_T1_groups_name{order};
    else
        fts_T1_all = [fts_T1_all; ft_T1_groups{order}];
        fts_T1_grp = [fts_T1_grp; fts_T1_groups_name{order}];
    end
end

[T1_vartest_p, T1_vartest_stats] = vartestn(fts_T1_all, fts_T1_grp);

[T1_anova1_p, T1_anova1_tbl, T1_anova1_stats] = anova1(fts_T1_all, fts_T1_grp);
[T1_mc_c, T1_mc_m, T1_mc_h, T1_mc_gnames] = multcompare(T1_anova1_stats); % 'CriticalValueType', 'bonferroni'
set(gca,'TickLabelInterpreter','none')

nComparisons = size(T1_mc_c, 1);
cohens_d = zeros(nComparisons, 1);
Grp = {'Target_1_T1_GM', 'Target_1_T1_GM', 'Target_1_T1_GM', 'Target_2_T1_GM', 'Target_2_T1_GM', 'Target_3_T1_GM'};
Grp_control = {'Target_2_T1_GM', 'Target_3_T1_GM', 'Target_4_T1_GM', 'Target_3_T1_GM', 'Target_4_T1_GM', 'Target_4_T1_GM'};

for i = 1:nComparisons
    group1_idx = find(strcmp(Grp{i}, fts_T1_grp));
    group2_idx = find(strcmp(Grp_control{i}, fts_T1_grp));

    group1_data = fts_T1_all(group1_idx);
    group2_data = fts_T1_all(group2_idx);

    mean_diff = mean(group1_data) - mean(group2_data);
    pooled_std = sqrt((var(group1_data) * (length(group1_data) - 1) + var(group2_data) * (length(group2_data) - 1)) / (length(group1_data) + length(group2_data) - 2));

    cohens_d(i, 1) = mean_diff / pooled_std;

    clear group1_idx group2_idx group1_data group2_data mean_diff pooled_std
end

T1_mc_tbl = array2table(T1_mc_c, "VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
T1_mc_tbl.("Group") = T1_mc_gnames(T1_mc_tbl.("Group"));
T1_mc_tbl.("Control Group") = T1_mc_gnames(T1_mc_tbl.("Control Group"));

% T1 WM
...
 
% T2 GM
...

% T2 WM
...

%% Spider plot
% color_set = [1, 0, 0.466, 0.494; 0, 0.447, 0.674, 0.184; 0, 0.741, 0.188, 0.556]';
color_set = [1, 1, 0.466, 1; 0, 0, 0.674, 0; 0, 0, 0.188, 0]';

Target_1_mean = [mean(Target_1_MRF_T1_ROI_GM_voxels), mean(Target_1_MRF_T2_ROI_GM_voxels), mean(Target_1_MRF_T1_ROI_WM_voxels), mean(Target_1_MRF_T2_ROI_WM_voxels)];
Target_2_mean = [mean(Target_2_MRF_T1_ROI_GM_voxels), mean(Target_2_MRF_T2_ROI_GM_voxels), mean(Target_2_MRF_T1_ROI_WM_voxels), mean(Target_2_MRF_T2_ROI_WM_voxels)];
Target_3_mean = [mean(Target_3_MRF_T1_ROI_GM_voxels), mean(Target_3_MRF_T2_ROI_GM_voxels), mean(Target_3_MRF_T1_ROI_WM_voxels), mean(Target_3_MRF_T2_ROI_WM_voxels)];
Target_4_mean = [mean(Target_4_MRF_T1_ROI_GM_voxels), mean(Target_4_MRF_T2_ROI_GM_voxels), mean(Target_4_MRF_T1_ROI_WM_voxels), mean(Target_4_MRF_T2_ROI_WM_voxels)];

Targets_mean = [Target_1_mean; Target_2_mean; Target_3_mean; Target_4_mean];

% Target_1_std = [std(Target_1_MRF_T1_ROI_GM_voxels), std(Target_1_MRF_T2_ROI_GM_voxels), std(Target_1_MRF_T1_ROI_WM_voxels), std(Target_1_MRF_T2_ROI_WM_voxels)];
% Target_2_std = [std(Target_2_MRF_T1_ROI_GM_voxels), std(Target_2_MRF_T2_ROI_GM_voxels), std(Target_2_MRF_T1_ROI_WM_voxels), std(Target_2_MRF_T2_ROI_WM_voxels)];
% Target_3_std = [std(Target_3_MRF_T1_ROI_GM_voxels), std(Target_3_MRF_T2_ROI_GM_voxels), std(Target_3_MRF_T1_ROI_WM_voxels), std(Target_3_MRF_T2_ROI_WM_voxels)];
% Target_4_std = [std(Target_4_MRF_T1_ROI_GM_voxels), std(Target_4_MRF_T2_ROI_GM_voxels), std(Target_4_MRF_T1_ROI_WM_voxels), std(Target_4_MRF_T2_ROI_WM_voxels)];
% 
% Targets_std = [Target_1_std; Target_2_std; Target_3_std; Target_4_std];
% 
% Target_1_stderr = [stderr_calc(Target_1_MRF_T1_ROI_GM_voxels), stderr_calc(Target_1_MRF_T2_ROI_GM_voxels), stderr_calc(Target_1_MRF_T1_ROI_WM_voxels), stderr_calc(Target_1_MRF_T2_ROI_WM_voxels)];
% Target_2_stderr = [stderr_calc(Target_2_MRF_T1_ROI_GM_voxels), stderr_calc(Target_2_MRF_T2_ROI_GM_voxels), stderr_calc(Target_2_MRF_T1_ROI_WM_voxels), stderr_calc(Target_2_MRF_T2_ROI_WM_voxels)];
% Target_3_stderr = [stderr_calc(Target_3_MRF_T1_ROI_GM_voxels), stderr_calc(Target_3_MRF_T2_ROI_GM_voxels), stderr_calc(Target_3_MRF_T1_ROI_WM_voxels), stderr_calc(Target_3_MRF_T2_ROI_WM_voxels)];
% Target_4_stderr = [stderr_calc(Target_4_MRF_T1_ROI_GM_voxels), stderr_calc(Target_4_MRF_T2_ROI_GM_voxels), stderr_calc(Target_4_MRF_T1_ROI_WM_voxels), stderr_calc(Target_4_MRF_T2_ROI_WM_voxels)];
% 
% Targets_stderr = [Target_1_stderr; Target_2_stderr; Target_3_stderr; Target_4_stderr];

% AxesLimits = [floor(min(P(:, 1))*1000)/1000-0.001, floor(min(P(:, 2))*1000)/1000-0.001, floor(min(P(:, 3))*1000)/1000-0.001, floor(min(P(:, 4))*1000)/1000-0.001; ...
%     ceil(max(P(:, 1))*1000)/1000+0.001, ceil(max(P(:, 2))*1000)/1000+0.001, ceil(max(P(:, 3))*1000)/1000+0.001, ceil(max(P(:, 4))*1000)/1000+0.001];
% , 'AxesLimits', AxesLimits , 'AxesPrecision', [1, 1, 1, 1]

s = spider_plot(Targets_mean, 'AxesLabels', {'T1_GM', 'T2_GM', 'T1_WM', 'T2_WM'}, 'Color', color_set, 'AxesPrecision', [2, 2, 2, 2]);
set(gcf, 'Position', [680 276 782 702]);

legend('Target_1', 'Target_2', 'Target_3', 'Target_4', 'Location', 'southoutside','Interpreter','none');
saveas(gcf, fullfile(FCD_target_folder, 'Spider_plot.png')); 

%%
% T1 GM
grp1 = Target_1_MRF_T1_ROI_GM_voxels;
grp2 = Target_2_MRF_T1_ROI_GM_voxels;
grp3 = Target_3_MRF_T1_ROI_GM_voxels;
grp4 = Target_4_MRF_T1_ROI_GM_voxels;

grp1_class = repmat({'Target_1'}, 1, length(Target_1_MRF_T1_ROI_GM_voxels));
grp2_class = repmat({'Target_2'}, 1, length(Target_2_MRF_T1_ROI_GM_voxels));
grp3_class = repmat({'Target_3'}, 1, length(Target_3_MRF_T1_ROI_GM_voxels));
grp4_class = repmat({'Target_4'}, 1, length(Target_4_MRF_T1_ROI_GM_voxels));

grouporder={'Target_1', 'Target_2', 'Target_3', 'Target_4'};
figure, violinplot2([grp1; grp2; grp3; grp4], [grp1_class'; grp2_class'; grp3_class'; grp4_class'], 'GroupOrder', grouporder);
title('MRF_T1_GM','Interpreter','none')
saveas(gcf, fullfile(FCD_EJ_folder, 'Violiin_plot_T1_GM.png')); 

% T1 WM
...

% T2 GM
...

% T2 WM
...