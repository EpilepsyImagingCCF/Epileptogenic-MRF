clear all; close all;

toolbox_path = '/path/ to/ toolbox';

addpath(fullfile(toolbox_path, 'NIfTI_20140122'))
addpath(fullfile(toolbox_path, 'Violin_plot'))
addpath(fullfile(toolbox_path, 'Spider_plot_github_repo'))
addpath(fullfile(toolbox_path, 'Silent_FCD'))

FCD_target_folder = '/path/ to/ target subjects';
HCs_folder = '/path/ to/ control subjects';

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
FCD_target_list = dir(fullfile(FCD_target_folder, 'P*'));
threshold = 500;

load(fullfile(HCs_folder, 'MRF_T1_imgset.mat'))
load(fullfile(HCs_folder, 'MRF_T2_imgset.mat'))
load(fullfile(HCs_folder, 'FAST_imgset.mat'))

count = 1;
for order = 1:length(FCD_target_list)
    if FCD_target_list(order).isdir == 1

        disp([num2str(order) '/'  num2str(length(FCD_target_list))]);

        T1_filename = dir(fullfile(FCD_target_list(order).folder, FCD_target_list(order).name, 'n_syN_T1_*_Warped.nii'));
        T2_filename = dir(fullfile(FCD_target_list(order).folder, FCD_target_list(order).name, 'n_syN_T2_*_Warped.nii'));
        MRF_T1_nii = load_untouch_nii(fullfile(T1_filename.folder, T1_filename.name));
        MRF_T2_nii = load_untouch_nii(fullfile(T2_filename.folder, T2_filename.name));
        MRF_T1_img = double(MRF_T1_nii.img);
        MRF_T2_img = double(MRF_T2_nii.img);
        
        % intensity truncation
        MRF_T2_img(find(MRF_T2_img > threshold)) = threshold;
        
        % lesion ROI
        ROI_nii = load_untouch_nii(fullfile(FCD_target_list(order).folder, FCD_target_list(order).name, 'n_ROI_bin.nii'));
        ROI_img = double(ROI_nii.img);
        ROI_WM_img = zeros(size(ROI_img));
        ROI_GM_img = zeros(size(ROI_img));
        
        Subco_Cereb_nii = load_untouch_nii(fullfile(FCD_target_list(order).folder, FCD_target_list(order).name, 'n_syN_T1w_data_brain_Warped_Subco_Cereb_mask.nii.gz'));
        Subco_Cereb_img = double(Subco_Cereb_nii.img);
        
        % CSF; remove CSF region from ROI, ROI_GM and ROI_WM
        FAST_nii = load_untouch_nii(fullfile(FCD_target_list(order).folder, FCD_target_list(order).name, 'n_syN_T1w_data_brain_Warped_pveseg.nii.gz'));
        FAST_img = double(FAST_nii.img);
        CSF_img = zeros(size(FAST_img));
        GM_img = zeros(size(FAST_img));
        WM_img = zeros(size(FAST_img));
        CSF_img(find(FAST_img == 1)) = 1;
        GM_img(find(FAST_img == 2)) = 1;
        WM_img(find(FAST_img == 3)) = 1; 
        
        if CSF_dil_flag == 1
            SE = strel("cube", 2);
            CSF_img = imdilate(CSF_img, SE);
        end
        
        ROI_GM_img(find(ROI_img + GM_img == 2)) = 1;
        ROI_GM_img(find(CSF_img)) = 0;
        ROI_WM_img(find(ROI_img + WM_img == 2)) = 1;
        ROI_WM_img(find(CSF_img)) = 0;

        if count == 1
            %
            [Norm_MRF_T1_ROI_GM_values, Norm_MRF_T2_ROI_GM_values, Norm_MRF_T1_ROI_WM_values, Norm_MRF_T2_ROI_WM_values] = Intersubject_normalization(MRF_T1_img, MRF_T2_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, MRF_T2_imgset, FAST_imgset, ori_prop_flag);
            Target_MRF_T1_ROI_GM_voxels = Norm_MRF_T1_ROI_GM_values;
            Target_MRF_T2_ROI_GM_voxels = Norm_MRF_T2_ROI_GM_values;
            
            Target_MRF_T1_ROI_WM_voxels = Norm_MRF_T1_ROI_WM_values;
            Target_MRF_T2_ROI_WM_voxels = Norm_MRF_T2_ROI_WM_values;
            clear Norm_MRF_T1_ROI_GM_values Norm_MRF_T2_ROI_GM_values Norm_MRF_T1_ROI_WM_values Norm_MRF_T2_ROI_WM_values
            
        else
            % 
            [Norm_MRF_T1_ROI_GM_values, Norm_MRF_T2_ROI_GM_values, Norm_MRF_T1_ROI_WM_values, Norm_MRF_T2_ROI_WM_values] = Intersubject_normalization(MRF_T1_img, MRF_T2_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, MRF_T2_imgset, FAST_imgset, ori_prop_flag);
            Target_MRF_T1_ROI_GM_voxels = cat(1, Target_MRF_T1_ROI_GM_voxels, Norm_MRF_T1_ROI_GM_values);
            Target_MRF_T2_ROI_GM_voxels = cat(1, Target_MRF_T2_ROI_GM_voxels, Norm_MRF_T2_ROI_GM_values);
            
            Target_MRF_T1_ROI_WM_voxels = cat(1, Target_MRF_T1_ROI_WM_voxels, Norm_MRF_T1_ROI_WM_values);
            Target_MRF_T2_ROI_WM_voxels = cat(1, Target_MRF_T2_ROI_WM_voxels, Norm_MRF_T2_ROI_WM_values);
            clear Norm_MRF_T1_ROI_GM_values Norm_MRF_T2_ROI_GM_values Norm_MRF_T1_ROI_WM_values Norm_MRF_T2_ROI_WM_values
            
        end
        
        count = count + 1;
        clear T1_filename T2_filename MRF_T1_nii MRF_T2_nii MRF_T1_img MRF_T2_img semiAct_ROI_nii semiAct_ROI_img semiAct_ROI_WM_img semiAct_ROI_GM_img ROI_nii ROI_img ROI_WM_img ROI_GM_img Subco_Cereb_nii Subco_Cereb_img FAST_nii FAST_img CSF_img GM_img WM_img Normal_ROI_GM_img Normal_ROI_WM_img SE
        clear ROI_WM_img_dil ROI_GM_img_dil semiAct_ROI_WM_img_dil semiAct_ROI_GM_img_dil Subco_Cereb_img_dil Lesion_ROI_dil
    end 
end

%
save(fullfile(FCD_target_folder, ['Target_MRF_T1_ROI_GM_voxels_' out_postfix '.mat']), 'Target_MRF_T1_ROI_GM_voxels')
save(fullfile(FCD_target_folder, ['Target_MRF_T2_ROI_GM_voxels_' out_postfix '.mat']), 'Target_MRF_T2_ROI_GM_voxels')
save(fullfile(FCD_target_folder, ['Target_MRF_T1_ROI_WM_voxels_' out_postfix '.mat']), 'Target_MRF_T1_ROI_WM_voxels')
save(fullfile(FCD_target_folder, ['Target_MRF_T2_ROI_WM_voxels_' out_postfix '.mat']), 'Target_MRF_T2_ROI_WM_voxels')
