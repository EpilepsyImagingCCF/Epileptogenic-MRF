function [Norm_MRF_T1_ROI_GM_values, Norm_MRF_T2_ROI_GM_values, Norm_MRF_T1_ROI_WM_values, Norm_MRF_T2_ROI_WM_values] = Intersubject_normalization(MRF_T1_img, MRF_T2_img, ROI_GM_img, ROI_WM_img, MRF_T1_imgset, MRF_T2_imgset, FAST_imgset, ori_prop_flag)
% ori_prop_flag: 1 for non-normalized processing; 0 for normalized processing
% 1. MRF_T1, MRF_T2 find the voxels in ROI(ROI, ROI_GM, ROI_WM)
% 2. divided by average from HCs mean map
if ori_prop_flag == 1
    HCs_MRF_T1_ROI_GM_mean = 0;
    HCs_MRF_T2_ROI_GM_mean = 0;
    HCs_MRF_T1_ROI_WM_mean = 0;
    HCs_MRF_T2_ROI_WM_mean = 0;

    HCs_MRF_T1_ROI_GM_std = 1;
    HCs_MRF_T1_ROI_WM_std = 1;
    HCs_MRF_T2_ROI_GM_std = 1;
    HCs_MRF_T2_ROI_WM_std = 1;
else
    HCs_MRF_T1_ROI_GM_values = [];
    HCs_MRF_T2_ROI_GM_values = [];
    HCs_MRF_T1_ROI_WM_values = [];
    HCs_MRF_T2_ROI_WM_values = [];

    for HC_order = 1:size(MRF_T1_imgset, 4)
        MRF_T1_imgset_per = MRF_T1_imgset(:, :, :, HC_order);
        MRF_T2_imgset_per = MRF_T2_imgset(:, :, :, HC_order);
        FAST_imgset_per = FAST_imgset(:, :, :, HC_order);

        FAST_GM_imgset_per = zeros(size(FAST_imgset_per));
        FAST_WM_imgset_per = zeros(size(FAST_imgset_per));
        FAST_CSF_imgset_per = zeros(size(FAST_imgset_per));

        FAST_CSF_imgset_per(find(FAST_imgset_per == 1)) = 1;
        FAST_GM_imgset_per(find(FAST_imgset_per == 2)) = 1;
        FAST_WM_imgset_per(find(FAST_imgset_per == 3)) = 1;

        SE = strel("cube", 2);
        FAST_CSF_imgset_per = imdilate(FAST_CSF_imgset_per, SE);

        FAST_GM_imgset_per(find(FAST_CSF_imgset_per)) = 0;
        FAST_WM_imgset_per(find(FAST_CSF_imgset_per)) = 0;

        HCs_MRF_T1_ROI_GM_values_per = MRF_T1_imgset_per(find(ROI_GM_img + FAST_GM_imgset_per == 2));
        HCs_MRF_T1_ROI_GM_values_per = HCs_MRF_T1_ROI_GM_values_per(find(HCs_MRF_T1_ROI_GM_values_per));
        HCs_MRF_T1_ROI_GM_values = cat(1, HCs_MRF_T1_ROI_GM_values, HCs_MRF_T1_ROI_GM_values_per);
        clear HCs_MRF_T1_ROI_GM_values_per 

        HCs_MRF_T2_ROI_GM_values_per = MRF_T2_imgset_per(find(ROI_GM_img + FAST_GM_imgset_per == 2));
        HCs_MRF_T2_ROI_GM_values_per = HCs_MRF_T2_ROI_GM_values_per(find(HCs_MRF_T2_ROI_GM_values_per));
        HCs_MRF_T2_ROI_GM_values = cat(1, HCs_MRF_T2_ROI_GM_values, HCs_MRF_T2_ROI_GM_values_per);
        clear HCs_MRF_T2_ROI_GM_values_per 

        HCs_MRF_T1_ROI_WM_values_per = MRF_T1_imgset_per(find(ROI_WM_img + FAST_WM_imgset_per == 2));
        HCs_MRF_T1_ROI_WM_values_per = HCs_MRF_T1_ROI_WM_values_per(find(HCs_MRF_T1_ROI_WM_values_per));
        HCs_MRF_T1_ROI_WM_values = cat(1, HCs_MRF_T1_ROI_WM_values, HCs_MRF_T1_ROI_WM_values_per);
        clear HCs_MRF_T1_ROI_WM_values_per 

        HCs_MRF_T2_ROI_WM_values_per = MRF_T2_imgset_per(find(ROI_WM_img + FAST_WM_imgset_per == 2));
        HCs_MRF_T2_ROI_WM_values_per = HCs_MRF_T2_ROI_WM_values_per(find(HCs_MRF_T2_ROI_WM_values_per));
        HCs_MRF_T2_ROI_WM_values = cat(1, HCs_MRF_T2_ROI_WM_values, HCs_MRF_T2_ROI_WM_values_per);
        clear HCs_MRF_T2_ROI_WM_values_per 

        clear MRF_T1_imgset_per MRF_T2_imgset_per 
        clear FAST_imgset_per FAST_GM_imgset_per FAST_WM_imgset_per FAST_CSF_imgset_per SE
    end

    HCs_MRF_T1_ROI_GM_mean = mean(HCs_MRF_T1_ROI_GM_values(find(HCs_MRF_T1_ROI_GM_values)));
    HCs_MRF_T2_ROI_GM_mean = mean(HCs_MRF_T2_ROI_GM_values(find(HCs_MRF_T2_ROI_GM_values)));
    HCs_MRF_T1_ROI_WM_mean = mean(HCs_MRF_T1_ROI_WM_values(find(HCs_MRF_T1_ROI_WM_values)));
    HCs_MRF_T2_ROI_WM_mean = mean(HCs_MRF_T2_ROI_WM_values(find(HCs_MRF_T2_ROI_WM_values)));

    HCs_MRF_T1_ROI_GM_std = std(HCs_MRF_T1_ROI_GM_values(find(HCs_MRF_T1_ROI_GM_values)));
    HCs_MRF_T2_ROI_GM_std = std(HCs_MRF_T2_ROI_GM_values(find(HCs_MRF_T2_ROI_GM_values)));
    HCs_MRF_T1_ROI_WM_std = std(HCs_MRF_T1_ROI_WM_values(find(HCs_MRF_T1_ROI_WM_values)));
    HCs_MRF_T2_ROI_WM_std = std(HCs_MRF_T2_ROI_WM_values(find(HCs_MRF_T2_ROI_WM_values)));
 
end

MRF_T1_ROI_GM_values = MRF_T1_img(find(ROI_GM_img)); MRF_T1_ROI_GM_values = MRF_T1_ROI_GM_values(find(MRF_T1_ROI_GM_values));
MRF_T2_ROI_GM_values = MRF_T2_img(find(ROI_GM_img)); MRF_T2_ROI_GM_values = MRF_T2_ROI_GM_values(find(MRF_T2_ROI_GM_values));

MRF_T1_ROI_WM_values = MRF_T1_img(find(ROI_WM_img)); MRF_T1_ROI_WM_values = MRF_T1_ROI_WM_values(find(MRF_T1_ROI_WM_values));
MRF_T2_ROI_WM_values = MRF_T2_img(find(ROI_WM_img)); MRF_T2_ROI_WM_values = MRF_T2_ROI_WM_values(find(MRF_T2_ROI_WM_values));

Norm_MRF_T1_ROI_GM_values = (MRF_T1_ROI_GM_values - HCs_MRF_T1_ROI_GM_mean)/HCs_MRF_T1_ROI_GM_std;
Norm_MRF_T2_ROI_GM_values = (MRF_T2_ROI_GM_values - HCs_MRF_T2_ROI_GM_mean)/HCs_MRF_T2_ROI_GM_std;

Norm_MRF_T1_ROI_WM_values = (MRF_T1_ROI_WM_values - HCs_MRF_T1_ROI_WM_mean)/HCs_MRF_T1_ROI_WM_std;
Norm_MRF_T2_ROI_WM_values = (MRF_T2_ROI_WM_values - HCs_MRF_T2_ROI_WM_mean)/HCs_MRF_T2_ROI_WM_std;