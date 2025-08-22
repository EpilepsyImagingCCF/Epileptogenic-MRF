# Utility of MR Fingerprinting in Differentiating Epileptogenic from Non-Epileptogenic Cortical Malformations

This repository accompanies the research paper:

Ryuzaburo Kochi, Ting-Yu Su, Vineet Punia, Spencer Morris, Hiroatsu Murakami, Xiaofeng Wang, Ingmar Blümcke, Stephen E. Jones, Imad Najm, Andreas V. Alexopoulos, Dan Ma, Zhong Irene Wang,  
*Utility of MR fingerprinting in differentiating epileptogenic from non-epileptogenic cortical malformations*,  
Journal of the Neurological Sciences, Volume 477, 2025, 123651,  
ISSN 0022-510X,  
[https://doi.org/10.1016/j.jns.2025.123651](https://doi.org/10.1016/j.jns.2025.123651)

---

## Pipeline Overview
This pipeline outlines the preprocessing and statistical analysis steps used in our study.

### 1. Data Organization
Organize your data in the following structure:
```
Data/
├── FCD_I/
│ ├── Patient01.nii
│ ├── Patient02.nii
│ └── ...
├── Target_subject/
│ └── Patient01.nii
├── HCs/
│ ├── V01.nii
│ ├── V02.nii
│ └── ...
├── ...
```

### 2. Image Registration
- Run `regis_2_MNI.sh` to perform ANTs registration.  
- This registers the images and ROIs into MNI space.

### 3. Tissue Segmentation
Run `Subco_Cereb_process_in_batch.sh` to:  
- Segment GM, WM, and CSF regions.  
- Segment the subcortical regions (`run_first_all`).  
- Generate the left and right cerebellum masks.  
- Combine the subcortical regions and cerebellum masks together.  

### 4. Normalization Data Generation
- Run `HCs_data_gen_4_norm.m` to generate the data for normalization.  
- ⚠️ Note: If there are multiple ROIs, please revise the code accordingly.  

### 5. Group Normalization
- Run `Target_group_norm.m` to normalize and collect data from a specific group.  
- ⚠️ Note: If there are multiple ROIs, please revise the code accordingly.  

### 6. Two-Group Statistical Analysis
- Run `Two_grps_stats_ana.m` to perform a **two-sample statistical analysis** between two targets.  
- Targets can be two different groups or two different ROIs.  

### 7. Multiple-Group Statistical Analysis
- Run `Multiple_grps_stats_ana.m` to perform a **multiple-group statistical analysis** between multiple targets.  
- Targets can be multiple groups or multiple ROIs.  

### 8. Visualization
- If needed, run `Violin_plot_python.ipynb` in Python to plot the violin plots.  

---
