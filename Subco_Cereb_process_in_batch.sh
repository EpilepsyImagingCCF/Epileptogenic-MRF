folder=$1

for sub in ${folder}/P*
do

filename=`basename ${sub}/n_syN_T1w_data_brain_Warped.nii`

echo "GM/WM/CSF seg"
fast -g ${sub}/${filename}.nii

echo "subco seg"
run_first_all -b -i ${sub}/${filename}.nii -o ${sub}/${filename}_SubcoSeg.nii

echo "first_flirt"
first_flirt ${sub}/${filename}.nii ${sub}/${filename}_to_std_sub -cort -b

echo "run_first"
run_first -i ${sub}/${filename}.nii -t ${sub}/${filename}_to_std_sub_cort.mat -n 320 -o ${sub}/${filename}_L_Cereb -m ${FSLDIR}/data/first/models_336_bin/intref_puta/L_Cereb.bmv -intref ${FSLDIR}/data/first/models_336_bin/05mm/L_Puta_05mm.bmv
run_first -i ${sub}/${filename}.nii -t ${sub}/${filename}_to_std_sub_cort.mat -n 320 -o ${sub}/${filename}_R_Cereb -m ${FSLDIR}/data/first/models_336_bin/intref_puta/R_Cereb.bmv -intref ${FSLDIR}/data/first/models_336_bin/05mm/R_Puta_05mm.bmv

echo "boundary corr"
first_boundary_corr -s ${sub}/${filename}_L_Cereb -i ${sub}/${filename}.nii -b none -o ${sub}/${filename}_L_Cereb_corr
first_boundary_corr -s ${sub}/${filename}_R_Cereb -i ${sub}/${filename}.nii -b none -o ${sub}/${filename}_R_Cereb_corr

echo "combination"
fslmaths ${sub}/${filename}_SubcoSeg.nii -add ${sub}/${filename}_R_Cereb_corr.nii.gz -add ${sub}/${filename}_L_Cereb_corr.nii.gz -thr 0.5 -bin ${sub}/${filename}_Subco_Cereb_mask.nii.gz

done
