#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N submit_all_combine

cd $PBS_O_WORKDIR

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-env

qiime feature-table merge \
  --i-tables sequences_in_work_ap4/all_BacV3_denoised_table_ap5.qza  --i-tables sequences_in_work_ap4/all_BacV5_denoised_table_ap5.qza \
  --o-merged-table sequences_in_work_ap4/all_bac_denoised_table_ap6.qza

qiime feature-table merge-seqs \
  --i-data sequences_in_work_ap4/all_BacV3_rep_seqs_ap5.qza --i-data sequences_in_work_ap4/all_BacV5_rep_seqs_ap5.qza \
  --o-merged-data sequences_in_work_ap4/all_bac_rep_seqs_ap6.qza



qiime feature-table merge \
  --i-tables sequences_in_work_ap4/all_FITS2_denoised_table_ap5.qza  --i-tables sequences_in_work_ap4/all_Ftrad_denoised_table_ap5.qza \
  --o-merged-table sequences_in_work_ap4/all_fun_denoised_table_ap6.qza

qiime feature-table merge-seqs \
  --i-data sequences_in_work_ap4/all_FITS2_rep_seqs_ap5.qza --i-data sequences_in_work_ap4/all_Ftrad_rep_seqs_ap5.qza \
  --o-merged-data sequences_in_work_ap4/all_fun_rep_seqs_ap6.qza



qiime feature-table merge \
  --i-tables sequences_in_work_ap4/all_OITS2_denoised_table_ap5.qza  --i-tables sequences_in_work_ap4/all_Otrad_denoised_table_ap5.qza \
  --o-merged-table sequences_in_work_ap4/all_oom_denoised_table_ap6.qza

qiime feature-table merge-seqs \
  --i-data sequences_in_work_ap4/all_OITS2_rep_seqs_ap5.qza --i-data sequences_in_work_ap4/all_Otrad_rep_seqs_ap5.qza \
  --o-merged-data sequences_in_work_ap4/all_oom_rep_seqs_ap6.qza



qiime feature-table merge \
  --i-tables sequences_in_work_ap4/all_PrV4_denoised_table_ap5.qza  --i-tables sequences_in_work_ap4/all_PrV9_denoised_table_ap5.qza \
  --o-merged-table sequences_in_work_ap4/all_euk_denoised_table_ap6.qza

qiime feature-table merge-seqs \
  --i-data sequences_in_work_ap4/all_PrV4_rep_seqs_ap5.qza --i-data sequences_in_work_ap4/all_PrV9_rep_seqs_ap5.qza \
  --o-merged-data sequences_in_work_ap4/all_euk_rep_seqs_ap6.qza