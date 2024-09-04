#PBS -l nodes=1:ppn=12
#PBS -l walltime=47:00:00
#PBS -l mem=120gb
#PBS -S /bin/bash
#PBS -N xfastspar_endo

cd $PBS_O_WORKDIR

source miniconda3/etc/profile.d/conda.sh
conda activate fastspar

INPUTFOLDER=$PWD/Input
OTU_TABLE=wl_otu_fp10.tsv

mkdir bootstrap_counts bootstrap_correlation

fastspar --otu_table $OTU_TABLE  --correlation median_correlation.tsv --covariance median_covariance.tsv > log1

fastspar_bootstrap --otu_table $OTU_TABLE --number 1000 --prefix bootstrap_counts/ > log2

parallel -j 12 fastspar --otu_table {} --correlation bootstrap_correlation/cor{/} --covariance bootstrap_correlation/cov{/} -i 50 ::: bootstrap_counts/*  > log3

fastspar_pvalues --otu_table $OTU_TABLE --correlation median_correlation.tsv --prefix bootstrap_correlation/cor --permutations 1000 --outfile pvalues.tsv > log4
