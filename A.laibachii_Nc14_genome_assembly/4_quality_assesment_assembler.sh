#PBS -l nodes=1:ppn=10
#PBS -l walltime=01:00:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xcoverage_table

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

mkdir -p quality_assesment_assembly
mkdir -p Logs

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate bioawk
echo "Running bioawk..."

STRAIN1=Hifiasm_Nc14_hap1
REF1=hifiasm_assembly/Nc14_asm.bp.hap1.p_ctg.fa
TYPE1=contig

LEN1=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF1` # Size of assembled genome

echo "line,length,type,coverage" > quality_assesment_assembly/length.csv

cat $REF1 | bioawk -c fastx -v line="$STRAIN1" '{print line","length($seq)","length($seq)}' | \
sort -k3rV -t "," | awk -F "," -v len="$LEN1" -v type="$TYPE1" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> quality_assesment_assembly/length.csv

STRAIN2=Hifiasm_Nc14_hap2
REF2=hifiasm_assembly/Nc14_asm.bp.hap2.p_ctg.fa
TYPE2=contig

LEN2=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF2` # Size of assembled genome

cat $REF2 | bioawk -c fastx -v line="$STRAIN2" '{print line","length($seq)","length($seq)}' | \
sort -k3rV -t "," | awk -F "," -v len="$LEN2" -v type="$TYPE2" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> quality_assesment_assembly/length.csv

STRAIN3=Canu_Nc14
REF3=canu_assembly/Nc14_asm.contigs.fasta
TYPE3=contig

LEN3=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF3` # Size of assembled genome

cat $REF3 | bioawk -c fastx -v line="$STRAIN3" '{print line","length($seq)","length($seq)}' | \
sort -k3rV -t "," | awk -F "," -v len="$LEN3" -v type="$TYPE3" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> quality_assesment_assembly/length.csv

conda deactivate
echo "Length script is finished"

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate busco38
echo "Running busco..."

busco -i $REF1 -o quality_assesment_assembly/$STRAIN1 -m genome --offline -c 10 -l stramenopiles_odb10
busco -i $REF2 -o quality_assesment_assembly/$STRAIN2 -m genome --offline -c 10 -l stramenopiles_odb10
busco -i $REF3 -o quality_assesment_assembly/$STRAIN3 -m genome --offline -c 10 -l stramenopiles_odb10

conda deactivate
echo "BUSCO is finished"

echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > quality_assesment_assembly/busco.csv

PREFIX1=quality_assesment_assembly/$STRAIN1
cat $PREFIX1/short_summary.specific.stramenopiles_odb10.$STRAIN1.txt | grep "(S)" | awk -v strain="$PREFIX1" '{print strain","$1}' > quality_assesment_assembly/complete_single.txt
cat $PREFIX1/short_summary.specific.stramenopiles_odb10.$STRAIN1.txt | grep "(D)" | awk '{print $1}' > quality_assesment_assembly/complete_duplicated.txt
cat $PREFIX1/short_summary.specific.stramenopiles_odb10.$STRAIN1.txt | grep "(F)" | awk '{print $1}' > quality_assesment_assembly/fragmented.txt
cat $PREFIX1/short_summary.specific.stramenopiles_odb10.$STRAIN1.txt | grep "(M)" | awk '{print $1}' > quality_assesment_assembly/missing.txt
paste -d "," quality_assesment_assembly/complete_single.txt quality_assesment_assembly/complete_duplicated.txt quality_assesment_assembly/fragmented.txt quality_assesment_assembly/missing.txt >> quality_assesment_assembly/busco.csv


PREFIX2=quality_assesment_assembly/$STRAIN2
cat $PREFIX2/short_summary.specific.stramenopiles_odb10.$STRAIN2.txt | grep "(S)" | awk -v strain="$PREFIX2" '{print strain","$1}' > quality_assesment_assembly/complete_single.txt
cat $PREFIX2/short_summary.specific.stramenopiles_odb10.$STRAIN2.txt | grep "(D)" | awk '{print $1}' > quality_assesment_assembly/complete_duplicated.txt
cat $PREFIX2/short_summary.specific.stramenopiles_odb10.$STRAIN2.txt | grep "(F)" | awk '{print $1}' > quality_assesment_assembly/fragmented.txt
cat $PREFIX2/short_summary.specific.stramenopiles_odb10.$STRAIN2.txt | grep "(M)" | awk '{print $1}' > quality_assesment_assembly/missing.txt
paste -d "," quality_assesment_assembly/complete_single.txt quality_assesment_assembly/complete_duplicated.txt quality_assesment_assembly/fragmented.txt quality_assesment_assembly/missing.txt >> quality_assesment_assembly/busco.csv

PREFIX3=quality_assesment_assembly/$STRAIN3
cat $PREFIX3/short_summary.specific.stramenopiles_odb10.$STRAIN3.txt | grep "(S)" | awk -v strain="$PREFIX3" '{print strain","$1}' > quality_assesment_assembly/complete_single.txt
cat $PREFIX3/short_summary.specific.stramenopiles_odb10.$STRAIN3.txt | grep "(D)" | awk '{print $1}' > quality_assesment_assembly/complete_duplicated.txt
cat $PREFIX3/short_summary.specific.stramenopiles_odb10.$STRAIN3.txt | grep "(F)" | awk '{print $1}' > quality_assesment_assembly/fragmented.txt
cat $PREFIX3/short_summary.specific.stramenopiles_odb10.$STRAIN3.txt | grep "(M)" | awk '{print $1}' > quality_assesment_assembly/missing.txt
paste -d "," quality_assesment_assembly/complete_single.txt quality_assesment_assembly/complete_duplicated.txt quality_assesment_assembly/fragmented.txt quality_assesment_assembly/missing.txt >> quality_assesment_assembly/busco.csv

# Delete temporary files
rm quality_assesment_assembly/complete_single.txt quality_assesment_assembly/complete_duplicated.txt quality_assesment_assembly/fragmented.txt quality_assesment_assembly/missing.txt

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate quast
echo "Quast environment activated"

quast.py -t 10 -o quality_assesment_assembly/quast_hifiasm1 -r ref/Nc14_A.laibachii.dna.toplevel.fa $REF1
quast.py -t 10 -o quality_assesment_assembly/quast_hifiasm2 -r ref/Nc14_A.laibachii.dna.toplevel.fa $REF2
quast.py -t 10 -o quality_assesment_assembly/quast_canu -r ref/Nc14_A.laibachii.dna.toplevel.fa $REF3

conda deactivate
echo "QUAST is finished"

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate r_env_stats
echo "R environment is running..."

Rscript 4_visualizing_asssembler_quality.R > Logs/4_visualizing_asssembler_quality.R.log 2> Logs/4_visualizing_asssembler_quality.R.err
Rscript 4_visualizing_busco_result.R > Logs/4_visualizing_busco_result.R.log 2> Logs/4_visualizing_busco_result.R.err

conda deactivate