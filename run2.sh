#!/bin/bash
#SBATCH --job-name=harm_run_2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH -o harm_run_2.o.log
#SBATCH -e harm_run_2.e.log
source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh
conda activate splicing
harm_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm"
harm_intron_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm_prep/intron_list.csv"
harm_exonlabel_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm_prep/exon_label.csv"
harm_sj_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm_prep/sj_annotation.csv"
python /home/ychen12/splicing_test/splicing_harmonization/DIRS/rmats_harm_prep.py -exonlabel $harm_exonlabel_output -i $harm_intron_output -sj $harm_sj_output -o $harm_output
python /home/ychen12/splicing_test/splicing_harmonization/DIRS/leafcutter_harm_prep.py -exonlabel $harm_exonlabel_output -i $harm_intron_output -sj $harm_sj_output -o $harm_output
python /home/ychen12/splicing_test/splicing_harmonization/DIRS/majiq_harm_prep.py -exonlabel $harm_exonlabel_output -i $harm_intron_output -sj $harm_sj_output -o $harm_output
relabel="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm_prep/exon_label.csv"
file_indir="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm"
file_outdir="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm/out"
module load openmpi
mpirun -np $SLURM_NTASKS python /home/ychen12/splicing_test/splicing_harmonization/DIRS/harm_assign_V2.py -exonlabel $harm_exonlabel_output -indir $file_indir -outdir $file_outdir
module unload openmpi
