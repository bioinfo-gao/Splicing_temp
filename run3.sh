#!/bin/bash
#SBATCH --job-name=harm_run_3
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH -o harm_run_3.o.log
#SBATCH -e harm_run_3.e.log
source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh
conda activate splicing
file_outdir="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/harm/out"
final_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/out"
final_allgene_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/out/all_gene"
python /home/ychen12/splicing_test/splicing_harmonization/DIRS/post_process.py -indir $file_outdir -outdir $final_output -allgene $final_allgene_output