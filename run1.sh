#!/bin/bash
#SBATCH --job-name=harm_run_1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --mem=192G
#SBATCH -o harm_run_1.o.log
#SBATCH -e harm_run_1.e.log
source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh
conda activate splicing
rmats_input="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/data/rmats/SH_BIO_2207180_3x_vs_SH_DMSO_*.txt"
leafcutter_input="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/data/leafcutter/SH_BIO_2207180_3x_vs_SH_DMSO*_cluster_significance.txt"
majiq_input="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/data/majiq/SH_BIO_2207180_3x_vs_SH_DMSO*voila.tsv"
rmats_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/junction_prep/rmats_junction_prep.csv"
leafcutter_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/junction_prep/leafcutter_junction_prep.csv"
majiq_output="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/junction_prep/majiq_junction_prep.csv"
Rscript /home/ychen12/splicing_test/splicing_harmonization/DIRS/rmats_junction_prep_V3.R -i "$rmats_input" -o "$rmats_output" -n "$comparisonRef_name"
Rscript /home/ychen12/splicing_test/splicing_harmonization/DIRS/leafcutter_junction_prep_V2.R -i "$leafcutter_input" -o "$leafcutter_output"
Rscript /home/ychen12/splicing_test/splicing_harmonization/DIRS/majiq_junction_prep.R -i "$majiq_input" -o "$majiq_output"
annotation="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/data/stringtie/SH_BIO_2207180_3x_vs_SH_DMSO.combined.gtf"
reference="/edgehpc/dept/compbio/projects/splice_pipeline_BSSI/mry_annot/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf"
junction_input="/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/SH_BIO_2207180_3x_vs_SH_DMSO/junction_prep/{method}_junction_prep.csv"
Rscript /home/ychen12/splicing_test/splicing_harmonization/DIRS/annotation_SJ.R -a "$annotation" -r "$reference" -i "$junction_input"
