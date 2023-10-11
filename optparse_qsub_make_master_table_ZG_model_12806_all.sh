#optparse options ===========<<<<<<<< in Yarrow
#prog_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/splicing/" #prog_dir="/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/splicing/"
#prog_dir="/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/BACK-contain_R4_bad_samples/Splicing-Detection-08-08" #prog_dir="/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/Splicing-Detection-08-08"
prog_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/2023_All_data_code_RNA3.06/" # Fatal error: cannot open file '/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code/TST12086/2023_modified_code_used/sub.make_master_table_v2.R': No such file or directory

majiq_cutoffs="0.1,0.2,0.3,0.4"

#======== ### Prefix ====== iPSC     1111111111111
out_prefix="iPSC_1755497_10x" # (1)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_1755497_10x/" # (1)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_1755497_3x"  # (2) NO "majiq_v1" Errorin "_reference_alternative.deltapsi.tsv" 12h
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_1755497_3x/"  # (2)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2006152_10x" # (3)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2006152_10x/" # (3)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2006152_3x"  # (4)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2006152_3x/"  # (4) # for leaf
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2189972_10x" # (5)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2189972_10x/" # (5)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2189972_3x"  # (6)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2189972_3x/"  # (6)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2194943_10x" # (7)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2194943_10x/" # (7)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2194943_3x"  # (8)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2194943_3x/"  # (8)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2195127_10x" # (9)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2195127_10x/" # (9)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2195127_3x"  # (10)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2195127_3x/"  # (10)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2195327_10x" # (11)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2195327_10x/" # (11)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2195327_3x"  # (12)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2195327_3x/"  # (12)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2197294_10x" # (13)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2197294_10x/" # (13)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2197294_3x"  # (14)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2197294_3x/"  # (14)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_1755497_10x_bridge"  # (15)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_1755497_10x_bridge/"  # (15)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_1755497_3x_bridge"   # (16)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_1755497_3x_bridge/"   # (16)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2006152_10x_bridge"  # (17)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2006152_10x_bridge/"  # (17)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_2006152_3x_bridge"   # (18)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2006152_3x_bridge/"   # (18)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="iPSC_DMSO_bridge"         # (19)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_DMSO_bridge/"         # (19)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

#======== ### Prefix ====== NGN2
out_prefix="NGN2_1755497_10x" # (1)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_1755497_10x/" # (1)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_1755497_3x"  # (2) NO "majiq_v1" Errorin "_reference_alternative.deltapsi.tsv" 12h
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_1755497_3x/"  # (2)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2006152_10x" # (3)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2006152_10x/" # (3)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2006152_3x"  # (4)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2006152_3x/"  # (4) # for leaf
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2189972_10x" # (5)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2189972_10x/" # (5)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2189972_3x"  # (6)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2189972_3x/"  # (6)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2194943_10x" # (7)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2194943_10x/" # (7)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2194943_3x"  # (8)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2194943_3x/"  # (8)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2195127_10x" # (9)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2195127_10x/" # (9)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2195127_3x"  # (10)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2195127_3x/"  # (10)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2195327_10x" # (11)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2195327_10x/" # (11)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2195327_3x"  # (12)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2195327_3x/"  # (12)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2197294_10x" # (13)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2197294_10x/" # (13)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="NGN2_2197294_3x"  # (14)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2197294_3x/"  # (14)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

#======== ### Prefix ====== Sy5Y
out_prefix="Sy5Y_1755497_10x" # (1)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_1755497_10x/" # (1)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_1755497_3x"  # (2) NO "majiq_v1" Errorin "_reference_alternative.deltapsi.tsv" 12h
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_1755497_3x/"  # (2)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2006152_10x" # (3)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2006152_10x/" # (3)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2006152_3x"  # (4)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2006152_3x/"  # (4) # for leaf
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2189972_10x" # (5)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2189972_10x/" # (5)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2189972_3x"  # (6)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2189972_3x/"  # (6)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2194943_10x" # (7)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2194943_10x/" # (7)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2194943_3x"  # (8)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2194943_3x/"  # (8)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2195127_10x" # (9)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2195127_10x/" # (9)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2195127_3x"  # (10)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2195127_3x/"  # (10)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2195327_10x" # (11)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2195327_10x/" # (11)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2195327_3x"  # (12)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2195327_3x/"  # (12)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2197294_10x" # (13)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2197294_10x/" # (13)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2197294_3x"  # (14)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2197294_3x/"  # (14)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_1755497_10x_bridge"  # (15)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_1755497_10x_bridge/"  # (15)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_1755497_3x_bridge"   # (16)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_1755497_3x_bridge/"   # (16)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2006152_10x_bridge"  # (17)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2006152_10x_bridge/"  # (17)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_2006152_3x_bridge"   # (18)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2006152_3x_bridge/"   # (18)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V

out_prefix="Sy5Y_DMSO_bridge"         # (19)
in_dir="/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_DMSO_bridge/"         # (19)
out_dir=$in_dir
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V 


# cpu.q long.q sucessfully finished in Yarrow
# Unable to run job: Job was rejected because job requests unknown queue 'cpu.q@ginseng18.hpc.biogen.com'. Exiting. chang to ++Yarrow @ginseng\d* regex delet all @ginseng01 -18
