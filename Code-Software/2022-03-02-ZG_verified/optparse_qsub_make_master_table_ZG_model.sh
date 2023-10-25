#optparse options ===========<<<<<<<< in Yarrow
prog_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/splicing/"

#======== ### Prefix ====== iPSC
#out_prefix="TST11872_iPSC-BIO-1755497-12nM"   #(1)
#out_prefix="TST11872_iPSC-BIO-1755497-40nM"   #(2)
#out_prefix="TST11872_iPSC-BIO-1949634-10uM"   #(3)
#out_prefix="TST11872_iPSC-BIO-1949634-3uM"    #(4)
#out_prefix="TST11872_iPSC-BIO-2006152-150nM"  #(5)
#out_prefix="TST11872_iPSC-BIO-2006152-45nM"   #(6)
#out_prefix="TST11872_iPSC-BIO-2059811-261nM"  #(7)
#out_prefix="TST11872_iPSC-BIO-2059811-870nM"  #(8)
#out_prefix="TST11872_iPSC-BIO-2060573-10uM"   #(9)
#out_prefix="TST11872_iPSC-BIO-2060573-3uM"    #(10) # analysis files are missing from /Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/ <<<<<<<<<<<<<<<<============ leafcutter =====failed #
#out_prefix="TST11872_iPSC-BIO-2060884-1000nM" #(11)
#out_prefix="TST11872_iPSC-BIO-2060884-3700nM" #(12) # analysis files are missing from /Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-3700nM/ <<<<<<<<<<<<<<<<============ majiq_v1 =====failed #
#out_prefix="TST11872_iPSC-BIO-2070692-1680nM" #(13)
#out_prefix="TST11872_iPSC-BIO-2070692-5600nM" #(14)
#out_prefix="TST11872_iPSC-BIO-2135644-1uM"    #(15)
#out_prefix="TST11872_iPSC-BIO-2135644-300nM"  #(16)
#out_prefix="TST11872_iPSC-BIO-2136770-1000nM" #(17)
#out_prefix="TST11872_iPSC-BIO-2136770-300nM"  #(18)
#======== ### Prefix ====== iPSC
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-12nM/"   #(1)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-40nM/"   #(2) 15:42 15:45
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-10uM/"   #(3)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-3uM/"    #(4)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-150nM/"  #(5)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-45nM/"   #(6)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-261nM/"  #(7)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-870nM/"  #(8)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-10uM/"   #(9)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/"    #(10)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-1000nM/" #(11)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-3700nM/" #(12)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-1680nM/" #(13)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-5600nM/" #(14)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-1uM/"    #(15)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-300nM/"  #(16)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-1000nM/" #(17)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-300nM/"  #(18)
# ######=================================##############
# ### Prefix ====== NGN2
#out_prefix="TST11872_NGN2-BIO-1755497-12nM"   #(1)
#out_prefix="TST11872_NGN2-BIO-1755497-40nM"   #(2)
#out_prefix="TST11872_NGN2-BIO-1949634-10uM"   #(3) Error: analysis files are missing from /Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1949634-10uM/ <<<<<<<<<<<<<<<============ majiq_v1 =====failed #
#out_prefix="TST11872_NGN2-BIO-1949634-3uM"    #(4)
#out_prefix="TST11872_NGN2-BIO-2006152-150nM"  #(5)
#out_prefix="TST11872_NGN2-BIO-2006152-45nM"   #(6)
out_prefix="TST11872_NGN2-BIO-2059811-261nM"  #(7) #/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2059811-261nM/       <<<<<<<<<<<<<<<============ majiq_v1 =====failed #
#out_prefix="TST11872_NGN2-BIO-2059811-870nM"  #(8)
#out_prefix="TST11872_NGN2-BIO-2060573-10uM"   #(9)
#out_prefix="TST11872_NGN2-BIO-2060573-3uM"    #(10)
#out_prefix="TST11872_NGN2-BIO-2060884-1000nM" #(11)
#out_prefix="TST11872_NGN2-BIO-2060884-3700nM" #(12)
#out_prefix="TST11872_NGN2-BIO-2070692-1680nM" #(13)
#out_prefix="TST11872_NGN2-BIO-2070692-5600nM" #(14)
#out_prefix="TST11872_NGN2-BIO-2135644-1uM"    #(15)
#out_prefix="TST11872_NGN2-BIO-2135644-300nM"  #(16) ??TST11872_NGN2-BIO-2135644-300nM TST11872_NGN2-BIO-2135644-1uM runned twice, but this missed. 
#out_prefix="TST11872_NGN2-BIO-2136770-1000nM" #(17)
#out_prefix="TST11872_NGN2-BIO-2136770-300nM"  #(18)
# ### Prefix ====== NGN2
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1755497-12nM/"   #(1)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1755497-40nM/"   #(2)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1949634-10uM/"   #(3)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1949634-3uM/"    #(4)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2006152-150nM/"  #(5)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2006152-45nM/"   #(6)
in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2059811-261nM/"  #(7)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2059811-870nM/"  #(8)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060573-10uM/"   #(9)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060573-3uM/"    #(10)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-1000nM/" #(11)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-3700nM/" #(12)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2070692-1680nM/" #(13)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2070692-5600nM/" #(14)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2135644-1uM/"    #(15)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2135644-300nM/"  #(16)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2136770-1000nM/" #(17)
#in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2136770-300nM/"  #(18)

out_dir=$in_dir
majiq_cutoffs="0.1,0.2,0.3,0.4"


#submit job
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q@ginseng18 -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V 
# cpu.q long.q sucessfully finished in Yarrow
# Unable to run job: Job was rejected because job requests unknown queue "cpu.q@ginseng18.hpc.biogen.com". Exiting. chang to ++Yarrow
