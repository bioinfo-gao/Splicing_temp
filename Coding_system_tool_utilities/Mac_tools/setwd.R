setwd()
# [1] "/camhpc/home/zgao1"
setwd("/Volumes/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/")
setwd("/Volumes/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/motif_code")
setwd("/Users/zgao1/Documents/Biogen_Project/Python_scripts/")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-300nM/majiq_v1/")
setwd("/camhpc/home/dhuh/project_RNAseq/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs/")
setwd("/home/jpoulo/SPLICE/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs") # joy splicing folder
setwd("/home/zgao1/NGS_projects/JoySplicing") # my soft link
#For the codes that Dan ran, you can refer to 
setwd("/home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing/")
#      /home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/run.00.download_from_DNAnexus_mehools_code.qsub

dir.create("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/path")



"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|-*_master_table.csv" # this is where results from all three algorithms are combined.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- leafcutter"   #  leafcutter results. Data from “leafcutter_ds_res_cluster_significance.txt” and “leafcutter_ds_res_effect_sizes.txt” are used for the master table.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- majiq_v1"   # majiq result. Data from “cutoff_*_reference_alternative.deltapsi.tsv” are used for the master table.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- rmats"        # rmats result. Data from *.MATS.JCEC.txt” are used for the master table.


Splice analysis example codes:
  For (1), please refer to 
“/camhpc/dept/compbio/project/splice_pipeline_BSSI/programs/dev/example_codes/” # , but I haven’t checked it myself.

For the codes that I ran, you can refer to 
“/home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing/”,

From “run.01.*” to  “run.06.*”

https://platform.dnanexus.com/projects/FPXy9Pj0G3b983Zv4x939bbz/monitor/analysis/G7yfpKQ0G3b9XFVY2KXY9Xfy
https://platform.dnanexus.com/projects/FPXy9Pj0G3b983Zv4x939bbz/data/analyses/TST11872


in the Edge folder

/edgehpc/dept/combio/zgao1 ??


# (base_ZG) /edgehpc/dept/compbio_old/projects z
# 
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11821
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11822
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11823
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11824
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11825
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11826
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11827
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11828
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11829
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11831
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11832
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11835
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11838
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11839
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11843
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11844
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11846
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11847
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11850
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11854
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11856
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11858
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11862
# dr-xr-xr-x  1 zouyang   compbio                     0 Mar 14  2022 TST11866
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11867
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST11879
# dr-xr-xr-x  1 zouyang   compbio                     0 Sep 26 18:07 TST11953
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TST5144
# dr-xr-sr-x  1 root      ngs                         0 Dec  3  2021 TSTUHRR
# dr-x--S---  1 root      ngs                         0 Dec  3  2021 TXP-of-RUUO-mouse-kidney
# dr-x--S---  1 root      ngs                         0 Dec  3  2021 UHRR
# dr-xr-x---  1 583172302 drdl_ukbiobank_imaging      0 Sep 14 11:05 UKBB_fundus
# -r-xr-x---  1 583172302              583172302  38071 Feb 24  2021 UKBB_fundus_du.out.txt
# dr-x--S---  1 root      ngs                         0 Dec  3  2021 UMRR
# dr-x--S---  1 root      ngs                         0 Dec  3  2021 Zoukhri-IL1-injury-model