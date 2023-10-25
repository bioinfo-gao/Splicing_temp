#Directories, argparse output directory This LEAFcutter only works in ================================Yarrow ==========<<<<<<<<<<<<<
prog_dir="/camhpc/dept/compbio/project/splice_pipeline_BSSI/programs/dev/"
anno_dir="/camhpc/dept/compbio/project/splice_pipeline_BSSI/data/annotation/Human.GRCh38.v34.l1_5.ERCC/"
bam_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/bam/"

# changes the names for all comparisons
##### iPSC (A-I)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-12nM/"   # (1) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-40nM/"   # (2)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-10uM/"   # (3)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-3uM/"    # (4)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-150nM/"  # (5)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-45nM/"   # (6)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-261nM/"  # (7)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-870nM/"  # (8)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-10uM/"   # (9)
out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/"    # (10) #(10) -3uM/ <<<<<== leafcutter=failed #
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-1000nM/" # (11)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-3700nM/" # (12)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-1680nM/" # (13)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-5600nM/" # (14)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-1uM/"    # (15)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-300nM/"  # (16) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-1000nM/" # (17)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-300nM/"  # (18)

#argparse select bam files
##### iPSC BAM (A-II)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1755497-?-12nM*.bam)   #(1)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1755497-?-40nM*.bam)   #(2)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1949634-?-10uM*.bam)   #(3)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1949634-?-3uM*.bam)    #(4)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2006152-?-150nM*.bam)  #(5)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2006152-?-45nM*.bam)   #(6)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2059811-?-261nM*.bam)  #(7)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2059811-?-870nM*.bam)  #(8)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060573-?-10uM*.bam)   #(9)
alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060573-?-3uM*.bam)    #(10)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060884-?-1000nM*.bam) #(11)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060884-?-3700nM*.bam) #(12)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2070692-?-1680nM*.bam) #(13)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2070692-?-5600nM*.bam) #(14)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2135644-?-1uM*.bam)    #(15)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2135644-?-300nM*.bam)  #(16)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2136770-?-1000nM*.bam) #(17)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2136770-?-300nM*.bam)  #(18)

# REF (A-III)
ref_bam_file=$(ls ${bam_dir}'TST11872_iPSC-DMSO-'*.bam)         # TST11742_HNDS-0030-DMSO-=TST11742_GM24468D-DMSO-

###########################================================================================================================================================================================
# ##### NGN2 (B-I)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1755497-12nM/"   # (1)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1755497-40nM/"   # (2) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1949634-10uM/"   # (3) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-1949634-3uM/"    # (4) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2006152-150nM/"  # (5)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2006152-45nM/"   # (6) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2059811-261nM/"  # (7) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2059811-870nM/"  # (8) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060573-10uM/"   # (9) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060573-3uM/"    # (10) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-1000nM/" # (11) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-3700nM/" # (12) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2070692-1680nM/" # (13) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2070692-5600nM/" # (14) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2135644-1uM/"    # (15) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2135644-300nM/"  # (16) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2136770-1000nM/" # (17) 
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2136770-300nM/"  # (18) 

# BAM (B-II)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-?-12nM*.bam)   #(1) #alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-[1-3]-12nM*.bam)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-?-40nM*.bam)   #(2)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1949634-?-10uM*.bam)   #(3)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1949634-?-3uM*.bam)    #(4)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2006152-?-150nM*.bam)  #(5)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2006152-?-45nM*.bam)   #(6)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2059811-?-261nM*.bam)  #(7)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2059811-?-870nM*.bam)  #(8)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060573-?-10uM*.bam)   #(9)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060573-?-3uM*.bam)    #(10)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060884-?-1000nM*.bam) #(11)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060884-?-3700nM*.bam) #(12)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2070692-?-1680nM*.bam) #(13)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2070692-?-5600nM*.bam) #(14)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2135644-?-1uM*.bam)    #(15)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2135644-?-300nM*.bam)  #(16)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2136770-?-1000nM*.bam) #(17)
#alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2136770-?-300nM*.bam)  #(18)

##### ref B-III
#ref_bam_file=$(ls ${bam_dir}'TST11872_NGN2-DMSO-'*.bam)         # TST11742_HNDS-0030-DMSO-=TST11742_GM24468D-DMSO-


##### ========================= log_dir SAME
log_dir="$out_dir/leafcut_logs/"; mkdir -p $log_dir

#argparse remaining options
read_type="paired"
read_length=100

anno_file="$anno_dir/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf"
gff3_file="$anno_dir/Human.GRCh38.v34.l1_5.ERCC.transcript.rtracklayer.gff3"
leafcutter_strandedness=1
leafcutter_anno_file="$anno_dir/Human.GRCh38.v34.l1_5.ERCC.transcript.txt.gz"
#No covariates available

#Prefer to make a multiline string to see the full command
multiline="source /etc/profile.d/modules_bash.sh; module load python/3.5.1;
python $prog_dir/biogen_splicing_pipeline_mry_alltools_v3.py
--leafcutter
--alt_bam_file $alt_bam_file
--ref_bam_file $ref_bam_file
--leafcutter_strandedness $leafcutter_strandedness
--leafcutter_anno_file $leafcutter_anno_file
--out_dir $out_dir
"

#Submit job
#may need additional module sourcing  #  ginseng02, ginseng03, ,.. up to ginseng18
echo $multiline | qsub -N out.leafcutter -o $log_dir -q cpu.q@ginseng02 -l h_rt=272:00:00 -l h_vmem=192G # Job rejected for unknown queue "cpu.q".  qhost -q by ZG
#echo $multiline | qsub -N out.leafcutter -o $log_dir -q cpu.q@ginseng02 -l h_rt=272:00:00 -l h_vmem=192G # Job rejected for unknown queue "cpu.q".  qhost -q by ZG

