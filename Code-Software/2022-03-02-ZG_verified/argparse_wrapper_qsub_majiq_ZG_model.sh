#Directories, argparse output directory
prog_dir="/camhpc/dept/compbio/project/splice_pipeline_BSSI/programs/dev/"
anno_dir="/camhpc/dept/compbio/project/splice_pipeline_BSSI/data/annotation/Human.GRCh38.v34.l1_5.ERCC/"
bam_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/bam/"

# changes the names for all comparisons
# ##### iPSC 
# 1)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-1755497-12nM/" # ; mkdir -p $out_dir
# 2)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-1755497-40nM/" # ; mkdir -p $out_dir
# 3)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-1949634-10uM/" # ; mkdir -p $out_dir
# 4)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-1949634-3uM/" # ; mkdir -p $out_dir
# 5)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2006152-150nM/" # ; mkdir -p $out_dir
# 6)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2006152-45nM/" # ; mkdir -p $out_dir
# 7)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2059811-261nM/" # ; mkdir -p $out_dir
# 8)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2059811-870nM/" # ; mkdir -p $out_dir
# 9)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2060573-10uM/" # ; mkdir -p $out_dir
# 10)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2060573-3uM/" # ; mkdir -p $out_dir
# 11)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2060884-1000nM/" # ; mkdir -p $out_dir
# 12)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2060884-3700nM/" # ; mkdir -p $out_dir
# 13)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2070692-1680nM/" # ; mkdir -p $out_dir
# 14)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2070692-5600nM/" # ; mkdir -p $out_dir
# 15)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2135644-1uM/" # ; mkdir -p $out_dir
# 16)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2135644-300nM/" # ; mkdir -p $out_dir
# 17)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2136770-1000nM/" # ; mkdir -p $out_dir
# 18)
out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_iPSC-BIO-2136770-300nM/" # ; mkdir -p $out_dir
# ##### NGN2
# 1)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-1755497-12nM/" # ; mkdir -p $out_dir
# 2)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-1755497-40nM/" # ; mkdir -p $out_dir
# 3)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-1949634-10uM/" # ; mkdir -p $out_dir
# 4)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-1949634-3uM/" # ; mkdir -p $out_dir
# 5)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2006152-150nM/" # ; mkdir -p $out_dir
# 6)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2006152-45nM/" # ; mkdir -p $out_dir
# 7)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2059811-261nM/" # ; mkdir -p $out_dir
# 8)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2059811-870nM/" # ; mkdir -p $out_dir
# 9)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2060573-10uM/" # ; mkdir -p $out_dir
# 10)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2060573-3uM/" # ; mkdir -p $out_dir
# 11)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2060884-1000nM/" # ; mkdir -p $out_dir
# 12)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2060884-3700nM/" # ; mkdir -p $out_dir
# 13)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2070692-1680nM/" # ; mkdir -p $out_dir
# 14)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2070692-5600nM/" # ; mkdir -p $out_dir
# 15)
# out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2135644-1uM/" # ; mkdir -p $out_dir
# 16)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2135644-300nM/" # ; mkdir -p $out_dir
# 17)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2136770-1000nM/" # ; mkdir -p $out_dir
# 18)
#out_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/majiq/TST11872_NGN2-BIO-2136770-300nM/" # ; mkdir -p $out_dir
##### log_dir
log_dir="$out_dir/logs/"; mkdir -p $log_dir


#argparse select bam files
##### iPSC 
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1755497-?-12nM*.bam)  #(1)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1755497-?-40nM*.bam) #(2)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1949634-?-10uM*.bam) #(3)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-1949634-?-3uM*.bam) #(4)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2006152-?-150nM*.bam) #(5)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2006152-?-45nM*.bam) #(6)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2059811-?-261nM*.bam) #(7)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2059811-?-870nM*.bam) #(8)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060573-?-10uM*.bam) #(9)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060573-?-3uM*.bam) #(10)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060884-?-1000nM*.bam) #(11)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2060884-?-3700nM*.bam) #(12)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2070692-?-1680nM*.bam) #(13)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2070692-?-5600nM*.bam) #(14)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2135644-?-1uM*.bam) # #(15)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2135644-?-300nM*.bam) #(16)
#alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2136770-?-1000nM*.bam) #(17)
alt_bam_file=$(ls ${bam_dir}TST11872_iPSC-BIO-2136770-?-300nM*.bam) #(18)
# 
ref_bam_file=$(ls ${bam_dir}'TST11872_iPSC-DMSO-'*.bam)         # TST11742_HNDS-0030-DMSO-=TST11742_GM24468D-DMSO-
# 
# 
# ##### NGN2 
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-?-12nM*.bam) ##(1) alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-[1-3]-12nM*.bam)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1755497-?-40nM*.bam) #(2)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1949634-?-10uM*.bam) #(3)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-1949634-?-3uM*.bam) ##(4)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2006152-?-150nM*.bam) #(5)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2006152-?-45nM*.bam) # #(6)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2059811-?-261nM*.bam) #(7)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2059811-?-870nM*.bam) #(8)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060573-?-10uM*.bam) #(9)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060573-?-3uM*.bam) ##(10)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060884-?-1000nM*.bam)#(11)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2060884-?-3700nM*.bam)#(12)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2070692-?-1680nM*.bam) #(13)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2070692-?-5600nM*.bam) #(14)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2135644-?-1uM*.bam) # #(15)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2135644-?-300nM*.bam) #(16)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2136770-?-1000nM*.bam)#(17)
# alt_bam_file=$(ls ${bam_dir}TST11872_NGN2-BIO-2136770-?-300nM*.bam) #(18)
# 
# ref_bam_file=$(ls ${bam_dir}'TST11872_NGN2-DMSO-'*.bam)         # TST11742_HNDS-0030-DMSO-=TST11742_GM24468D-DMSO-





#argparse remaining options
read_type="paired"
read_length=100
genome_file="$anno_dir/hg19as_SMN1_masked.fa"
majiq_cutoffs="$(seq 0.1 0.1 0.4)"

anno_file="$anno_dir/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf"
gff3_file="$anno_dir/Human.GRCh38.v34.l1_5.ERCC.transcript.rtracklayer.gff3"

#Prefer to make a multiline string to see the full command
multiline="source /etc/profile.d/modules_bash.sh; module load python/3.5.1;
python $prog_dir/biogen_splicing_pipeline_mry_alltools_v3.py
--majiq_v1
--alt_bam_file $alt_bam_file
--ref_bam_file $ref_bam_file
--read_type $read_type
--read_length $read_length
--gff3_file $gff3_file
--genome_file $genome_file
--majiq_cutoffs $majiq_cutoffs
--out_dir $out_dir
"

#Submit job
#may need additional module sourcing  
echo $multiline | qsub -o $log_dir -N out.majiq -q all.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 12

