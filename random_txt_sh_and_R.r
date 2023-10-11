module load anaconda3
/edgehpc/dept/compbio/users/zgao1/Bash_and_sub_code zgao1@edge-hpc-log-101 $ 

$  conda env list
# conda environments:
#
base                  *  /edgehpc/apps/gb/anaconda3/4.9.2
agat                     /edgehpc/apps/gb/anaconda3/4.9.2/envs/agat
awscli                   /edgehpc/apps/gb/anaconda3/4.9.2/envs/awscli
clone_for_testing        /edgehpc/apps/gb/anaconda3/4.9.2/envs/clone_for_testing
compchem                 /edgehpc/apps/gb/anaconda3/4.9.2/envs/compchem
compchem_backup          /edgehpc/apps/gb/anaconda3/4.9.2/envs/compchem_backup
epi2melabs-wf-cas9       /edgehpc/apps/gb/anaconda3/4.9.2/envs/epi2melabs-wf-cas9
gpu_playground           /edgehpc/apps/gb/anaconda3/4.9.2/envs/gpu_playground
jptr_pure                /edgehpc/apps/gb/anaconda3/4.9.2/envs/jptr_pure
jptrclone                /edgehpc/apps/gb/anaconda3/4.9.2/envs/jptrclone
juphubbackup             /edgehpc/apps/gb/anaconda3/4.9.2/envs/juphubbackup
jupyter_dcv              /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyter_dcv
jupyterhub               /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub
jupyterhub-dev           /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub-dev
jupyterhub-dev.old       /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub-dev.old
jupyterhub.bkp           /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub.bkp
jupyterhub.bkp.20340324     /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub.bkp.20340324
jupyterhub_old           /edgehpc/apps/gb/anaconda3/4.9.2/envs/jupyterhub_old
leafcutter               /edgehpc/apps/gb/anaconda3/4.9.2/envs/leafcutter
majiq-1.1.4              /edgehpc/apps/gb/anaconda3/4.9.2/envs/majiq-1.1.4
multiqc                  /edgehpc/apps/gb/anaconda3/4.9.2/envs/multiqc
r-kernel-3.6             /edgehpc/apps/gb/anaconda3/4.9.2/envs/r-kernel-3.6
r-kernel-4.1             /edgehpc/apps/gb/anaconda3/4.9.2/envs/r-kernel-4.1
seminar                  /edgehpc/apps/gb/anaconda3/4.9.2/envs/seminar
toolsetv1                /edgehpc/apps/gb/anaconda3/4.9.2/envs/toolsetv1
base_ZG                  /home/zgao1/.conda/envs/base_ZG
base_ZG_py27             /home/zgao1/.conda/envs/base_ZG_py27
                         /home/zgao1/.local/share/r-miniconda
                         /home/zgao1/.local/share/r-miniconda/envs/r-reticulate

(base) /edgehpc/dept/compbio/users/zgao1/Bash_and_sub_code zgao1@edge-hpc-log-101 $ 
$  conda activate base_ZG
(base_ZG) /edgehpc/dept/compbio/users/zgao1/Bash_and_sub_code zgao1@edge-hpc-log-101 $ 
$  tree 
.
├── run.01_1.make_optparse_mastertable.bash
template=./template_script/optparse_qsub_make_master_table.sh
	
currdir=$(pwd)
for din in $currdir/TST*; do
	doutname=${din##*/}
	cp $template $din
	sed -i "s/INDIR/$doutname/g" $doutname/optparse_qsub_make_master_table.sh
done


├── run.01.make_dir_bash_per_treament.bash

template=./template_script/*.bash

while IFS="" read -r line; do
	lin=($line)
       	alt=${lin[0]} # alt prefix
       	ref=${lin[1]} # ref prefix
	
	temp1=${alt//X-7/X7}
	temp2=${temp1//C-1/C1}
	doutname=${temp2/S-0/S0} #output directory names

	mkdir -p $doutname # make the output directory per treatment
	
	cp $template $doutname # copy the template scripts
	sed -i "s/OUTDIR/$doutname/g" $doutname/*.bash
	sed -i "s/ALTBAMFILE/$alt/g" $doutname/*.bash # modify the scripts per treatment
	sed -i "s/REFBAMFILE/$ref/g" $doutname/*.bash

done < ./config/config.txt

├── run.02.make_dir_bash_per_treament_exclude_rep4.bash
for din in TST*; do
	dout=${din}_excrep4
	mkdir $dout

	cp $din/*.bash $dout

	sed -i 's/.bam)/.bam | grep -v rep4 )/g' $dout/*.bash
	sed -i "s/$din/$dout/g" $dout/*.bash
done

├── run.03.run_rmats_camhpc.bash
currdir=$(pwd)
for din in $currdir/TST*excrep4; do
	echo $din
        bash $din/argparse_wrapper_qsub_rmats.bash
done

├── run.04.run_leafcutter_yarrow.bash
currdir=$(pwd)

for din in $currdir/TST*; do
	if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist, 
       		 bash $din/argparse_wrapper_qsub_leafcutter.bash # then run the leafcutter job
	fi
done

├── run.05.run_majiq_yarrow.bash
currdir=$(pwd)

for din in $currdir/TST*; do
	if [ ! -d ${din}/majiq_v1 ]; then
        	bash $din/argparse_wrapper_qsub_majiq.bash
	fi		
done
├── run.06.make_mastertable.bash
currdir=$(pwd)

for din in $currdir/TST*; do
	bash $din/optparse_qsub_make_master_table.sh
done

for din in TST*; do
	if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist,
	  sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_leafcutter.bash
        fi
done
├── run.z.run.z.makeover_leafcutter.bash
#
#for din in TST*excrep4; do
#        if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist,
#          sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_leafcutter.bash
#		echo 
#        fi
#done

├── run.z2.makeover_majiq.bash
for din in TST*; do
	if [ ! -d ${din}/majiq_v1 ]; then # if leafcutter directory does not exist,
	  #sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_majiq.bash
	  sed -i "s/cpu.q@ginseng03/cpu.q@ginseng04/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/cpu.q@ginseng01/cpu.q/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/-V cpu.q@ginseng20/cpu.q/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/thread 12/thread 12 -V/g"  $din/argparse_wrapper_qsub_majiq.bash
        fi
done

├── run.z3.sortout_majiq.bash
for din in TST*; do
	if compgen -G "$din/majiq_v1/cutoff_0.1*" > /dev/null; then
#            ls $din/majiq_v1/cutoff_0.1*
		echo "ex"
        else 
            rm -r $din/majiq_v1/
#		echo $din	
        fi
done

├── run.z4.sortout_leafcuttuer_with_covar.bash
for din in TST*; do
	if compgen -G "$din/majiq_v1/cutoff_0.1*" > /dev/null; then
#            ls $din/majiq_v1/cutoff_0.1*
		echo "ex"
        else 
            rm -r $din/majiq_v1/
#		echo $din	
        fi
done

├── run.z.makeover_leafcutter.bash

└── template_script
    ├── argparse_wrapper_qsub_leafcutter.bash
    ├── argparse_wrapper_qsub_leafcutter_covar.bash
    ├── argparse_wrapper_qsub_majiq.bash
    ├── argparse_wrapper_qsub_rmats.bash
    └── optparse_qsub_make_master_table.sh
### 13 ============######### failed :
[1]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/pkKfK9Q67yX90kjpgPgb5B6KZK51v2BPpyQKyq9B/NGN2-DMSO-3.Aligned.sortedByCoord.out.bam
[7]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/fKkG93xZ6YZPpv6KgbqpjPxF2PYjvq7vY17j17YZ/Sy5Y-BIO-1755497-3x-2.Aligned.sortedByCoord.out.bam
[9]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/3V0x1qXj63gQp43Py1j96j25FBgy783J1z8g12qY/Sy5Y-BIO-1755497-3x-3.Aligned.sortedByCoord.out.bam
[11]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/gkYbzv5YpFVVGqvKbGpz45xVb85qZbKxkxV26yyq/Sy5Y-BIO-1755497-3x-4.Aligned.sortedByCoord.out.bam
[13]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/ky7YK19K7jPQb00Kp4Z7VfyBbQ6X888kKY6F8ZG7/Sy5Y-BIO-1755497-3x-bridge-b.Aligned.sortedByCoord.out.bam
[17]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/x47b32q1V18K62Pk899GKKbpzbG0ZVg0YP62ZJYJ/Sy5Y-BIO-1755497-10x-1.Aligned.sortedByCoord.out.bam
[19]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/94x1BgZVxfqpG3bXygq22gZvy3Ppj0fQxjqq3VjV/Sy5Y-BIO-1755497-10x-2.Aligned.sortedByCoord.out.bam
[21]   Exit 18                 nohup curl -O https://dl.dnanex.us/F/D/16gzy5P07vQJ2V57JPbb14Bv9yJ3jjf5vp9ZkpB6/Sy5Y-BIO-1755497-10x-3.Aligned.sortedByCoord.out.bam


 /camhpc/home/zgao1/NGS_projects/Generally_Procedure/camhpc/home/zgao1/NGS_projects/Generally_Procedure
1) obtain the TST number , here TST 11955
1.1A) directly click the link in 
http://ngs.biogen.com/ngs_one/app//core/app_task_chat.php?ID=2628#Message_972

1.1B) Goto
http://ngs.biogen.com/ngs_one/  search TST11872, read the basic information


2) Goto https://fastr.biogen.com/analysis
2.1) seleect the TST11872

2.2) in the samples 
click "select all"


2.3) Review Reference
Species
homo_sapiens

2.4) Version
Human.GRCh38.v34.l1_5.ERCC


2.5) clcik the bottom "start run analysis" 

a bubule came out in the topright.

the process takes 18 hours.

# Here’s the leafviz example in “https://rstudio-prod.hpc.biogen.com/”
# Please be sure to use “R4.2.0”
########################
rm(list=ls())
# ‘rstan’, ‘Hmisc’, ‘DirichletMultinomial’, ‘TailRank’, ‘StanHeaders’, ‘RcppEigen’ are not available for package ‘leafcutter’
# https://stackoverflow.com/questions/52902174/how-to-install-hmisc-in-rhttps://stackoverflow.com/questions/52902174/how-to-install-hmisc-in-r
# install.packages("survival")
# install.packages("lattice")
# install.packages("ggplot2")
library("survival")
library("lattice")
library("ggplot2")
# devtools::install_version ("RcppEigen", "0.3.3.7.0") # 4install.packages("RcppEigen", dependencies = T) # F
# install.packages("interp", dependencies = T)
# install.packages("latticeExtra", dependencies = T)
# install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz") #install.packages("foreign", dependencies = T)
# install.packages("Hmisc")
# install.packages('rpart') # https://stackoverflow.com/questions/7937282/rpart-package-installation-in-r
# library("rpart")
library("Hmisc")

# install.packages(c("latticeExtra", "copula", "ellipse", "gridBase", "locfit", "logspline", "mapproj", "maps", "MEMSS", "mlmRev", "RColorBrewer" ))
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("flowCore", "flowViz", "hexbin"))
# install.packages("devtools") # if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org') # devtools::install_github("stan-dev/rstantools") # not working
# install.packages("remotes")
# remotes::install_github("RcppCore/RcppParallel")
# devtools::install_github("stan-dev/rstantools")
# install.packages("rstan", dependencies = T)
# install.packages("shinystan", dependencies = T)
# install.packages("shiny")
# install.packages("ggrepel")
# library("leafcutter") # use R4.2.0
# if (!require("BiocManager", quietly = TRUE))   install.packages("BiocManager")
# BiocManager::install("DirichletMultinomial")            # ERROR: dependency ‘DirichletMultinomial’ is not available for package ‘leafcutter’
#devtools::install_github("davidaknowles/leafcutter/leafcutter") # ignore update by ZG for R 4.2
library("leafcutter") # use R4.2.0
library(shiny)
library(tidyverse)

# setwd("/home/dhuh/dhuh_edgehpc/SMA/20181210_TST11354_sma_fibroblast_compound_treated/leafcutter") 
fin="/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/LMI070_01x/leafcutter/leafviz.Rdata";

#dir.create("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/leafviz")
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/leafviz") 
# fin="/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2197294_10x/leafcutter/leafviz.Rdata" 


fin="/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2197294_3x/leafcutter/leafviz.Rdata" 
#fin="/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/NGN2_2197294_3x/leafcutter/leafviz.Rdata" 
#fin="/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/Sy5Y_2197294_3x/leafcutter/leafviz.Rdata" 
load(fin)

ls()
annotation_code
cluster_ids
cluster_summary
clusters
code
counts
exons_table
fin
intron_summary 
introns
introns_to_plot
# LIMI070_01x
#      chr    start      end       clu     middle 
# 1   chr2   219001   219747  clu_1_NA   219374.0
meta            
pca             
sample_table   

shiny::runApp("/edgehpc/dept/compbio/users/emarshal/scicomp/public/splicing_work/updated_leafcutter_install/leafcutter/leafviz/")


# compareInfo.csv
CompareName,Subsetting_group,Model,Covariate_levels,Group_name,Group_test,Group_ctrl,Analysis_method,Shrink_logFC,LFC_cutoff
iPSC_1755497_10x_VS_DMSO,,Treatment,,Treatment,iPSC_1755497_10x,iPSC_DMSO,DESeq2,Yes,0.2


# multiqc.sh
conda
cd      /data1/ZG/ATACseq/2020-12-09-Moro-ATAC-Set2/QC/
multiqc /data1/ZG/ATACseq/2020-12-09-Moro-ATAC-Set2/QC/
cd      /data1/ZG/ATACseq/2020-12-09-Moro-ATAC-Set56/QC/
multiqc /data1/ZG/ATACseq/2020-12-09-Moro-ATAC-Set56/QC/
  
https://dl.dnanex.us/F/D/jp5qfZXv6FY3b1fgv179vq4yfj9yvZ66vg58vy15/multiqc_report.html?inline=true#general_stats

if [ "${MODULE_VERSION:-}" = "" ]; then
        MODULE_VERSION_STACK="3.2.10"
        MODULE_VERSION="3.2.10"
        export MODULE_VERSION
else
        MODULE_VERSION_STACK="$MODULE_VERSION"
fi
export MODULE_VERSION_STACK

module() { eval `/camhpc/pkg/Modules/$MODULE_VERSION/bin/modulecmd bash $*`; }
export -f module

MODULESHOME=/camhpc/pkg/Modules/3.2.10
export MODULESHOME

if [ "${LOADEDMODULES:-}" = "" ]; then
  LOADEDMODULES=
  export LOADEDMODULES
fi

if [ "${MODULEPATH:-}" = "" ]; then
  MODULEPATH=`sed -n 's/[       #].*$//; /./H; $ { x; s/^\n//; s/\n/:/g; p; }' ${MODULESHOME}/init/.modulespath | tr -d '\t'`
  export MODULEPATH
fi

if [ ${BASH_VERSINFO:-0} -ge 3 ] && [ -r ${MODULESHOME}/init/bash_completion ]; then
 . ${MODULESHOME}/init/bash_completion
fi

# modules_bash.sh
if [ "${MODULE_VERSION:-}" = "" ]; then
        MODULE_VERSION_STACK="3.2.10"
        MODULE_VERSION="3.2.10"
        export MODULE_VERSION
else
        MODULE_VERSION_STACK="$MODULE_VERSION"
fi
export MODULE_VERSION_STACK

module() { eval `/camhpc/pkg/Modules/$MODULE_VERSION/bin/modulecmd bash $*`; }
export -f module

MODULESHOME=/camhpc/pkg/Modules/3.2.10
export MODULESHOME

if [ "${LOADEDMODULES:-}" = "" ]; then
  LOADEDMODULES=
  export LOADEDMODULES
fi

if [ "${MODULEPATH:-}" = "" ]; then
  MODULEPATH=`sed -n 's/[       #].*$//; /./H; $ { x; s/^\n//; s/\n/:/g; p; }' ${MODULESHOME}/init/.modulespath | tr -d '\t'`
  export MODULEPATH
fi

if [ ${BASH_VERSINFO:-0} -ge 3 ] && [ -r ${MODULESHOME}/init/bash_completion ]; then
 . ${MODULESHOME}/init/bash_completion
fi

# splice result discrepancy 
Zhen Gao
Mon 3/28/2022 3:22 PM
Hi, Dear Dann: I will do it from this afternoon and will report the results and discrepancy (if found) to you. Best, Zhen
Dann Huh
Zhen Gao
Hi Zhen,

I hope you had a good weekend.

I need your help on troubleshooting the splice pipeline that we migrated to Edge cluster from Camhpc.

I see slight discrepancy in results, and I wonder if you could try a few things to confirm this.


Location of scripts and results for my example:
  a.       Camhpc:  
  
  /camhpc/home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing/all5pre_w_and_wo_covar/TST11742_GM24468D-RGX71083-10x/
  
  argparse_wrapper_qsub_leafcutter.bash

argparse_wrapper_qsub_majiq.bash

argparse_wrapper_qsub_rmats.bash

optparse_qsub_make_master_table.sh

b.       Edge:
  
  /edgehpc/dept/compbio/users/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing_Edge_test/TST11742_GM24468D-RGX71083-10x/edge/
  
  argparse_wrapper_qsub_leafcutter_edge_test.bash

argparse_wrapper_qsub_majiq_edge_test.bash

argparse_wrapper_qsub_rmats_edge_test.bash

optparse_qsub_make_master_table_edge_test.sh



as you can see, a. and b. are having the same data as input (please check, though)



Problem:
  rMATS seems to be OK.



However, Leafcutter and Majiq are not returning exact same numbers, e.g.,

Camhpc: chr1,-,WASH7P,,chr1:clu_16300,leafcutter,0.957199121932822,0.792311182184021,-0.00459556493355229,Success,0.44426073959744,0.439665174663887,logef:0.000734201137315743|loglr:0.84561197237781|df:4|N:5|coord:chr1:17055-17915|chr1:17055:17233:clu_16300:-,cryptic,,,,

Edge:       chr1,-,WASH7P,,chr1:clu_16300,leafcutter,0.954724752137984,0.786470436526946,-0.00489946041609252,Success,0.444511085532613,0.43961162511652,logef:0.00367472043072459|loglr:0.861678815403138|df:4|N:5|coord:chr1:17055-17915|chr1:17055:17233:clu_16300:-,cryptic,,,,



Probably this is acceptable overall, but still bothers me if there’s obvious thing that we are missing.



Request:
  Could you please run a test on different dataset that you are comfortable with, both in camhpc and Edge, and see if you see the same discrepancy?
  
  Also, can you run the same set twice in Edge, and see if the results are always the same? (and do the same on camhpc).



Thanks!
  
  
  
  -Dann

# Biogen_automatic_pipeline_ZG_data_11955_record_real2.R
############################################### 0 ###############################################  
https://wiki.biogen.com/pages/viewpage.action?pageId=185896568

which EAinit
/edgehpc/dept/compbio/edge_tools/RNASequest/EAinit


############################################### 1 ###############################################  
cd /edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955


EAinit /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao
loading resource ...


***** 2023-01-22 22:40:14 *****
  ###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
$prj_name
[1] "initPrjName"

$prj_title
[1] "initPrjTitle"

$prj_counts
[1] "initCounts"

$prj_effLength
[1] "initEffLength"

$prj_seqQC
[1] "initSeqQC"

$prj_TPM
[1] "initTPM"

$sample_meta
[1] "initPrjMeta"

$sample_factor
[1] "initPrjFactor"

$sample_name
[1] "Sample_Name"

$sample_alias
NULL

$split_meta
NULL

$species
[1] "initSpecies"

$gene_annotation
[1] "initGeneAnnotation"

$min_count
[1] 1

$min_sample
[1] 1

$count_prior
[1] 0.25

$output
[1] "initOutput"

$DA_file_outpath
[1] "initOutput/DA_Import_Files"

$core
[1] 2

$parallel
[1] FALSE

$qsubTime
[1] 180

$min_median_effective_length
[1] 5

$seqQC_replot
[1] FALSE

$geneLength_QC
[1] FALSE

$rmGeneStartWith
list()

$covariates_check
[1] "initCovariates"

$covariates_check_PCcutoff
[1] 5

$covariates_check_FDRcutoff
[1] 0.1

$covariates_check_plotNcol
[1] 3

$covariates_adjust
list()

$covariates_method
[1] "limma"

$comparison_file
[1] "initPrjComp"

$sample_group
list()

$gene_network_high_variable_N
[1] 10000

$gene_network_cor_cutoff
[1] 0.7

$gene_network_p_cutoff
[1] 0.05

$gene_network_max_edge
[1] "2e6"

$gene_network_min_edge
[1] "2e3"

$shinyOne_Title
NULL

$shinyOne_Description
NULL

$shinyOne_Data_Generated_By
NULL

$shinyOne_Disease_Model
NULL

$shinyOne_Biological_System
NULL

$shinyOne_Sample_Type
NULL

$shinyOne_Technology_Platform
[1] "RNA-Seq"

$shinyOne_Data_Source
[1] "Internal"

$shinyOne_Data_Cleaning
NULL

using /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/analysis-GFJjf6Q0G3b1GQ0GGk7fBFVQ
Create gene annotation ...
/edgehpc/dept/compbio/reference/DNAnexus_references/rnaseq/homo_sapiens/Human.GRCh38.v34.l1_5.ERCC/Human.GRCh38.v34.l1_5.ERCC.transcript.gene_info.csv
getting count file by genes.estcount_table or genes.expected_count
/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/combine_rsem_outputs/genes.expected_count.tsv
Extracting effective length by genes.effective_length
/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/genes.effective_length.tsv
getting seqQC file by combined.metrics or merged_metrics
/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/combine_rnaseqc/merged_metrics.tsv
Appending sequencing QC into sample meta file ...
Creating project folder
Create empty comparison template ...
saving initialization ...
==========================================
  ExpressionAnalysis project folder is created at /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0
-----> 'EAqc' can be used to identify the covariates to be adjusted as:
  EAqc /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml


----->'EArun' can be used to obtain the QuickOmics objects after comparison definition file is updated:
  /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv
EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml


-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.


Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

====== 1h
############################################### 3 ###############################################  

----->'EArun' can be used to obtain the QuickOmics objects after comparison definition file is updated:
  #/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230116_0/data/compareInfo.csv
  EAqc /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
EAqc /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 00:32:29 *****
  ###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
  1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
reading sample meta
checking against config file
reading sample counts
Filtering genes (46064) with minimal counts (>=) 1 in at least (>=) 1 samples
reading effective length
The nominal length of following 921 genes will be used:
  MIR345, MIRLET7F1, MIR302C, MIR148B, MIRLET7D, MIR302D, MIRLET7G, MIRLET7I, AC016394.1, RNU6-834P, RNA5SP46, AC006023.1, RNU6-661P, RNU6-1094P, AL354797.1, RNU5E-1, AL591846.1, RNA5SP432, AC010853.1, RNA5SP453, RNU5F-1, RNY4P19, RNA5SP250, RNA5SP301, AL929554.1, RNA5SP370, SNORD9, AC005520.1, AC105267.1, SNORD115-25, RNU6-212P, AL080243.1, AC021105.1, RNA5SP195, AC121247.1, RNU6-37P, RNU5A-1, RNA5SP174, SNORD114-14, AL353692.1, SNORD33, AL390882.1, RNU6-1330P, RNU6-1272P, RNU6-1266P, SNORD16, SNORD115-2, RNU6-1079P, AC008121.1, SNORD36A, SNORD104, AC016692.1, RNA5SP383, RNA5SP358, RNA5SP375, RNU6-942P, RNA5SP452, AP003499.1, RNU6-888P, AL390195.1, RNY1P16, AP000967.1, VTRNA1-1, AC087638.1, AC025260.1, RNA5SP65, RNA5SP218, AC009065.1, RNY4P9, SNORD68, RNU6-433P, AC021443.1, RNU6-984P, AC023891.1, AC013452.1, AC096649.1, RNU5B-1, AL607035.1, RNU6-238P, AC009309.1, RNU6-529P, RNU6-97P, SNORD35A, AL158824.1, AC006963.1, RNA5SP215, RNU6-485P, AC023356.1, RNU6-1099P, RNA5SP74, SNORA63B, AC090527.1, AC078818.1, RNU6-112P, RNU6-1017P, AL135938.1, SNORD35B, RNU6-137P, RNU6-103P, AC144536.1, RNY1P11, AC006116.1, AC111170.1, RNU6-379P, AC100814.1, SNORD115-8, AL355490.1, RNA5SP161, RNA5SP68, SNORD8, AL445222.1, RNU6-770P, SNORD114-2, AL139113.1, RNU6-82P, AC023105.1, AP003307.1, RNU6-457P, SNORD14E, RNU6-598P, SNORD46, RNA5SP435, AL138751.1, RNU6-501P, SNORD115-32, AC087884.1, RNU5A-8P, AC107390.1, AL049542.1, RNA5SP242, RNU6-268P, AC100840.1, RNY1, RNU6-11P, RNY1P12, AL450992.1, RNU6-353P, RNU6-853P, RNA5SP202, AC089999.2, AC026347.1, RNU6-1305P, AC022080.1, AC060812.2, RNU4-59P, RNA5S9, RNU6-361P, SNORD115-23, RNU6-1251P, RNY4P23, RNU6-558P, AC007375.1, SNORD14B, AL109936.1, RNA5SP141, AC006960.1, AL772161.1, RNA5SP379, RNY4P7, SNORD45B, AL356299.1, RNU6-312P, AC093605.1, RNA5SP513, RNU4-85P, AL662899.2, AC084824.1, AL118556.1, RNU6-343P, RNU6-593P, AL732366.1, RNA5SP84, AC099782.1, AC104828.1, RNY4P34, RNU6-7, SNORD32A, AL136303.1, RNA5SP479, RNU6-828P, RNU6-534P, SNORD52, RNU6-384P, AC104411.1, SNORD117, RNA5SP488, RNA5SP149, SNORD48, RNA5SP298, AC098591.1, AL590434.1, AL732366.2, RNU6-967P, RNY3P1, AL662797.1, RNU6-934P, RNU6-1240P, SNORD38A, RNA5SP152, AL513282.1, AC104335.1, SNORD58C, AC064871.1, RNU6-407P, AC092663.1, RNA5SP128, AL033527.1, AP005061.1, RNU6-1138P, RNA5SP85, SNORD14C, AL031662.1, RNA5SP22, RNA5SP77, SNORD6, RNA5SP366, RNU6-8, RNY3, AC024575.1, RNU6-652P, AC062037.1, RNU5E-6P, AL096840.1, AC084879.1, RNA5SP283, AC116348.1, RNA5SP151, VTRNA1-3, RNU6-395P, AC116366.1, RNU6-1329P, RNU6-354P, RNU6-132P, RNY1P5, SNORD45C, SNORD116-14, RNU6-1, SNORD60, U72788.1, RNU6-1279P, AC006583.1, AL139412.1, SNORD116-17, Z99297.1, RNU6-945P, AC092802.1, SNORD21, AC106802.1, RNU6-1189P, SNORD116-18, AC138969.2, RNU6-26P, RNU6-1056P, SNORD116-9, AC104741.1, AC004461.1, AC022217.1, SNORD101, RNU6-10P, RNU6-116P, RNU6-926P, AC026790.1, AP000751.1, AC008006.1, RNU6-190P, RNU6-48P, AC092610.1, RNU6-36P, RNU6-481P, RNU6-80P, AC126755.2, RNU6-4P, AL133244.1, RNU6-1201P, AC012358.2, RNU6-5P, RNU6-610P, RNU6-574P, SNORD116-2, RNU6-611P, SNORA54, AC016773.1, RNU6-318P, SNORD116-3, AL162740.1, RNU6-975P, AP001453.1, SNORD59A, AL356494.1, AP001011.1, RNU6-3P, SNORD116-1, RNU6-463P, AC010619.1, RNU6-171P, AC073195.1, SNORD116-8, RNU6-31P, SNORD14D, RNA5SP187, SNORD116-7, RNU6-106P, AC079601.1, RNY1P14, RNU6-549P, RNU6-905P, RNU6-1340P, SNORD116-15, RNU6-1157P, RNA5SP122, SNORD116-5, AC131953.1, AP000590.1, RNU6-790P, AC090543.1, AC106897.1, RNU6-125P, AC099850.1, RNU6-1005P, RNU6-342P, SNORD116-16, RNU6-15P, RNA5SP284, RNU6-831P, SNORD116-24, RNU6-30P, AC139256.1, SNORD7, AL033528.1, RNU6-1263P, RNU6-146P, RNU6-680P, RNU6-658P, AL133238.1, AC131235.1, RNU6-2, RNU6-925P, RNU6-14P, RNU6-665P, SNORD116-23, AC012442.1, RNU6-310P, AC027237.2, RNU6-1011P, AP000704.1, AC103952.1, RNU6-813P, AC138932.2, SNORD116-6, RNU6-520P, SNORD116-19, RNU6-9, RNU6-59P, RNU6-33P, AC009053.1, MIR25, MIR217, MIR647, MIR635, MIR598, MIR191, MIR181C, MIR619, MIR7-3, MIR641, MIR570, MIR621, MIR558, MIR592, MIR659, MIR573, MIR597, MIR218-1, MIR590, MIR648, MIR26A2, MIR640, MIR27B, MIR328, MIR140, MIR645, MIR618, SNORD12C, MT-TL1, SNORD83B, SNORD83A, SNORD41, MT-TF, MT-TV, MT-TI, MT-TQ, MT-TM, MT-TA, MT-TN, MT-TC, MT-TY, MT-TS1, MT-TD, MT-TK, MT-TH, MT-TS2, MT-TL2, MT-TE, MT-TT, MIR766, MIR769, MIR765, MIR762, TRAJ24, SNORD67, RNA5SP372, SNORD66, AC096637.1, RNA5SP244, SNORD89, RNU6-821P, SNORD12, RNU6-316P, RNU6-244P, RNU6-1177P, RNU6-482P, SNORD115-45, RNU6-817P, RNY4P36, RNU6-1111P, SNORD90, AC008128.1, RNU6-1158P, SNORD19, SNORD86, RNA5SP300, AL121924.1, RNA5SP212, SNORD70, RNA5SP474, RNA5SP467, SNORA26, MIR1250, MIR1231, SNORD110, MIR1224, AL513534.2, MIR1253, RNU6ATAC27P, SNORD88A, MIR663B, MIR1238, SNORD100, SNORD99, RNU6ATAC10P, RNU6ATAC42P, MIR1249, MIR1225, MIR1282, RNU6-1165P, RNA5SP129, RNA5SP431, RNA5SP201, RNU6-101P, RNU6-757P, RNA5SP425, RNA5SP193, RNA5SP118, RNA5SP237, SNORD12B, AC015563.1, RNU6-554P, RNA5SP511, AC073520.1, RNU6-1077P, AL121929.1, RNU6-519P, RNU6-705P, RNA5SP60, RNA5SP345, AC063943.1, RNU2-42P, RNU6-1267P, RNU6-1190P, RNU4-14P, RNU6-321P, RNA5SP21, RNU6-920P, AC080112.1, RNA5SP104, RNA5SP265, RNA5SP344, RNU6-130P, AC024580.1, RNU6-1245P, RNA5SP450, AL392105.1, SNORD71, RNA5SP508, RNA5SP492, RNU6-195P, RNA5SP107, RNA5SP270, AC073261.1, RPL35AP4, AC137499.1, AC079150.1, SNORD57, AC005326.1, RPL41P1, AL162726.3, AC244023.1, SNORD56, AC123900.1, AC093162.1, Z96811.1, AC233266.1, RPS29P23, SNORD62B, AL161935.2, AF228730.1, AL021068.1, SNORD62A, AC113174.1, AL627311.1, GAGE12B, AL353691.2, SNORD121B, AC122179.1, SNORD13P1, SNORD4A, AL133211.1, SNORD42A, AC096757.1, SNORD124, SNORD121A, AC004912.1, RNU7-1, RNA5SP246, SNORD13, AC069200.2, SNORD127, SNORD125, AC091849.1, AC108729.3, AC018645.1, AC010598.1, RNU2-13P, AC004849.1, RNU6-1053P, RNU6-999P, RNA5SP217, RNU6-377P, RNU6-202P, RNU6-322P, RNA5SP162, RNU6-513P, RNA5SP216, RNU6-1143P, RNA5SP180, RNY4P37, RNA5SP464, RNU6-415P, RNU6-469P, RNY3P15, AC002543.2, RNU6-313P, RNU6-1045P, AC090227.1, RNY3P16, RNA5SP318, RNU6-781P, AC008393.2, RNU6-764P, RNA5SP92, RNA5SP61, RNU4ATAC12P, AC079173.1, RNY4, SNORD116-25, AC022916.1, RNU6-503P, RNU6-1061P, AC005722.2, RNU6ATAC12P, RNU6-118P, RNA5SP438, RNU6-358P, RNA5SP395, RNU6-1004P, RNU6-100P, AC115088.1, AC092567.2, RNU6-405P, RNA5SP168, RNA5SP434, RNU6-126P, RNA5SP33, RNU6-1016P, RNU6-531P, RNA5SP82, RNU6-965P, RNA5SP466, RNU6-759P, RNU6-307P, RNU6-914P, AC009509.1, RNU6-807P, RNA5SP441, RNA5SP481, RNU6-826P, RNA5SP268, RNU6-1136P, AC012640.3, RNU6-579P, RNU6-850P, RNA5SP143, RNU6-143P, RNU6-703P, RNU6-577P, AC005037.2, RNU6-341P, AC104446.1, RN7SKP233, RNU6-388P, RNU2-46P, RNU6-570P, AC021443.2, RNU6-548P, RNU6-1003P, AL139099.1, RNU6-1223P, RNA5SP340, AC016959.1, RNA5SP70, RNU6-667P, RNU6-902P, AL356356.1, RNA5SP373, AL049840.3, RNU6-323P, AC007991.2, AC024995.1, SNORD87, RNY1P9, RPL41P5, AC068831.3, CYP4A43P, MIR3187, MIR3936, AC008670.1, MIR4737, MIR4688, AC016596.3, SNORD55, MIR4433B, MIR4653, MIR4656, AC090616.4, SNORD95, MIR23C, MIR5194, MIR5191, MIR4441, AC099677.4, MIR5094, MIR4301, MIR4312, SNORD84, MIR3935, MIR3195, MIR4263, MIR5188, AC016601.1, MIR4502, MIR4664, SNORD53B, MIR4292, MIR3685, MIR4519, MIR4659A, MIR4440, MIR1260B, MIR5091, MIR744, MIR3665, MIR4273, MIR3934, MIR4297, MIR3945, MIR3661, MIR4482, MIR212, AC011447.4, AC008737.2, AL096840.2, RPL12P50, AC242498.1, RNU6-94P, MIR98, RNA5SP108, SNORD58B, SNORD14A, RNU6-6P, AC103810.6, AC002395.1, AC103810.7, SNORD96A, RNU6-90P, AC103810.8, MIR6129, AC005562.2, MIR6792, AL731566.1, MIR6765, AP003499.3, AC004386.3, AP000446.2, SNORD1C, AC007027.1, RNA5SP439, AC007351.3, MIR6864, AC004079.1, AC026318.1, MIR1273H, AL132655.4, MIR6832, AC012354.2, SNORD28, MIR378J, AL157834.3, AL049766.1, AL133351.3, AL096828.4, MIR6750, AL133467.5, SNORD25, SNORD50B, AC211486.7, SNORD116-22, AC243994.1, AP003499.4, MIR6740, AC007351.4, AL121875.1, AC233263.2, SNORD116-4, MIR6505, AC084368.1, AL669914.2, AL021155.1, MIR6807, AC211486.8, Z99774.2, AC108065.1, AC004542.4, AL356356.2, AC242852.1, AC211476.9, AC087284.2, AC015726.2, RNU6-467P, SNORD64, MIR6769A, AL445423.2, AL033538.2, MIR7111, AL132655.6, AC243829.3, SNORD26, AL132655.7, RNA5SP530, SNORD49B, AC013734.1, AL021938.1, SNORD65, AC009014.3, AC106873.7, MIR7845, AC009812.5, MIR6856, MIR6719, SNORD116-21, SNORD30, MIR8075, AL049766.3, AC010332.1, AC073476.4, AL139099.4, AC117415.2, AP000944.3, AL669914.3, AP000944.4, AC108065.2, AP000769.5, AC084123.1, MIR6892, AC002542.4, MIR7161, MIR6080, MIR6077, MIR6775, SNORD116-20, AC114498.2, AP001362.1, AC067863.2, AC005480.2, AC090498.1, AC006435.5, RNA5SP196, AC002059.2, AC009065.10, MIR3651, AC007383.4, SNORD38B, AL356095.3, MIR3120, MIR4256, MIR6841, MIR219B, AL358790.2, FO624990.1, AC024933.2, AL139128.2, MIR4454, MIR6759, MIR5087, AC010487.2, MIR555, MIR6831, MIR6884, MIR5193, MIR6755, MIR1178, MIR622, MIR612, MIR4453, MIR3198-2, MIR636, MIR7844, MIR639, MIR378I, MIR4751, MIR214, MIR4721, MIR130B, MIR3911, MIR3916, MIR6132, MIR637, MIR5010, MIR423, MIR600, MIR4709, MIR3610, MIR1248, MIR4442, MIR6084, MIR4697, MIR6834, MIR1199, MIR6073, MIR611, MIR198, MIR6789, MIR4700, MIR922, MIR4784, MIR6758, MIR3605, MIR7-1, MIR6506, MIR671, MIR1257, MIR658, MIR29C, MIR2682, MIR6516, MIR5004, MIR4683, MIR4321, MIR4800, MIR939, MIR124-1, MIR3655, MIR4741, MIR940, MIR142, MIR5047, MIR147B, MIR5001, MIR3191, MIR663A, MIR4724, MIR4647, MIR2110, MIR1236, MIR6501, MIR1304, MIR1306, MIR4720, MIR564, MIR4687, MIR4723, MIR3614, MIR2861, MIR3064, MIR7705, MIR4793, MIR5006, AP000553.4, AP000553.7, AC073869.6, RNA5SP343, AC119676.1, AC135721.2, AC022506.1, AC026740.2, AC011994.1, AC104805.2, AC104843.3, AC023049.1, AL359092.3, AC015871.7, AC007599.3
reading sequence QC
Finished TPM estimation!
  plotting sequencing QC @/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/sequenceQC.pdf
Please check output files at: /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0 
-----> PC analysis without covariate adjusted:
  /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/notAdjusted
Warning message:
  In plotPCanlaysis(config, D$logTPM, D$meta, estC = D$counts, effL = D$effLength) :
  < covariates_adjust is NOT set in the config file, no covariate was adjusted! >
  ==========================================
  ----->'EArun' can be used to obtain the QuickOmics object after necessary 'covariates_adjust' is set and comparison definition file is filled:
  /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv
EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml


-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.

Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

====5min


############################################### 3 ############################################### 
#EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230116_0/config.yml

cp 	/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/compareInfo.csv /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv
EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml


$  EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 00:43:43 *****
  ###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
  1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
Warning message:
  In readLines(config[[one]]) :
  incomplete final line found on '/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv'
reading sample meta
checking against config file
Error in checkComparisonInfo(read_file(config$comparison_file, T), D$meta,  : 
                               The DEG 'Model' didn't include 'Group_name' in the comparison file for iPSC_4088_10x_vs_iPSC_DMSO
Calls: getEAData -> checkComparisonInfo
Execution halted

##  change the Group_name  value from "iPSC_4088_10x" to 
Treatment
Treatment
Treatment
.
.
.

EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 01:08:19 *****
###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
Warning message:
In readLines(config[[one]]) :
  incomplete final line found on '/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv'
reading sample meta
checking against config file
Error in checkComparisonInfo(read_file(config$comparison_file, T), D$meta,  : 
  /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv should have column header(s): Group_test
Calls: getEAData -> checkComparisonInfo
Execution halted

### accidently changed the value back to "Group_test", early "Group-test" when keep the Group in sample sheet "iPSC-4088-10x"
(Prepaired by Joon)

/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao zgao1@edge-hpc-log-101 $ 
$  EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 01:09:31 *****
###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
Warning message:
In readLines(config[[one]]) :
  incomplete final line found on '/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv'
reading sample meta
checking against config file
Error in checkComparisonInfo(read_file(config$comparison_file, T), D$meta,  : 
  'Group_test' entry iPSC-4088-10x for comparison iPSC_4088_10x_vs_iPSC_DMSO contains characters other than letters, numbers, and delimiters '_' or '.'. Please use only letters, numbers, '_' or '.', as these are safe characters for column names in R.
Calls: getEAData -> checkComparisonInfo
Execution halted


update samplesheet.tsv column "Treatment" from  "iPSC-4088-10x" to "iPSC_4088_10x"

############################################### 3 ###############################################  
 EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_3/config.yml
loading resource ...


***** 2023-01-23 07:04:23 *****
###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
Warning message:
In readLines(config[[one]]) :
  incomplete final line found on '/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_3/data/compareInfo.csv'
reading sample meta
checking against config file
	Checking model for iPSC_4088_10x_vs_iPSC_DMSO
Error in checkComparisonModel(comp_info, meta) : 
  Error in iPSC_4088_10x_vs_iPSC_DMSO : ctrl_group  iPSC_DMSO  is NOT defined in the sample meta file.
Calls: getEAData -> checkComparisonInfo -> checkComparisonModel
Execution halted


iPSC_DMSO_1A
iPSC_DMSO_2A
iPSC_DMSO_3A

should changed to 

iPSC_DMSO
iPSC_DMSO
iPSC_DMSO

############################################### 3 ###############################################  
 EAqc /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 20:37:47 *****
###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
reading sample meta
checking against config file
reading sample counts
	Filtering genes (46064) with minimal counts (>=) 1 in at least (>=) 1 samples
reading effective length
		The nominal length of following 921 genes will be used:
		MIR345, MIRLET7F1, MIR302C, MIR148B, MIRLET7D, MIR302D, MIRLET7G, MIRLET7I, AC016394.1, RNU6-834P, RNA5SP46, AC006023.1, RNU6-661P, RNU6-1094P, AL354797.1, RNU5E-1, AL591846.1, RNA5SP432, AC010853.1, RNA5SP453, RNU5F-1, RNY4P19, RNA5SP250, RNA5SP301, AL929554.1, RNA5SP370, SNORD9, AC005520.1, AC105267.1, SNORD115-25, RNU6-212P, AL080243.1, AC021105.1, RNA5SP195, AC121247.1, RNU6-37P, RNU5A-1, RNA5SP174, SNORD114-14, AL353692.1, SNORD33, AL390882.1, RNU6-1330P, RNU6-1272P, RNU6-1266P, SNORD16, SNORD115-2, RNU6-1079P, AC008121.1, SNORD36A, SNORD104, AC016692.1, RNA5SP383, RNA5SP358, RNA5SP375, RNU6-942P, RNA5SP452, AP003499.1, RNU6-888P, AL390195.1, RNY1P16, AP000967.1, VTRNA1-1, AC087638.1, AC025260.1, RNA5SP65, RNA5SP218, AC009065.1, RNY4P9, SNORD68, RNU6-433P, AC021443.1, RNU6-984P, AC023891.1, AC013452.1, AC096649.1, RNU5B-1, AL607035.1, RNU6-238P, AC009309.1, RNU6-529P, RNU6-97P, SNORD35A, AL158824.1, AC006963.1, RNA5SP215, RNU6-485P, AC023356.1, RNU6-1099P, RNA5SP74, SNORA63B, AC090527.1, AC078818.1, RNU6-112P, RNU6-1017P, AL135938.1, SNORD35B, RNU6-137P, RNU6-103P, AC144536.1, RNY1P11, AC006116.1, AC111170.1, RNU6-379P, AC100814.1, SNORD115-8, AL355490.1, RNA5SP161, RNA5SP68, SNORD8, AL445222.1, RNU6-770P, SNORD114-2, AL139113.1, RNU6-82P, AC023105.1, AP003307.1, RNU6-457P, SNORD14E, RNU6-598P, SNORD46, RNA5SP435, AL138751.1, RNU6-501P, SNORD115-32, AC087884.1, RNU5A-8P, AC107390.1, AL049542.1, RNA5SP242, RNU6-268P, AC100840.1, RNY1, RNU6-11P, RNY1P12, AL450992.1, RNU6-353P, RNU6-853P, RNA5SP202, AC089999.2, AC026347.1, RNU6-1305P, AC022080.1, AC060812.2, RNU4-59P, RNA5S9, RNU6-361P, SNORD115-23, RNU6-1251P, RNY4P23, RNU6-558P, AC007375.1, SNORD14B, AL109936.1, RNA5SP141, AC006960.1, AL772161.1, RNA5SP379, RNY4P7, SNORD45B, AL356299.1, RNU6-312P, AC093605.1, RNA5SP513, RNU4-85P, AL662899.2, AC084824.1, AL118556.1, RNU6-343P, RNU6-593P, AL732366.1, RNA5SP84, AC099782.1, AC104828.1, RNY4P34, RNU6-7, SNORD32A, AL136303.1, RNA5SP479, RNU6-828P, RNU6-534P, SNORD52, RNU6-384P, AC104411.1, SNORD117, RNA5SP488, RNA5SP149, SNORD48, RNA5SP298, AC098591.1, AL590434.1, AL732366.2, RNU6-967P, RNY3P1, AL662797.1, RNU6-934P, RNU6-1240P, SNORD38A, RNA5SP152, AL513282.1, AC104335.1, SNORD58C, AC064871.1, RNU6-407P, AC092663.1, RNA5SP128, AL033527.1, AP005061.1, RNU6-1138P, RNA5SP85, SNORD14C, AL031662.1, RNA5SP22, RNA5SP77, SNORD6, RNA5SP366, RNU6-8, RNY3, AC024575.1, RNU6-652P, AC062037.1, RNU5E-6P, AL096840.1, AC084879.1, RNA5SP283, AC116348.1, RNA5SP151, VTRNA1-3, RNU6-395P, AC116366.1, RNU6-1329P, RNU6-354P, RNU6-132P, RNY1P5, SNORD45C, SNORD116-14, RNU6-1, SNORD60, U72788.1, RNU6-1279P, AC006583.1, AL139412.1, SNORD116-17, Z99297.1, RNU6-945P, AC092802.1, SNORD21, AC106802.1, RNU6-1189P, SNORD116-18, AC138969.2, RNU6-26P, RNU6-1056P, SNORD116-9, AC104741.1, AC004461.1, AC022217.1, SNORD101, RNU6-10P, RNU6-116P, RNU6-926P, AC026790.1, AP000751.1, AC008006.1, RNU6-190P, RNU6-48P, AC092610.1, RNU6-36P, RNU6-481P, RNU6-80P, AC126755.2, RNU6-4P, AL133244.1, RNU6-1201P, AC012358.2, RNU6-5P, RNU6-610P, RNU6-574P, SNORD116-2, RNU6-611P, SNORA54, AC016773.1, RNU6-318P, SNORD116-3, AL162740.1, RNU6-975P, AP001453.1, SNORD59A, AL356494.1, AP001011.1, RNU6-3P, SNORD116-1, RNU6-463P, AC010619.1, RNU6-171P, AC073195.1, SNORD116-8, RNU6-31P, SNORD14D, RNA5SP187, SNORD116-7, RNU6-106P, AC079601.1, RNY1P14, RNU6-549P, RNU6-905P, RNU6-1340P, SNORD116-15, RNU6-1157P, RNA5SP122, SNORD116-5, AC131953.1, AP000590.1, RNU6-790P, AC090543.1, AC106897.1, RNU6-125P, AC099850.1, RNU6-1005P, RNU6-342P, SNORD116-16, RNU6-15P, RNA5SP284, RNU6-831P, SNORD116-24, RNU6-30P, AC139256.1, SNORD7, AL033528.1, RNU6-1263P, RNU6-146P, RNU6-680P, RNU6-658P, AL133238.1, AC131235.1, RNU6-2, RNU6-925P, RNU6-14P, RNU6-665P, SNORD116-23, AC012442.1, RNU6-310P, AC027237.2, RNU6-1011P, AP000704.1, AC103952.1, RNU6-813P, AC138932.2, SNORD116-6, RNU6-520P, SNORD116-19, RNU6-9, RNU6-59P, RNU6-33P, AC009053.1, MIR25, MIR217, MIR647, MIR635, MIR598, MIR191, MIR181C, MIR619, MIR7-3, MIR641, MIR570, MIR621, MIR558, MIR592, MIR659, MIR573, MIR597, MIR218-1, MIR590, MIR648, MIR26A2, MIR640, MIR27B, MIR328, MIR140, MIR645, MIR618, SNORD12C, MT-TL1, SNORD83B, SNORD83A, SNORD41, MT-TF, MT-TV, MT-TI, MT-TQ, MT-TM, MT-TA, MT-TN, MT-TC, MT-TY, MT-TS1, MT-TD, MT-TK, MT-TH, MT-TS2, MT-TL2, MT-TE, MT-TT, MIR766, MIR769, MIR765, MIR762, TRAJ24, SNORD67, RNA5SP372, SNORD66, AC096637.1, RNA5SP244, SNORD89, RNU6-821P, SNORD12, RNU6-316P, RNU6-244P, RNU6-1177P, RNU6-482P, SNORD115-45, RNU6-817P, RNY4P36, RNU6-1111P, SNORD90, AC008128.1, RNU6-1158P, SNORD19, SNORD86, RNA5SP300, AL121924.1, RNA5SP212, SNORD70, RNA5SP474, RNA5SP467, SNORA26, MIR1250, MIR1231, SNORD110, MIR1224, AL513534.2, MIR1253, RNU6ATAC27P, SNORD88A, MIR663B, MIR1238, SNORD100, SNORD99, RNU6ATAC10P, RNU6ATAC42P, MIR1249, MIR1225, MIR1282, RNU6-1165P, RNA5SP129, RNA5SP431, RNA5SP201, RNU6-101P, RNU6-757P, RNA5SP425, RNA5SP193, RNA5SP118, RNA5SP237, SNORD12B, AC015563.1, RNU6-554P, RNA5SP511, AC073520.1, RNU6-1077P, AL121929.1, RNU6-519P, RNU6-705P, RNA5SP60, RNA5SP345, AC063943.1, RNU2-42P, RNU6-1267P, RNU6-1190P, RNU4-14P, RNU6-321P, RNA5SP21, RNU6-920P, AC080112.1, RNA5SP104, RNA5SP265, RNA5SP344, RNU6-130P, AC024580.1, RNU6-1245P, RNA5SP450, AL392105.1, SNORD71, RNA5SP508, RNA5SP492, RNU6-195P, RNA5SP107, RNA5SP270, AC073261.1, RPL35AP4, AC137499.1, AC079150.1, SNORD57, AC005326.1, RPL41P1, AL162726.3, AC244023.1, SNORD56, AC123900.1, AC093162.1, Z96811.1, AC233266.1, RPS29P23, SNORD62B, AL161935.2, AF228730.1, AL021068.1, SNORD62A, AC113174.1, AL627311.1, GAGE12B, AL353691.2, SNORD121B, AC122179.1, SNORD13P1, SNORD4A, AL133211.1, SNORD42A, AC096757.1, SNORD124, SNORD121A, AC004912.1, RNU7-1, RNA5SP246, SNORD13, AC069200.2, SNORD127, SNORD125, AC091849.1, AC108729.3, AC018645.1, AC010598.1, RNU2-13P, AC004849.1, RNU6-1053P, RNU6-999P, RNA5SP217, RNU6-377P, RNU6-202P, RNU6-322P, RNA5SP162, RNU6-513P, RNA5SP216, RNU6-1143P, RNA5SP180, RNY4P37, RNA5SP464, RNU6-415P, RNU6-469P, RNY3P15, AC002543.2, RNU6-313P, RNU6-1045P, AC090227.1, RNY3P16, RNA5SP318, RNU6-781P, AC008393.2, RNU6-764P, RNA5SP92, RNA5SP61, RNU4ATAC12P, AC079173.1, RNY4, SNORD116-25, AC022916.1, RNU6-503P, RNU6-1061P, AC005722.2, RNU6ATAC12P, RNU6-118P, RNA5SP438, RNU6-358P, RNA5SP395, RNU6-1004P, RNU6-100P, AC115088.1, AC092567.2, RNU6-405P, RNA5SP168, RNA5SP434, RNU6-126P, RNA5SP33, RNU6-1016P, RNU6-531P, RNA5SP82, RNU6-965P, RNA5SP466, RNU6-759P, RNU6-307P, RNU6-914P, AC009509.1, RNU6-807P, RNA5SP441, RNA5SP481, RNU6-826P, RNA5SP268, RNU6-1136P, AC012640.3, RNU6-579P, RNU6-850P, RNA5SP143, RNU6-143P, RNU6-703P, RNU6-577P, AC005037.2, RNU6-341P, AC104446.1, RN7SKP233, RNU6-388P, RNU2-46P, RNU6-570P, AC021443.2, RNU6-548P, RNU6-1003P, AL139099.1, RNU6-1223P, RNA5SP340, AC016959.1, RNA5SP70, RNU6-667P, RNU6-902P, AL356356.1, RNA5SP373, AL049840.3, RNU6-323P, AC007991.2, AC024995.1, SNORD87, RNY1P9, RPL41P5, AC068831.3, CYP4A43P, MIR3187, MIR3936, AC008670.1, MIR4737, MIR4688, AC016596.3, SNORD55, MIR4433B, MIR4653, MIR4656, AC090616.4, SNORD95, MIR23C, MIR5194, MIR5191, MIR4441, AC099677.4, MIR5094, MIR4301, MIR4312, SNORD84, MIR3935, MIR3195, MIR4263, MIR5188, AC016601.1, MIR4502, MIR4664, SNORD53B, MIR4292, MIR3685, MIR4519, MIR4659A, MIR4440, MIR1260B, MIR5091, MIR744, MIR3665, MIR4273, MIR3934, MIR4297, MIR3945, MIR3661, MIR4482, MIR212, AC011447.4, AC008737.2, AL096840.2, RPL12P50, AC242498.1, RNU6-94P, MIR98, RNA5SP108, SNORD58B, SNORD14A, RNU6-6P, AC103810.6, AC002395.1, AC103810.7, SNORD96A, RNU6-90P, AC103810.8, MIR6129, AC005562.2, MIR6792, AL731566.1, MIR6765, AP003499.3, AC004386.3, AP000446.2, SNORD1C, AC007027.1, RNA5SP439, AC007351.3, MIR6864, AC004079.1, AC026318.1, MIR1273H, AL132655.4, MIR6832, AC012354.2, SNORD28, MIR378J, AL157834.3, AL049766.1, AL133351.3, AL096828.4, MIR6750, AL133467.5, SNORD25, SNORD50B, AC211486.7, SNORD116-22, AC243994.1, AP003499.4, MIR6740, AC007351.4, AL121875.1, AC233263.2, SNORD116-4, MIR6505, AC084368.1, AL669914.2, AL021155.1, MIR6807, AC211486.8, Z99774.2, AC108065.1, AC004542.4, AL356356.2, AC242852.1, AC211476.9, AC087284.2, AC015726.2, RNU6-467P, SNORD64, MIR6769A, AL445423.2, AL033538.2, MIR7111, AL132655.6, AC243829.3, SNORD26, AL132655.7, RNA5SP530, SNORD49B, AC013734.1, AL021938.1, SNORD65, AC009014.3, AC106873.7, MIR7845, AC009812.5, MIR6856, MIR6719, SNORD116-21, SNORD30, MIR8075, AL049766.3, AC010332.1, AC073476.4, AL139099.4, AC117415.2, AP000944.3, AL669914.3, AP000944.4, AC108065.2, AP000769.5, AC084123.1, MIR6892, AC002542.4, MIR7161, MIR6080, MIR6077, MIR6775, SNORD116-20, AC114498.2, AP001362.1, AC067863.2, AC005480.2, AC090498.1, AC006435.5, RNA5SP196, AC002059.2, AC009065.10, MIR3651, AC007383.4, SNORD38B, AL356095.3, MIR3120, MIR4256, MIR6841, MIR219B, AL358790.2, FO624990.1, AC024933.2, AL139128.2, MIR4454, MIR6759, MIR5087, AC010487.2, MIR555, MIR6831, MIR6884, MIR5193, MIR6755, MIR1178, MIR622, MIR612, MIR4453, MIR3198-2, MIR636, MIR7844, MIR639, MIR378I, MIR4751, MIR214, MIR4721, MIR130B, MIR3911, MIR3916, MIR6132, MIR637, MIR5010, MIR423, MIR600, MIR4709, MIR3610, MIR1248, MIR4442, MIR6084, MIR4697, MIR6834, MIR1199, MIR6073, MIR611, MIR198, MIR6789, MIR4700, MIR922, MIR4784, MIR6758, MIR3605, MIR7-1, MIR6506, MIR671, MIR1257, MIR658, MIR29C, MIR2682, MIR6516, MIR5004, MIR4683, MIR4321, MIR4800, MIR939, MIR124-1, MIR3655, MIR4741, MIR940, MIR142, MIR5047, MIR147B, MIR5001, MIR3191, MIR663A, MIR4724, MIR4647, MIR2110, MIR1236, MIR6501, MIR1304, MIR1306, MIR4720, MIR564, MIR4687, MIR4723, MIR3614, MIR2861, MIR3064, MIR7705, MIR4793, MIR5006, AP000553.4, AP000553.7, AC073869.6, RNA5SP343, AC119676.1, AC135721.2, AC022506.1, AC026740.2, AC011994.1, AC104805.2, AC104843.3, AC023049.1, AL359092.3, AC015871.7, AC007599.3
reading sequence QC
Finished TPM estimation!
plotting sequencing QC @/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/sequenceQC.pdf
Please check output files at: /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0 
-----> PC analysis without covariate adjusted:
	/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/notAdjusted
Warning message:
In plotPCanlaysis(config, D$logTPM, D$meta, estC = D$counts, effL = D$effLength) :
  < covariates_adjust is NOT set in the config file, no covariate was adjusted! >
==========================================
----->'EArun' can be used to obtain the QuickOmics object after necessary 'covariates_adjust' is set and comparison definition file is filled:
	/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv
		EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml


-----> (additional) 'EAsplit' can be used to split into sub-project according to one column (split_meta) defined in the sample meta file.

Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

############################################### 3 ###############################################  

		cp	/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/compareInfo.csv 			/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv
		cp	/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao//compareInfo.csv 		
############################################### 3 ###############################################  
 EArun /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...


***** 2023-01-23 23:07:00 *****
###########
## ExpressionAnalysis: https://github.com/interactivereport/RNASequest.git
## Pipeline Path: /mnt/depts/dept04/compbio/edge_tools/RNASequest
## Pipeline Date: 2023-01-20 23:45:59 -0500
## git HEAD: e472604b20c76438641a00c00ae817686189e6ee
###########

Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
Warning message:
In readLines(config[[one]]) :
  incomplete final line found on '/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/data/compareInfo.csv'
reading sample meta
checking against config file
	Checking model for iPSC_4088_10x_vs_iPSC_DMSO
	Checking model for iPSC_4088_3x_vs_iPSC_DMSO
	Checking model for iPSC_4090_10x_vs_iPSC_DMSO
	Checking model for iPSC_4090_3x_vs_iPSC_DMSO
	Checking model for iPSC_4714_10x_vs_iPSC_DMSO
	Checking model for iPSC_4714_3x_vs_iPSC_DMSO
	Checking model for iPSC_4741_10x_vs_iPSC_DMSO
	Checking model for iPSC_4741_3x_vs_iPSC_DMSO
	Checking model for iPSC_4748_10x_vs_iPSC_DMSO
	Checking model for iPSC_4748_3x_vs_iPSC_DMSO
	Checking model for iPSC_5420_10x_vs_iPSC_DMSO
	Checking model for iPSC_5420_3x_vs_iPSC_DMSO
	Checking model for iPSC_6152_10x_vs_iPSC_DMSO
	Checking model for iPSC_6152_3x_vs_iPSC_DMSO
	Checking model for iPSC_6527_10x_vs_iPSC_DMSO
	Checking model for iPSC_6527_3x_vs_iPSC_DMSO
	Checking model for iPSC_6866_10x_vs_iPSC_DMSO
	Checking model for iPSC_6866_3x_vs_iPSC_DMSO
	Checking model for iPSC_6960_10x_vs_iPSC_DMSO
	Checking model for iPSC_6960_3x_vs_iPSC_DMSO
	Checking model for iPSC_9701_10x_vs_iPSC_DMSO
	Checking model for iPSC_9701_3x_vs_iPSC_DMSO
	Checking model for NGN2_4088_10x_vs_NGN2_DMSO
	Checking model for NGN2_4088_3x_vs_NGN2_DMSO
	Checking model for NGN2_4090_10x_vs_NGN2_DMSO
	Checking model for NGN2_4090_3x_vs_NGN2_DMSO
	Checking model for NGN2_4714_10x_vs_NGN2_DMSO
	Checking model for NGN2_4714_3x_vs_NGN2_DMSO
	Checking model for NGN2_4741_10x_vs_NGN2_DMSO
	Checking model for NGN2_4741_3x_vs_NGN2_DMSO
	Checking model for NGN2_4748_10x_vs_NGN2_DMSO
	Checking model for NGN2_4748_3x_vs_NGN2_DMSO
	Checking model for NGN2_5420_10x_vs_NGN2_DMSO
	Checking model for NGN2_5420_3x_vs_NGN2_DMSO
	Checking model for NGN2_6152_10x_vs_NGN2_DMSO
	Checking model for NGN2_6152_3x_vs_NGN2_DMSO
	Checking model for NGN2_6527_10x_vs_NGN2_DMSO
	Checking model for NGN2_6527_3x_vs_NGN2_DMSO
	Checking model for NGN2_6866_10x_vs_NGN2_DMSO
	Checking model for NGN2_6866_3x_vs_NGN2_DMSO
	Checking model for NGN2_6960_10x_vs_NGN2_DMSO
	Checking model for NGN2_6960_3x_vs_NGN2_DMSO
	Checking model for NGN2_9701_10x_vs_NGN2_DMSO
	Checking model for NGN2_9701_3x_vs_NGN2_DMSO
reading sample counts
	Filtering genes (46064) with minimal counts (>=) 1 in at least (>=) 1 samples
reading effective length
		The nominal length of following 921 genes will be used:
		MIR345, MIRLET7F1, MIR302C, MIR148B, MIRLET7D, MIR302D, MIRLET7G, MIRLET7I, AC016394.1, RNU6-834P, RNA5SP46, AC006023.1, RNU6-661P, RNU6-1094P, AL354797.1, RNU5E-1, AL591846.1, RNA5SP432, AC010853.1, RNA5SP453, RNU5F-1, RNY4P19, RNA5SP250, RNA5SP301, AL929554.1, RNA5SP370, SNORD9, AC005520.1, AC105267.1, SNORD115-25, RNU6-212P, AL080243.1, AC021105.1, RNA5SP195, AC121247.1, RNU6-37P, RNU5A-1, RNA5SP174, SNORD114-14, AL353692.1, SNORD33, AL390882.1, RNU6-1330P, RNU6-1272P, RNU6-1266P, SNORD16, SNORD115-2, RNU6-1079P, AC008121.1, SNORD36A, SNORD104, AC016692.1, RNA5SP383, RNA5SP358, RNA5SP375, RNU6-942P, RNA5SP452, AP003499.1, RNU6-888P, AL390195.1, RNY1P16, AP000967.1, VTRNA1-1, AC087638.1, AC025260.1, RNA5SP65, RNA5SP218, AC009065.1, RNY4P9, SNORD68, RNU6-433P, AC021443.1, RNU6-984P, AC023891.1, AC013452.1, AC096649.1, RNU5B-1, AL607035.1, RNU6-238P, AC009309.1, RNU6-529P, RNU6-97P, SNORD35A, AL158824.1, AC006963.1, RNA5SP215, RNU6-485P, AC023356.1, RNU6-1099P, RNA5SP74, SNORA63B, AC090527.1, AC078818.1, RNU6-112P, RNU6-1017P, AL135938.1, SNORD35B, RNU6-137P, RNU6-103P, AC144536.1, RNY1P11, AC006116.1, AC111170.1, RNU6-379P, AC100814.1, SNORD115-8, AL355490.1, RNA5SP161, RNA5SP68, SNORD8, AL445222.1, RNU6-770P, SNORD114-2, AL139113.1, RNU6-82P, AC023105.1, AP003307.1, RNU6-457P, SNORD14E, RNU6-598P, SNORD46, RNA5SP435, AL138751.1, RNU6-501P, SNORD115-32, AC087884.1, RNU5A-8P, AC107390.1, AL049542.1, RNA5SP242, RNU6-268P, AC100840.1, RNY1, RNU6-11P, RNY1P12, AL450992.1, RNU6-353P, RNU6-853P, RNA5SP202, AC089999.2, AC026347.1, RNU6-1305P, AC022080.1, AC060812.2, RNU4-59P, RNA5S9, RNU6-361P, SNORD115-23, RNU6-1251P, RNY4P23, RNU6-558P, AC007375.1, SNORD14B, AL109936.1, RNA5SP141, AC006960.1, AL772161.1, RNA5SP379, RNY4P7, SNORD45B, AL356299.1, RNU6-312P, AC093605.1, RNA5SP513, RNU4-85P, AL662899.2, AC084824.1, AL118556.1, RNU6-343P, RNU6-593P, AL732366.1, RNA5SP84, AC099782.1, AC104828.1, RNY4P34, RNU6-7, SNORD32A, AL136303.1, RNA5SP479, RNU6-828P, RNU6-534P, SNORD52, RNU6-384P, AC104411.1, SNORD117, RNA5SP488, RNA5SP149, SNORD48, RNA5SP298, AC098591.1, AL590434.1, AL732366.2, RNU6-967P, RNY3P1, AL662797.1, RNU6-934P, RNU6-1240P, SNORD38A, RNA5SP152, AL513282.1, AC104335.1, SNORD58C, AC064871.1, RNU6-407P, AC092663.1, RNA5SP128, AL033527.1, AP005061.1, RNU6-1138P, RNA5SP85, SNORD14C, AL031662.1, RNA5SP22, RNA5SP77, SNORD6, RNA5SP366, RNU6-8, RNY3, AC024575.1, RNU6-652P, AC062037.1, RNU5E-6P, AL096840.1, AC084879.1, RNA5SP283, AC116348.1, RNA5SP151, VTRNA1-3, RNU6-395P, AC116366.1, RNU6-1329P, RNU6-354P, RNU6-132P, RNY1P5, SNORD45C, SNORD116-14, RNU6-1, SNORD60, U72788.1, RNU6-1279P, AC006583.1, AL139412.1, SNORD116-17, Z99297.1, RNU6-945P, AC092802.1, SNORD21, AC106802.1, RNU6-1189P, SNORD116-18, AC138969.2, RNU6-26P, RNU6-1056P, SNORD116-9, AC104741.1, AC004461.1, AC022217.1, SNORD101, RNU6-10P, RNU6-116P, RNU6-926P, AC026790.1, AP000751.1, AC008006.1, RNU6-190P, RNU6-48P, AC092610.1, RNU6-36P, RNU6-481P, RNU6-80P, AC126755.2, RNU6-4P, AL133244.1, RNU6-1201P, AC012358.2, RNU6-5P, RNU6-610P, RNU6-574P, SNORD116-2, RNU6-611P, SNORA54, AC016773.1, RNU6-318P, SNORD116-3, AL162740.1, RNU6-975P, AP001453.1, SNORD59A, AL356494.1, AP001011.1, RNU6-3P, SNORD116-1, RNU6-463P, AC010619.1, RNU6-171P, AC073195.1, SNORD116-8, RNU6-31P, SNORD14D, RNA5SP187, SNORD116-7, RNU6-106P, AC079601.1, RNY1P14, RNU6-549P, RNU6-905P, RNU6-1340P, SNORD116-15, RNU6-1157P, RNA5SP122, SNORD116-5, AC131953.1, AP000590.1, RNU6-790P, AC090543.1, AC106897.1, RNU6-125P, AC099850.1, RNU6-1005P, RNU6-342P, SNORD116-16, RNU6-15P, RNA5SP284, RNU6-831P, SNORD116-24, RNU6-30P, AC139256.1, SNORD7, AL033528.1, RNU6-1263P, RNU6-146P, RNU6-680P, RNU6-658P, AL133238.1, AC131235.1, RNU6-2, RNU6-925P, RNU6-14P, RNU6-665P, SNORD116-23, AC012442.1, RNU6-310P, AC027237.2, RNU6-1011P, AP000704.1, AC103952.1, RNU6-813P, AC138932.2, SNORD116-6, RNU6-520P, SNORD116-19, RNU6-9, RNU6-59P, RNU6-33P, AC009053.1, MIR25, MIR217, MIR647, MIR635, MIR598, MIR191, MIR181C, MIR619, MIR7-3, MIR641, MIR570, MIR621, MIR558, MIR592, MIR659, MIR573, MIR597, MIR218-1, MIR590, MIR648, MIR26A2, MIR640, MIR27B, MIR328, MIR140, MIR645, MIR618, SNORD12C, MT-TL1, SNORD83B, SNORD83A, SNORD41, MT-TF, MT-TV, MT-TI, MT-TQ, MT-TM, MT-TA, MT-TN, MT-TC, MT-TY, MT-TS1, MT-TD, MT-TK, MT-TH, MT-TS2, MT-TL2, MT-TE, MT-TT, MIR766, MIR769, MIR765, MIR762, TRAJ24, SNORD67, RNA5SP372, SNORD66, AC096637.1, RNA5SP244, SNORD89, RNU6-821P, SNORD12, RNU6-316P, RNU6-244P, RNU6-1177P, RNU6-482P, SNORD115-45, RNU6-817P, RNY4P36, RNU6-1111P, SNORD90, AC008128.1, RNU6-1158P, SNORD19, SNORD86, RNA5SP300, AL121924.1, RNA5SP212, SNORD70, RNA5SP474, RNA5SP467, SNORA26, MIR1250, MIR1231, SNORD110, MIR1224, AL513534.2, MIR1253, RNU6ATAC27P, SNORD88A, MIR663B, MIR1238, SNORD100, SNORD99, RNU6ATAC10P, RNU6ATAC42P, MIR1249, MIR1225, MIR1282, RNU6-1165P, RNA5SP129, RNA5SP431, RNA5SP201, RNU6-101P, RNU6-757P, RNA5SP425, RNA5SP193, RNA5SP118, RNA5SP237, SNORD12B, AC015563.1, RNU6-554P, RNA5SP511, AC073520.1, RNU6-1077P, AL121929.1, RNU6-519P, RNU6-705P, RNA5SP60, RNA5SP345, AC063943.1, RNU2-42P, RNU6-1267P, RNU6-1190P, RNU4-14P, RNU6-321P, RNA5SP21, RNU6-920P, AC080112.1, RNA5SP104, RNA5SP265, RNA5SP344, RNU6-130P, AC024580.1, RNU6-1245P, RNA5SP450, AL392105.1, SNORD71, RNA5SP508, RNA5SP492, RNU6-195P, RNA5SP107, RNA5SP270, AC073261.1, RPL35AP4, AC137499.1, AC079150.1, SNORD57, AC005326.1, RPL41P1, AL162726.3, AC244023.1, SNORD56, AC123900.1, AC093162.1, Z96811.1, AC233266.1, RPS29P23, SNORD62B, AL161935.2, AF228730.1, AL021068.1, SNORD62A, AC113174.1, AL627311.1, GAGE12B, AL353691.2, SNORD121B, AC122179.1, SNORD13P1, SNORD4A, AL133211.1, SNORD42A, AC096757.1, SNORD124, SNORD121A, AC004912.1, RNU7-1, RNA5SP246, SNORD13, AC069200.2, SNORD127, SNORD125, AC091849.1, AC108729.3, AC018645.1, AC010598.1, RNU2-13P, AC004849.1, RNU6-1053P, RNU6-999P, RNA5SP217, RNU6-377P, RNU6-202P, RNU6-322P, RNA5SP162, RNU6-513P, RNA5SP216, RNU6-1143P, RNA5SP180, RNY4P37, RNA5SP464, RNU6-415P, RNU6-469P, RNY3P15, AC002543.2, RNU6-313P, RNU6-1045P, AC090227.1, RNY3P16, RNA5SP318, RNU6-781P, AC008393.2, RNU6-764P, RNA5SP92, RNA5SP61, RNU4ATAC12P, AC079173.1, RNY4, SNORD116-25, AC022916.1, RNU6-503P, RNU6-1061P, AC005722.2, RNU6ATAC12P, RNU6-118P, RNA5SP438, RNU6-358P, RNA5SP395, RNU6-1004P, RNU6-100P, AC115088.1, AC092567.2, RNU6-405P, RNA5SP168, RNA5SP434, RNU6-126P, RNA5SP33, RNU6-1016P, RNU6-531P, RNA5SP82, RNU6-965P, RNA5SP466, RNU6-759P, RNU6-307P, RNU6-914P, AC009509.1, RNU6-807P, RNA5SP441, RNA5SP481, RNU6-826P, RNA5SP268, RNU6-1136P, AC012640.3, RNU6-579P, RNU6-850P, RNA5SP143, RNU6-143P, RNU6-703P, RNU6-577P, AC005037.2, RNU6-341P, AC104446.1, RN7SKP233, RNU6-388P, RNU2-46P, RNU6-570P, AC021443.2, RNU6-548P, RNU6-1003P, AL139099.1, RNU6-1223P, RNA5SP340, AC016959.1, RNA5SP70, RNU6-667P, RNU6-902P, AL356356.1, RNA5SP373, AL049840.3, RNU6-323P, AC007991.2, AC024995.1, SNORD87, RNY1P9, RPL41P5, AC068831.3, CYP4A43P, MIR3187, MIR3936, AC008670.1, MIR4737, MIR4688, AC016596.3, SNORD55, MIR4433B, MIR4653, MIR4656, AC090616.4, SNORD95, MIR23C, MIR5194, MIR5191, MIR4441, AC099677.4, MIR5094, MIR4301, MIR4312, SNORD84, MIR3935, MIR3195, MIR4263, MIR5188, AC016601.1, MIR4502, MIR4664, SNORD53B, MIR4292, MIR3685, MIR4519, MIR4659A, MIR4440, MIR1260B, MIR5091, MIR744, MIR3665, MIR4273, MIR3934, MIR4297, MIR3945, MIR3661, MIR4482, MIR212, AC011447.4, AC008737.2, AL096840.2, RPL12P50, AC242498.1, RNU6-94P, MIR98, RNA5SP108, SNORD58B, SNORD14A, RNU6-6P, AC103810.6, AC002395.1, AC103810.7, SNORD96A, RNU6-90P, AC103810.8, MIR6129, AC005562.2, MIR6792, AL731566.1, MIR6765, AP003499.3, AC004386.3, AP000446.2, SNORD1C, AC007027.1, RNA5SP439, AC007351.3, MIR6864, AC004079.1, AC026318.1, MIR1273H, AL132655.4, MIR6832, AC012354.2, SNORD28, MIR378J, AL157834.3, AL049766.1, AL133351.3, AL096828.4, MIR6750, AL133467.5, SNORD25, SNORD50B, AC211486.7, SNORD116-22, AC243994.1, AP003499.4, MIR6740, AC007351.4, AL121875.1, AC233263.2, SNORD116-4, MIR6505, AC084368.1, AL669914.2, AL021155.1, MIR6807, AC211486.8, Z99774.2, AC108065.1, AC004542.4, AL356356.2, AC242852.1, AC211476.9, AC087284.2, AC015726.2, RNU6-467P, SNORD64, MIR6769A, AL445423.2, AL033538.2, MIR7111, AL132655.6, AC243829.3, SNORD26, AL132655.7, RNA5SP530, SNORD49B, AC013734.1, AL021938.1, SNORD65, AC009014.3, AC106873.7, MIR7845, AC009812.5, MIR6856, MIR6719, SNORD116-21, SNORD30, MIR8075, AL049766.3, AC010332.1, AC073476.4, AL139099.4, AC117415.2, AP000944.3, AL669914.3, AP000944.4, AC108065.2, AP000769.5, AC084123.1, MIR6892, AC002542.4, MIR7161, MIR6080, MIR6077, MIR6775, SNORD116-20, AC114498.2, AP001362.1, AC067863.2, AC005480.2, AC090498.1, AC006435.5, RNA5SP196, AC002059.2, AC009065.10, MIR3651, AC007383.4, SNORD38B, AL356095.3, MIR3120, MIR4256, MIR6841, MIR219B, AL358790.2, FO624990.1, AC024933.2, AL139128.2, MIR4454, MIR6759, MIR5087, AC010487.2, MIR555, MIR6831, MIR6884, MIR5193, MIR6755, MIR1178, MIR622, MIR612, MIR4453, MIR3198-2, MIR636, MIR7844, MIR639, MIR378I, MIR4751, MIR214, MIR4721, MIR130B, MIR3911, MIR3916, MIR6132, MIR637, MIR5010, MIR423, MIR600, MIR4709, MIR3610, MIR1248, MIR4442, MIR6084, MIR4697, MIR6834, MIR1199, MIR6073, MIR611, MIR198, MIR6789, MIR4700, MIR922, MIR4784, MIR6758, MIR3605, MIR7-1, MIR6506, MIR671, MIR1257, MIR658, MIR29C, MIR2682, MIR6516, MIR5004, MIR4683, MIR4321, MIR4800, MIR939, MIR124-1, MIR3655, MIR4741, MIR940, MIR142, MIR5047, MIR147B, MIR5001, MIR3191, MIR663A, MIR4724, MIR4647, MIR2110, MIR1236, MIR6501, MIR1304, MIR1306, MIR4720, MIR564, MIR4687, MIR4723, MIR3614, MIR2861, MIR3064, MIR7705, MIR4793, MIR5006, AP000553.4, AP000553.7, AC073869.6, RNA5SP343, AC119676.1, AC135721.2, AC022506.1, AC026740.2, AC011994.1, AC104805.2, AC104843.3, AC023049.1, AL359092.3, AC015871.7, AC007599.3
reading sequence QC
Finished TPM estimation!
Plotting correlation between covarites and comparison groups:
	Index_ID, RIN, Well_Row, Sapio_Plate_Name, Well_Column v.s. Treatment
	Treatment .vs. Index_ID
	Treatment .vs. RIN
	Treatment .vs. Well_Row
	Treatment .vs. Sapio_Plate_Name
	Treatment .vs. Well_Column
Warning message:
In anova.lm(lm(as.formula(paste(c(y, x), collapse = "~")), data = D)) :
  ANOVA F-tests on an essentially perfect fit are unreliable
====== Starting DE analysis ...
serial DEG process ...
Comparison for iPSC_4088_10x_vs_iPSC_DMSO; iPSC_4088_3x_vs_iPSC_DMSO; iPSC_4090_10x_vs_iPSC_DMSO; iPSC_4090_3x_vs_iPSC_DMSO; iPSC_4714_10x_vs_iPSC_DMSO; iPSC_4714_3x_vs_iPSC_DMSO; iPSC_4741_10x_vs_iPSC_DMSO; iPSC_4741_3x_vs_iPSC_DMSO; iPSC_4748_10x_vs_iPSC_DMSO; iPSC_4748_3x_vs_iPSC_DMSO; iPSC_5420_10x_vs_iPSC_DMSO; iPSC_5420_3x_vs_iPSC_DMSO; iPSC_6152_10x_vs_iPSC_DMSO; iPSC_6152_3x_vs_iPSC_DMSO; iPSC_6527_10x_vs_iPSC_DMSO; iPSC_6527_3x_vs_iPSC_DMSO; iPSC_6866_10x_vs_iPSC_DMSO; iPSC_6866_3x_vs_iPSC_DMSO; iPSC_6960_10x_vs_iPSC_DMSO; iPSC_6960_3x_vs_iPSC_DMSO; iPSC_9701_10x_vs_iPSC_DMSO; iPSC_9701_3x_vs_iPSC_DMSO; NGN2_4088_10x_vs_NGN2_DMSO; NGN2_4088_3x_vs_NGN2_DMSO; NGN2_4090_10x_vs_NGN2_DMSO; NGN2_4090_3x_vs_NGN2_DMSO; NGN2_4714_10x_vs_NGN2_DMSO; NGN2_4714_3x_vs_NGN2_DMSO; NGN2_4741_10x_vs_NGN2_DMSO; NGN2_4741_3x_vs_NGN2_DMSO; NGN2_4748_10x_vs_NGN2_DMSO; NGN2_4748_3x_vs_NGN2_DMSO; NGN2_5420_10x_vs_NGN2_DMSO; NGN2_5420_3x_vs_NGN2_DMSO; NGN2_6152_10x_vs_NGN2_DMSO; NGN2_6152_3x_vs_NGN2_DMSO; NGN2_6527_10x_vs_NGN2_DMSO; NGN2_6527_3x_vs_NGN2_DMSO; NGN2_6866_10x_vs_NGN2_DMSO; NGN2_6866_3x_vs_NGN2_DMSO; NGN2_6960_10x_vs_NGN2_DMSO; NGN2_6960_3x_vs_NGN2_DMSO; NGN2_9701_10x_vs_NGN2_DMSO; NGN2_9701_3x_vs_NGN2_DMSO
converting counts to integer mode
estimating size factors
estimating dispersions
gene-wise dispersion estimates: 2 workers
mean-dispersion relationship
final dispersion estimates, fitting model and testing: 2 workers
	--- iPSC_4088_10x_vs_iPSC_DMSO
found results columns, replacing these
		relevel
	--- iPSC_4088_3x_vs_iPSC_DMSO
	--- iPSC_4090_10x_vs_iPSC_DMSO
	--- iPSC_4090_3x_vs_iPSC_DMSO
	--- iPSC_4714_10x_vs_iPSC_DMSO
	--- iPSC_4714_3x_vs_iPSC_DMSO
	--- iPSC_4741_10x_vs_iPSC_DMSO
	--- iPSC_4741_3x_vs_iPSC_DMSO
	--- iPSC_4748_10x_vs_iPSC_DMSO
	--- iPSC_4748_3x_vs_iPSC_DMSO
	--- iPSC_5420_10x_vs_iPSC_DMSO
	--- iPSC_5420_3x_vs_iPSC_DMSO
	--- iPSC_6152_10x_vs_iPSC_DMSO
	--- iPSC_6152_3x_vs_iPSC_DMSO
	--- iPSC_6527_10x_vs_iPSC_DMSO
	--- iPSC_6527_3x_vs_iPSC_DMSO
	--- iPSC_6866_10x_vs_iPSC_DMSO
	--- iPSC_6866_3x_vs_iPSC_DMSO
	--- iPSC_6960_10x_vs_iPSC_DMSO
	--- iPSC_6960_3x_vs_iPSC_DMSO
	--- iPSC_9701_10x_vs_iPSC_DMSO
	--- iPSC_9701_3x_vs_iPSC_DMSO
	--- NGN2_4088_10x_vs_NGN2_DMSO
found results columns, replacing these
		relevel
	--- NGN2_4088_3x_vs_NGN2_DMSO
	--- NGN2_4090_10x_vs_NGN2_DMSO
	--- NGN2_4090_3x_vs_NGN2_DMSO
	--- NGN2_4714_10x_vs_NGN2_DMSO
	--- NGN2_4714_3x_vs_NGN2_DMSO
	--- NGN2_4741_10x_vs_NGN2_DMSO
	--- NGN2_4741_3x_vs_NGN2_DMSO
	--- NGN2_4748_10x_vs_NGN2_DMSO
	--- NGN2_4748_3x_vs_NGN2_DMSO
	--- NGN2_5420_10x_vs_NGN2_DMSO
	--- NGN2_5420_3x_vs_NGN2_DMSO
	--- NGN2_6152_10x_vs_NGN2_DMSO
	--- NGN2_6152_3x_vs_NGN2_DMSO
	--- NGN2_6527_10x_vs_NGN2_DMSO
	--- NGN2_6527_3x_vs_NGN2_DMSO
	--- NGN2_6866_10x_vs_NGN2_DMSO
	--- NGN2_6866_3x_vs_NGN2_DMSO
	--- NGN2_6960_10x_vs_NGN2_DMSO
	--- NGN2_6960_3x_vs_NGN2_DMSO
	--- NGN2_9701_10x_vs_NGN2_DMSO
	--- NGN2_9701_3x_vs_NGN2_DMSO
     user    system   elapsed 
18389.751   268.139 10088.410 
There were 44 warnings (use warnings() to see them)
Obtaining networks ...
   user  system elapsed 
 54.524   4.348  59.017 
saving QuickOmics object ...
	Formating the expression data
	Formating the DEG results
	Formating the sample meta information
	saving...
=================================================
Results are saved in /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0

-----> Please visit: http://ngs.biogen.com:3838/Quickomics/?testfile=TST11955
Please carefully review the results before publishing:
----->'EApub' can be used to publish the project into ShinyOne project manager: http://ngs.biogen.com/shinyone

Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

############################################### 3 ###############################################  
 EApub /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0/config.yml
loading resource ...
Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
preparing information for ShinyOne
submitting to ShinyOne manager ...
=================================================
ShinyOne access: http://ngs.biogen.com/shinyone/app/core/app_project_review.php?ID=237

Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

############################################### 3 ###############################################  
$  EApub /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA20230123_0/config.yml
loading resource ...
Warning messages:
1: package ‘R.oo’ was built under R version 4.1.3 
2: package ‘R.methodsS3’ was built under R version 4.1.3 
3: package ‘matrixStats’ was built under R version 4.1.3 
preparing information for ShinyOne
submitting to ShinyOne manager ...
=================================================
ShinyOne access: http://ngs.biogen.com/shinyone/app/core/app_project_review.php?ID=238

Powered by the Research Data Sciences Group [fergal.casey@biogen.com;zhengyu.ouyang@biogen.com]

############################################### 3 ###############################################  


# sgailong 
module load anaconda3/4.9.2
conda activate /edgehpc/dept/compbio/users/scao/software/anaconda3/r_knnwgcna2
jupyter lab --no-browser --ip=`/usr/bin/hostname -I | /usr/bin/grep -oE "10.106.192.[0-9]{1,3}"`

# conda yml file
#sample_environments.yml

# cd /edgehpc/dept/compbio/users/zgao1/Coding_system_tool_utilities/conda 
# module load anaconda3
# conda create -n ggsashimi_py3  --file=sample_environments.yml
# conda env create -f environment.yml
name: ggsashimi_py3
channels:      ## Add any other channels below if necessary
- conda-forge
- bioconda
dependencies:  ## Prioritize conda packages
- python=3.10
- jupyter
- conda
- mamba
- ipython
- ipykernel
- numpy
- matplotlib
- scipy
- pandas
- pip
- pre-commit
- black
- nbstripout
- mypy
- flake8
- pycodestyle
- pydocstyle
- pytest
- pytest-cov
- pytest-xdist
- pysam
- r-base=4.2.0
- r-ggplot2
- r-data.table
- pip:  ## Add in pip packages if necessary
  - mkdocs
  - mkdocs-material
  - mkdocstrings
  - mknotebooks


# conda create -n ggsashimi_py3 
# Kejie Li <kejie.li@biogen.com>
# Baohong Zhang;DL-Research Data Sciences
# Zhengyu Ouyang <oyoung@bioinforx.com>;Mehool Patel;Eric Marshall;Lingyu Zhou
# Trick is to use mamba update not create.
# conda install -c conda-forge mamba
# conda create -n <your envname>
# mamba env update -n <your envname> -f VIP_conda_R.yml
# conda activate <your envname>
# 
# Best,
# 
# -Kejie


file locations for SMN2 program RASLseq panel

Hi Zhen,
Thanks for helping on this project. The aim is to select 20 genes and the cooridantes for RASLseq probe for SMN2 off-targets.
We should focus on fibroblast data set (20190409..pptx attached), Risdiplam-treated cells, 03x and 10x EC50 (please ignore 30x).
You can pick the genes with highest abs(dPSI) with FDR<0.05 or P(dPSI)>0.9.

Also cross-reference with Branaplam (LMI070) data, and you’re likely to find something overlapping.

Also cross-reference with iPSC-MN data (TST11451..pptx attached), but this is not absolutely necessary.


We will discuss with tox group if there’s other considerations.

Location of data files:

Fibroblast: this one does not have the “master table”, but separate results are at:
                                /edgehpc/dept/compbio/users/dhuh/SMASM/20181210_TST11354_sma_fibroblast_compound_treated/
                                               majiq/res.03.deltapsi_voila_thrdpsi_0.3/
                                               rmats/
                                               Leafcutter/results_m1_p0.001/
Compound names:  
“RG7916” = Risdiplam
“LMI070” = Branaplam
“907272” = analogous compound from a chemical library

                                Also, DSG lists I complied can be found at (if it helps):
/edgehpc/dept/compbio/users/dhuh/SMASM/20181210_TST11354_sma_fibroblast_compound_treated/analysis_02_compare_combine_lc_rm_mj/DSG_list.xlsx   ,     tab “DSG_per_dose_compound”

                iPSC-MN:

/edgehpc/dept/compbio/users/dhuh/SMASM/20190913_TST11451_sma_iPSCMN_compound_treated/master_table/

res_comb_GM24468D-RG7961*/master_all.txt    # “GM..” is the name of the cell line used.
res_comb_HNDS0030_01-RG7961*/master_all.txt
res_comb_HNDS006_01C2-RG7961*/master_all.txt

                                DSG that’s overlapping in all three cell lines:
/edgehpc/dept/compbio/users/dhuh/SMASM/20190913_TST11451_sma_iPSCMN_compound_treated/analysis_03_DSG/analysis_02_combined_DSG/res.02.compare_plot_DSG_abThr/DSGoverlap_genes_per_dose.txt

Please don’t hesitate to email me if you have any questions.
Have a nice weekend!
-regards,
Dann


========= 2023-03-02 updated!
Dann Huh
Zhen Gao
Hi Zhen,
I sorted out a few things, and if you are searching for the following data, please use the new location:
TST11354_sma_fibroblast_compound_treated -> /edgehpc/dept/compbio/projects/TST11354/
TST11451_sma_iPSCMN_compound_treated -> /edgehpc/dept/compbio/projects/TST11451
-regards,
Dann


## mv_files_to_subfolder1.sh
#! /bin/bash

# https://bash-intro.rsquaredacademy.com/r-command-line.html
# ===> In RStudio, commands can be executed from shell scripts by pressing Ctrl + Enter , 
# ===> OR just select a line and click the Run the right corne 
cd /home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.01.filter_events_dPSI_padj/
mkdir iPSC_dPSI_0.1/ iPSC_dPSI_0.2/ iPSC_dPSI_0.3/ iPSC_dPSI_0.4/
mkdir NGN2_dPSI_0.1/ NGN2_dPSI_0.2/ NGN2_dPSI_0.3/ NGN2_dPSI_0.4/

#for f in aula-??.?.mp4 ; do # https://unix.stackexchange.com/questions/527073/move-files-to-specific-folders-based-on-name filteredevents_TST11872_NGN2-BIO-2060573-3uM_dPSI0.4_padj0.05_M50_EdPSI0.?_PdPSI0.9.mp4

# ls     filteredevents_TST11872_NGN2-BIO-2060573-3uM_dPSI0.*_padj0.05_M50_EdPSI0.*_PdPSI0.9.csv 

for f in filteredevents_TST11872_iPSC-BIO-*_dPSI0.?_padj0.05_M50_EdPSI0.?_PdPSI0.9.csv ; do
    if   [[ $f == *0.1_PdPSI0.9.csv ]] ; then   # * is used for pattern matching 
        echo $f 
        mv "$f"   iPSC_dPSI_0.1/
    elif [[ $f == *0.2_PdPSI0.9.csv ]] ; then   # * is used for pattern matching 
        echo $f 
        mv "$f"   iPSC_dPSI_0.2/
    elif [[ $f == *0.3_PdPSI0.9.csv ]] ; then   # * is used for pattern matching
        echo $f 
        mv "$f"   iPSC_dPSI_0.3/
    elif [[ $f == *0.4_PdPSI0.9.csv ]] ; then   # * is used for pattern matching
        echo $f 
        mv "$f"  iPSC_dPSI_0.4/       #    optional ===> else # NO "then"" here <========== the format
    fi
done


for f in filteredevents_TST11872_NGN2-BIO-*_dPSI0.?_padj0.05_M50_EdPSI0.?_PdPSI0.9.csv ; do
    if   [[ $f == *0.1_PdPSI0.9.csv ]] ; then   # * is used for pattern matching 
        echo $f 
        mv "$f"   NGN2_dPSI_0.1/
    elif [[ $f == *0.2_PdPSI0.9.csv ]] ; then   # * is used for pattern matching 
        echo $f 
        mv "$f"   NGN2_dPSI_0.2/
    elif [[ $f == *0.3_PdPSI0.9.csv ]] ; then   # * is used for pattern matching
        echo $f 
        mv "$f"   NGN2_dPSI_0.3/
    elif [[ $f == *0.4_PdPSI0.9.csv ]] ; then   # * is used for pattern matching
        echo $f 
        mv "$f"  NGN2_dPSI_0.4/       #    optional ===> else # NO "then"" here <========== the format
    fi
done


# modifly a genome

# modify a genome 
# https://fastr.biogen.com/resource
# Manage Resources Section  ==> rnaseq_genome_fasta ==> Download 
# type ==> modify ; species ==> human

## ======= rnaseq_genome_fasta
File_Name                             type   species                  version                 catalogue               description               uploader_email              Creation_Timestamp
hg19as_SMN1_masked.ref.fa.gz	      modified	human	        hg19.gencode.v28lift37.noSMN1	    rnaseq_genome_fasta	   human genome fasta	  svc-DNANexusSciComp@biogen.com     2/13/2019, 11:54:26 AM


## ====== rnaseq_rsem_reference
File_Name                             type   species                  version                 catalogue               description                 uploader_email             Creation_Timestamp
hg19as_SMN1_masked.ref.RSEM.tar.gz	modified	human	        hg19.gencode.v28lift37.noSMN1  	  rnaseq_rsem_reference          null	          svc-DNANexusSciComp@biogen.com	  10/8/2021, 5:03:40 PM

## ======rnaseq_star_index
File_Name                             type   species                  version                 catalogue          description                 uploader_email             Creation_Timestamp
Human.GRCh38.v34.l1_5.ERCC.transcript.STMN2_ce.star-index.tar.gz	modified	human	Human.GRCh38.v34.l1_5.ERCC.transcript.STMN2_ce	rnaseq_star_index	generated by official workflow	lingyu.zhou@biogen.com	4/15/2022, 11:36:46 AM
## ======

## When Run RNAseq, the only selection is:
hg19.gencode.v28lift37.noSMN1

/edgehpc/dept/compbio/genomes/human/gencode.v34

GeneCards: 
https://www.genecards.org/cgi-bin/carddisp.pl?gene=SMN1&keywords=SMN1

2023-08-15

Latest Assembly
chr5:70,924,941-70,966,375
(GRCh38/hg38)
Size: 41,435 bases Orientation: Plus strand

seqkit seq Human.GRCh38.r34.l1_5.ERCC.fa
seqkit seq -i Human.GRCh38.r34.l1_5.ERCC.fa
seqkit seq stat Human.GRCh38.r34.l1_5.ERCC.fa
qkit seq stats Human.GRCh38.r34.l1_5.ERCC.fa
seqkit seq Human.GRCh38.r34.l1_5.ERCC.fa -n
seqkit seq Human.GRCh38.r34.l1_5.ERCC.fa -n -i


cat hsa.fa | seqkit mutate -d 1:5 -s chr1,chr2 -o hsa_del.fa


cat Human.GRCh38.r34.l1_5.ERCC.fa | seqkit stats
file  format  type  num_seqs        sum_len  min_len       avg_len      max_len
-     FASTA   DNA        286  3,099,833,474      273  10,838,578.6  248,956,422

cat Human.GRCh38.r34.l1_5.ERCC.fa | seqkit mutate -d 70924941:70966375 -s chr5 -o Human.GRCh38.r34.l1_5.ERCC_noSMN1.fa #####################


-rwxrwxr-- 1 zgao1 zgao1 3151503000 Aug 16 09:22 Human.GRCh38.r34.l1_5.ERCC.fa
-rw-rw-r-- 1 zgao1 zgao1 3151460606 Aug 16 09:45 Human.GRCh38.r34.l1_5.ERCC_noSMN1.fa



cp -R /edgehpc/dept/compbio/genomes/human/gencode.v34/fasta/ /edgehpc/dept/compbio/users/zgao1/Documents/Genome/
cp -R /edgehpc/dept/compbio/genomes/human/gencode.v34/rsem/ /edgehpc/dept/compbio/users/zgao1/Documents/Genome/
cd /edgehpc/dept/compbio/users/zgao1/Documents/Genome/
mv fasta/ fasta_raw/
cp -R fasta_raw/ fasta_smn1_deleted/

  
head Human.GRCh38.r34.l1_5.ERCC.fa
>chr1 1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
..
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

tail Human.GRCh38.r34.l1_5.ERCC.fa
GGACGGGCACGCTCATATCAGGCTATATTTGGTCCGGGTTATTATCGTCG

/
>chr5
>chr5 5
14661008 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

# google "delete gene sequence start and end from whole genome fasta file"

# shenwei.me
# https://bioinf.shenwei.me › seqkit › usage
# SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. ... insertion, deletion) pair match up paired-end reads from two fastq files ...

# https://bioinf.shenwei.me/seqkit/usage/

# SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation
# 
# Version: 2.5.0
# 
# Author: Wei Shen <shenwei356@gmail.com>
#   
#   Documents  : http://bioinf.shenwei.me/seqkit
# Source code: https://github.com/shenwei356/seqkit
# Please cite: https://doi.org/10.1371/journal.pone.0163962
# 
# 
# Seqkit utlizies the pgzip (https://github.com/klauspost/pgzip) package to
# read and write gzip file, and the outputted gzip file would be slighty
# larger than files generated by GNU gzip.
# 
# Seqkit writes gzip files very fast, much faster than the multi-threaded pigz,
# therefore there's no need to pipe the result to gzip/pigz.
# 
# Seqkit also supports reading and writing xz (.xz) and zstd (.zst) formats since v2.2.0.
# Bzip2 format is supported since v2.4.0.
# 
# Available Commands:
#   amplicon        extract amplicon (or specific region around it) via primer(s)
#   ....
#   grep            search sequences by ID/name/sequence/sequence motifs, mismatch allowed
#   head            print first N FASTA/Q records
#   head-genome     print sequences of the first genome with common prefixes in name
#   locate          locate subsequences/motifs, mismatch allowed
#   merge-slides    merge sliding windows generated from seqkit sliding
#   mutate          edit sequence (point mutation, insertion, deletion)
#
#

#=== mutate
Usage
edit sequence (point mutation, insertion, deletion)

#
# https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta
# This is what samtools faidx is intended for.
# 
# It needs to be called twice.

module load anaconda3
conda activate base_ZG
conda install -c "bioconda/label/cf201901" seqkit

Collecting package metadata (current_repodata.json): done
Solving environment: done


==> WARNING: A newer version of conda exists. <==
current version: 22.11.1
latest version: 23.7.2

## Package Plan ##

environment location: /home/zgao1/.conda/envs/base_ZG

added / updated specs:
  - seqkit

The following packages will be downloaded:
  
  package                    |            build
---------------------------|-----------------
  seqkit-0.10.0              |                1         2.6 MB  bioconda/label/cf201901
------------------------------------------------------------
  Total:         2.6 MB

The following NEW packages will be INSTALLED:
  
  seqkit             bioconda/label/cf201901/linux-64::seqkit-0.10.0-1 

Proceed ([y]/n)? y
Downloading and Extracting Packages
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
############################################################################################### 

cat hsa.fa | seqkit mutate -p -1:X -s chr1,chr2
cat hsa.fa | seqkit mutate -p -1:X -s chr1,chr2 -o hsa_append_X.fa
#cp hsa.fa hsa_back.fa
cat hsa.fa | seqkit mutate -p -1:X -s chr1,chr2 -o hsa.fa
cat hsa.fa | seqkit mutate -d 1:5 -s chr1,chr2 -o hsa_del.fa

echo -ne ">1\nACTGNactgn\n>2\nactgnACTGN\n" | seqkit mutate -d -3:-1 -s 2

echo -ne ">1\nACTGNactgn\n>2\nactgnACTGN\n" | seqkit mutate -d -3:-1 -s 2
echo -ne ">1\nACTGNactgn\n>2\nactgnACTGN\n" | seqkit mutate -d -3:-1 --quiet 


############################################################################################### 
How to Modify Reference Sequence and Annotation Files

https://www.biostars.org/p/9516986/
  
#   if it's only for this one time (or limited number of times) best and most easy is to do this manually.
# 
# Open your sequence & annotation file in a genome-editing tool, eg. jBrowse, Apollo, IGV(?), artemis, Genomeview, ...
# make your edit to the genome
# save your fasta and annotation file
# export in desired output format.
# Normally if you make the edit the annotations will change accordingly (and especially so if you change it in a non genic region, genic regions might be a bit more cumbersome).
# 
# Perhaps not all of those browsers I mentioned do allow to edit the genome sequence though (this is rather a special feature of some of those) but among that list there should be some that can do this.
# 

qsub issue:


#optparse options
in_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/INDIR"
#; mkdir -p $in_dir      # CHANGE EVERY TIME ======<<<<<<
out_dir=$in_dir
out_prefix="INDIR"
majiq_cutoffs="0.1,0.2,0.3,0.4"

#script location
#prog_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/"
prog_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/"

#submit job
echo "source /etc/profile.d/modules_bash.sh; module load xz/5.2.2; module load hdf5/1.8.17; module load R/3.6.1; Rscript $prog_dir/sub.make_master_table_v2_ZG.R --in_dir $in_dir --out_dir $out_dir --out_prefix $out_prefix --majiq_cutoffs $majiq_cutoffs" | qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V


qsub -N out.master_${out_prefix} -q cpu.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V



qstat -u zgao1 # http://web.mit.edu/longjobs/www/status.html


qstat
job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
------------------------------------------------------------------------------------------------------------------------------------------------
   4556838 0.57333 out.majiq  zgao1        r     03/01/2022 17:32:33 all.q@camhpcpe01                                                 12        
   4556839 0.57333 out.majiq  zgao1        r     03/01/2022 17:33:01 all.q@camhpcpe26.hpc.biogen.co                                   12        
   4556859 0.57333 out.majiq  zgao1        r     03/01/2022 20:12:49 all.q@camhpcpe32.hpc.biogen.co                                   12        
   4556860 0.57333 out.majiq  zgao1        r     03/01/2022 20:15:50 all.q@camhpcpe31.hpc.biogen.co                                   12        
   4556861 0.57333 out.majiq  zgao1        r     03/01/2022 20:22:12 all.q@camhpcpe22.hpc.biogen.co                                   12        
   4556862 0.57333 out.majiq  zgao1        r     03/01/2022 20:28:27 all.q@camhpcpe10.hpc.biogen.co                                   12        
   
qstat -j 4556867
#Following jobs do not exist or permissions are not sufficient: 
#4556867
#
 qstat -j 4556838
==============================================================
job_number:                 4556838
jclass:                     NONE
exec_file:                  job_scripts/4556838
submission_time:            03/01/2022 17:32:33.131
owner:                      zgao1
uid:                        583371011
group:                      zgao1
gid:                        583371011
supplementary group:        compbio, ngs, camhpcusers, zgao1
sge_o_home:                 /home/zgao1
sge_o_log_name:             zgao1
sge_o_path:                 /usr/lib64/qt-3.3/bin:/camhpc/admin/sge/8.6.6/bin/lx-amd64:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/home/zgao1/bin
sge_o_shell:                /bin/bash
sge_o_workdir:              /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x
sge_o_host:                 camhpcps01
account:                    sge
hard resource_list:         h_rt=979200,h_vmem=192G
mail_list:                  zgao1@camhpcps01.hpc.biogen.com
notify:                     FALSE
job_name:                   out.majiq
stdout_path_list:           NONE:NONE:/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-1000nM//logs/
priority:                   0
jobshare:                   0
hard_queue_list:            all.q
env_list:                   
script_file:                STDIN
parallel environment:       thread range: 12
department:                 defaultdepartment
binding:                    NONE
mbind:                      NONE
submit_cmd:                 qsub -o /home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_NGN2-BIO-2060884-1000nM//logs/ -N out.majiq -q all.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 12
category_id:                1261
request_dispatch_info:      FALSE
start_time            1:    03/01/2022 17:32:33.193
job_state             1:    r
exec_host_list        1:    camhpcpe01:12
usage                 1:    wallclock=03:56:52, cpu=17:11:02, mem=142461.48881 GBs, io=2375.46641 GB, iow=448.680 s, ioops=168493217, vmem=6.478G, maxvmem=25.987G
scheduling info:            -

#5) show avaiable queue # https://stackoverflow.com/questions/7152884/how-do-i-find-a-complete-list-of-available-torque-pbs-queues
qhost -q
HOSTNAME                ARCH         NCPU NSOC NCOR NTHR NLOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
----------------------------------------------------------------------------------------------
global                  -               -    -    -    -     -       -       -       -       -
# camhpcpe01              lx-amd64       56    2   28   56  0.05  504.5G    5.9G   16.0G    2.3M
#    long.q               BIP   0/0/1         
#    all.q                BP    0/12/56       
# camhpcpe10              lx-amd64       28    2   28   28  0.08  504.7G   11.2G   16.0G  341.7M
#    medium.q             BIP   0/0/28        
#    all.q                BP    0/12/28       
# camhpcpe12              lx-amd64       64    4   32   64  0.01  378.5G    4.9G   16.0G  229.9M
#    short.q              BIP   0/0/64        
#    all.q                BP    0/0/64        

 bash optparse_qsub_make_master_table_ZG.sh 
Your job 4556868 ("out.master_TST11872_iPSC-BIO-1755497-12nM") has been submitted
~/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x zgao1@camhpcps01 $ 
$  qstat # http://docs.adaptivecomputing.com/torque/4-1-3/Content/topics/commands/qstat.htm
job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
------------------------------------------------------------------------------------------------------------------------------------------------
   4556838 0.57333 out.majiq  zgao1        r     03/01/2022 17:32:33 all.q@camhpcpe01                                                 12        
   4556839 0.57333 out.majiq  zgao1        r     03/01/2022 17:33:01 all.q@camhpcpe26.hpc.biogen.co                                   12        
   4556859 0.57333 out.majiq  zgao1        r     03/01/2022 20:12:49 all.q@camhpcpe32.hpc.biogen.co                                   12        
   4556860 0.57333 out.majiq  zgao1        r     03/01/2022 20:15:50 all.q@camhpcpe31.hpc.biogen.co                                   12        
   4556861 0.57333 out.majiq  zgao1        r     03/01/2022 20:22:12 all.q@camhpcpe22.hpc.biogen.co                                   12        
   4556862 0.57333 out.majiq  zgao1        r     03/01/2022 20:28:27 all.q@camhpcpe10.hpc.biogen.co                                   12        
   4556719 0.57333 out.majiq  zgao1        Eqw   02/28/2022 19:51:18                                                                  12        
   4556726 0.57333 out.majiq  zgao1        Eqw   02/28/2022 19:58:23                                                                  12        
   4556727 0.57333 out.majiq  zgao1        Eqw   02/28/2022 19:58:38                                                                  12        
   4556728 0.57333 out.majiq  zgao1        Eqw   02/28/2022 19:58:53                                                                  12        
   4556729 0.57333 out.majiq  zgao1        Eqw   02/28/2022 19:59:29                                                                  12        
   4556730 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:06:08                                                                  12        
   4556731 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:06:48                                                                  12        
   4556732 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:07:10                                                                  12        
   4556733 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:07:28                                                                  12        
   4556734 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:07:44                                                                  12        
   4556735 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:07:59                                                                  12        
   4556736 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:08:16                                                                  12        
   4556737 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:08:31                                                                  12        
   4556738 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:08:47                                                                  12        
   4556739 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:09:07                                                                  12        
   4556740 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:09:22                                                                  12        
   4556741 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:09:36                                                                  12        
   4556742 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:09:51                                                                  12        
   4556743 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:11:50                                                                  12        
   4556744 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:12:07                                                                  12        
   4556745 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:12:25                                                                  12        
   4556746 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:12:41                                                                  12        
   4556747 0.57333 out.rmats  zgao1        Eqw   02/28/2022 20:13:02                                                                  12        
   4556709 0.50000 out.leafcu zgao1        Eqw   02/28/2022 19:33:04                                                                   1        
   4556710 0.50000 out.leafcu zgao1        Eqw   02/28/2022 19:33:19                                                                   1        
   4556711 0.50000 out.leafcu zgao1        Eqw   02/28/2022 19:33:35                                                                   1        
   4556712 0.50000 out.leafcu zgao1        Eqw   02/28/2022 19:33:50                                                                   1        
~/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x zgao1@camhpcps01 $ 
$  qstat -j 4556868
Following jobs do not exist or permissions are not sufficient: 
4556868

# https://stackoverflow.com/questions/28857807/use-qdel-to-delete-all-my-jobs-at-once-not-one-at-a-time
qdel -u zgao1 
zgao1 has registered the job 335860 for deletion
zgao1 has registered the job 335869 for deletion
zgao1 has registered the job 335870 for deletion

2022-03-02
qstat
job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
------------------------------------------------------------------------------------------------------------------------------------------------
    336705 0.57833 out.majiq  zgao1        r     03/03/2022 18:50:06 all.q@ginseng07.hpc.biogen.com                                   12        
    336707 0.57833 out.majiq  zgao1        r     03/03/2022 18:52:03 all.q@ginseng05.hpc.biogen.com                                   12        
    336708 0.57833 out.majiq  zgao1        r     03/03/2022 19:06:49 all.q@yarrow07.hpc.biogen.com   


# install old packages
# https://cran.r-project.org/web/packages/rowr/index.html
# https://cran.r-project.org/src/contrib/Archive/rowr/
require(devtools)
install_version("rowr", version = "1.1.3", repos = "http://cran.us.r-project.org")
library(rowr)

# biomaRt.R
## ----setup, echo = FALSE--------------------------------------------------------------------------
# knitr::opts_chunk$set(error = TRUE, cache = FALSE, eval = TRUE, out.width = "100%")
# httr::set_config(httr::config(ssl_verifypeer = FALSE))
# options(width=100)

## ----useEnsembl-----------------------------------------------------------------------------------
library(biomaRt)
listEnsembl()
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl

## ----listDatasets---------------------------------------------------------------------------------
datasets <- listDatasets(ensembl)
head(datasets)

mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listAttributes( mart) # ====<<<< # ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# rank                                                      Exon rank in transcript    structure

genes = c("PDE7A",  "POMT2",  "EVC",  "DENND1A",  "WDR27",  "KCNT2",  "DLGAP4",   "DHFR",  "ASAP1",  "ZFP82",  "KLHL1",  "LIMCH1",  "STRN3",  "TAF2",  "ABCB8",  "MADD",  "KDM6A",  "VPS29",
          "MELTF",  "WDR17",  "ERGIC3",  "TRIM65",  "CLASP2",  "AXIN1",  "STXBP5L",  "RAF1",  "EXOC6B",  "ARL15",  "HDX",  "LARP7",  "LRRC42",  "SKP1",  "PPIP5K1",  "LINC02210",  "PHOX2B-AS1", 
          "VDAC3",  "AL391807.1",  "AL589743.1",  "AC136604.3",  "PITPN",  "XRN2",  "PDXDC1",  "FHOD3",  "KCNT2",  "CIP2A",  "TENT2",  "HLTF",  "HTR3A",  "GNAQ",  "KDSR",  "L3MBTL2",  "SLC7A6",
          "ADAMTS19",  "RDX",  "DLG5",  "BTBD10",  "FBL",  "ASNS",  "AGPS",  "MLLT10",  "RAD21",  "ELP4",  "VPS41")

getAtt <- c('external_gene_name','chromosome_name', 'start_position', 'end_position',  'strand','exon_chrom_start','exon_chrom_end', 'ensembl_exon_id', 'rank')

elocs=getBM(attributes=getAtt,filters="hgnc_symbol",value=genes,mart=mart)
print(elocs) 

elocs #rownames(elocs)<-NULL # elocs_order = elocs [ with(elocs, order(external_gene_name, rank) ),  ] 

elocs_order = elocs [ with(elocs, order(external_gene_name, exon_chrom_start,  exon_chrom_end, ensembl_exon_id, rank) ),  ]  # unique_elocs_order = unique_elocs [ with(unique_elocs, order(external_gene_name, rank) ),  ] 
elocs_order 


elocs_order[!duplicated(elocs_order), ]
unique_elocs_order  <- elocs_order[!duplicated(elocs_order), ]
unique_elocs_order 


# R seletion
coi_ar=dir(din,pattern="_.*x") #coi_ar= dir(din) coi_ar=dir(din,pattern="TST11955") # this one remove DMSO, contain files only have 3x and 10x , # I decicde to remove this one to keep  "iPSC_DMSO_bridge"   and in later analysis to remove DMSO_bride
coi_ar1=grep("excrep4",coi_ar,value=T) #round4 excluded 
coi_ar2=coi_ar[!coi_ar %in% coi_ar1] # all 5. 
# https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
coi_ar = coi_ar[grep('iPSC_1755497_10x|iPSC_1755497_3x',  coi_ar)] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
coi_ar = coi_ar[grep('1755497',  coi_ar)] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
coi_ar = coi_ar[grep('bridge',invert = TRUE,  coi_ar)] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
coi_ar = grep( coi_ar , pattern = "1755497") # code in leafviz7.R  orkding_dirs = grep( workding_dirs , pattern = 'backup-',invert = TRUE, value = TRUE)   # exclude "backupXXX" folder not needed


# search
ll -R | grep  biogen_splicing_pipeline_mry_alltools
-rwxrwxr-x. 1 zgao1 ngs 23347 Jun 20 01:17 sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py
-rwxrwxr-x. 1 zgao1 ngs 23322 May 27 00:47 biogen_splicing_pipeline_mry_alltools_v3.py
-rwxrwx---. 1 zgao1 ngs 24150 May 30 02:06 biogen_splicing_pipeline_mry_alltools_v3_ZG.py
-rwxrwx---. 1 zgao1 ngs 24283 Jun 18 02:33 biogen_splicing_pipeline_mry_alltools_v3_ZG_MAJIQ.py
-rwxrwx---. 1 zgao1 ngs 24150 Jun 18 01:24 biogen_splicing_pipeline_mry_alltools_v3_ZG.py
-rwxrwxr-x. 1 zgao1 ngs 23347 Mar  1 15:12 sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py
(base) /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code zgao1@camhpcps01 $ 
$  find -name sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py
./splicing/sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py
./example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/2022-03-02-ZG_verified/Explore/sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py

find -name sub.biogen_splicing_pipeline*
./splicing/sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py
/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/2022-03-02-ZG_verified/Explore/sub.biogen_splicing_pipeline_mry_alltools_v3_DHmod.py

(base) /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code zgao1@camhpcps01 $ 
$  find -name biogen_splicing_pipeline_mry_alltools_v3_ZG_MAJIQ.py
/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/MAJIQ_Viola_code_explore/biogen_splicing_pipeline_mry_alltools_v3_ZG_MAJIQ.py

# tool_untilities.sh
# rename 001 '' *.gz 
 rename TST11872_ 'rmat_TST11872_' TST* 
 
 
 #(1)summarized sizes of directories and their subdirectories?
du -sh /*
# What this means:
# 
# -s to give only the total for each command line argument.
# -h for human-readable suffixes like M for megabytes and G for gigabytes (optional).
# /* simply expands to all directories (and files) in /.
# 
# Note: dotfiles are not included; run shopt -s dotglob to include those too.
# 
# Also useful is sorting by size:
du -sh /* | sort -h

# (2)
 tree -L 1 TST11872_*
# TST11872_iPSC-BIO-1755497-12nM
# ├── leafcutter
# ├── logs
# └── TST11872_iPSC-BIO-1755497-12nM
# TST11872_iPSC-BIO-1755497-40nM
# ├── leafcutter
# ├── logs
# └── TST11872_iPSC-BIO-1755497-40nM

# http://www.hypexr.org/linux_scp_help.php 
#scp -r foo your_username@remotehost.edu:/some/remote/directory/bar
# leafcutter_files = file.path(opt$in_dir, "leafcutter",  c("leafcutter_ds_res_cluster_significance.txt","leafcutter_ds_res_effect_sizes.txt","leafviz.Rdata"))
#scp -r foo/ your_username@remotehost.edu:/some/remote/directory/bar
cd /edgehpc/dept/compbio/users/zgao1/ZG_AS_code/Analysis_2_Edge/TST11872_iPSC-BIO-2060573-3uM/leafcutter
scp leafcutter_ds_res_cluster_significance.txt leafcutter_ds_res_effect_sizes.txt leafviz.Rdata zgao1@camhpcps01.hpc.biogen.com:/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM-Edge-1st/
# Warning: Permanently added 'camhpcps01.hpc.biogen.com,10.9.76.162' (RSA) to the list of known hosts.
# zgao1@camhpcps01.hpc.biogen.com's password: 
# tput: No value for $TERM and no -T specified
# leafcutter_ds_res_cluster_significance.txt                                                                                 100% 1154KB   7.7MB/s   00:00    
# leafcutter_ds_res_effect_sizes.txt                                                                                         100% 5039KB  37.3MB/s   00:00    
# leafviz.Rdata                               

cd /edgehpc/dept/compbio/users/zgao1/ZG_AS_code/Analysis_2_Edge/TST11872_iPSC-BIO-2060573-3uM_2nd_for_consistance/leafcutter
scp leafcutter_ds_res_cluster_significance.txt leafcutter_ds_res_effect_sizes.txt leafviz.Rdata zgao1@camhpcps01.hpc.biogen.com:/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM-Edge-2nd/
# zgao1@camhpcps01.hpc.biogen.com's password: 
# tput: No value for $TERM and no -T specified
# leafcutter_ds_res_cluster_significance.txt                                                                                 100% 1154KB   7.1MB/s   00:00    
# leafcutter_ds_res_effect_sizes.txt                                                                                         100% 5039KB  37.6MB/s   00:00    
# leafviz.Rdata                                                                                                              100%   15MB  98.0MB/s   00:00 

# COPY  CAMHPC folder to Edge
cd ~/Working_code
scp -r * zgao1@edge.hpc.biogen.com:/home/zgao1/CompBio_zgao1_LK/Working_R_code/ # copy all files and subfolders
#or
scp -r * zgao1@edge.hpc.biogen.com:/edgehpc/dept/compbio/users/zgao1/Working_R_code/

2022-08-13

cd ~/NGS_projects
scp -r JoySplicing/ zgao1@edge.hpc.biogen.com:/camhpc/home/zgao1/NGS_projects/ # copy all files and subfolders
#or
scp -r JoySplicing/ zgao1@edge.hpc.biogen.com:/home/camhpc_zgao1/NGS_projects/ # copy all files and subfolders

# ==== HPC address
# Dann Huh
# Tue 1/11/2022 4:11 PM
# 
# CAMHPC: cambridge.hpc.biogen.com
# Yarrow: camhpcve23.hpc.biogen.com
# Edge: edge.hpc.biogen.com

scp -r /home/jpoulo/pacbio/TST11702/ TST11702 zgao1@edge.hpc.biogen.com:/edgehpc/dept/compbio/users/zgao1/PacBio

/edgehpc/dept/compbio/users/zgao1 
#./CompBio_zgao1_LK -> /edgehpc/dept/compbio/users/zgao1
/edgehpc/dept/compbio/users/zgao1/PacBio

Pacbio Raw_data
/edgehpc/dept/instruments/PacBio/SequelII/raw_data/TST11702/r64138_20210714_191607/1_A01
-rw-r--r-- 1 svc-bings svc-bings 700287738915 Jul 16  2021 m64138_210714_193827.subreads.bam
-rw-r--r-- 1 svc-bings svc-bings   1480500261 Jul 16  2021 m64138_210714_193827.subreads.bam.pbi

/edgehpc/dept/instruments/PacBio/SequelII/raw_data/TST11702/r64138_20210714_191607/2_B01

/edgehpc/dept/instruments/PacBio/SequelII/smrtlink/userdata/jobs_root/cromwell-executions/
drwxrwsr-x 1 svc-bings svc-bings 0 Oct 21  2020 pb_ccs
drwxrwsr-x 1 svc-bings svc-bings 0 Sep 10  2021 pb_ccs_demux
drwxrwsr-x 1 svc-bings svc-bings 0 May 28  2021 pb_ccs_demux_auto
drwxrwsr-x 1 svc-bings svc-bings 0 May 14  2021 pb_demux_subreads_auto
drwxrwsr-x 1 svc-bings svc-bings 0 Sep 10  2021 sl_import_subreads
drwxrwsr-x 1 svc-bings svc-bings 0 Oct 21  2020 sl_merge_datasets

/edgehpc/dept/compbio/users/jpoulo/jpoulo
# CAMHPC /home/jpoulo/pacbio/TST11702/

    Welcome to Edge HPC cluster
    Please find documentation about the system here:
    https://wiki.biogen.com/display/HPC/EDGE+Documentation

    OS Image stage: "prod"
    OS Image build: "2022-04-28"
    Repo snapshot : "20220321"

# grep 'chr4:3212709'     /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/Result/*/*.csv  >  HTT_eoi_start_only_loose.csv #ZG 2022-08-24:2min.  #ZG 2022-04-07:5min
#/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/RG7916_10x
cp /camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/*/*.csv  "/Users/zgao1/OneDrive - Biogen/"
cd "/Users/zgao1/OneDrive - Biogen "
-bash: cd: /Users/zgao1/OneDrive - Biogen : No such file or directory
~/OneDrive - Biogen :  
zgao1@ Zhens-MacBook-Pro $ cd "/Users/zgao1/OneDrive - Biogen"

## tool_R
# (1) remove empty
install.packages("janitor") 
library(janitor) 
remove_empty(x) # remove empty line in dataframe, not in file


# (2) Selecting data frame rows based on partial string match in a column
#https://stackoverflow.com/questions/13043928/selecting-data-frame-rows-based-on-partial-string-match-in-a-column 
library(data.table)
mtcars[rownames(mtcars) %like% "Merc", ]
iris[iris$Species %like% "osa", ]

# dplyr
# select(iris,contains("Sepal"))
# See the Selection section in ?select for numerous other helpers like starts_with, ends_with, etc.

DSG_3v3_2nd_6152 = select( DSG_3v3_2nd , contains("2006152")  )



QC_raw_table_data_set1 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "iPSC"))
QC_raw_table_data_set2 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "NGN2")) 

# grep also ok
mtcars[grep("Merc", rownames(mtcars)), ]


# (3) select a column
QC_raw_table_data      %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set1 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set2 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()

Mapped_reads_lower_threshold  = quantile( ( QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ), 0.25)  * 0.5 

QC_raw_table_data[      which( (      QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set2[ which( ( QC_raw_table_data_set2$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 

QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , c( "Total.Purity.Filtered.Reads.Sequenced")] 


# (4) - order the rows of data frame to display
QC_raw_table_data[ order(rownames(QC_raw_table_data)), c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_Total.passed_Reads_percent", "rRNA", "rRNA.rate","Alternative.Aligments") ] # 


QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(Sample, "4088")) 
QC_raw_table_data_4088[ order(QC_raw_table_data_4088$Sample), c("Sample", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs", "rRNA", "rRNA.rate","Alternative.Aligments")] # QC_raw_table_data_set1[ order(QC_raw_table_data_set1$Sample), ]


# (5) change all column of Data fraome to numeric
# https://r-lang.com/how-to-convert-column-to-numeric-in-r/#:~:text=To%20convert%20all%20the%20columns,the%20columns%20were%20a%20factor.
QC_raw_table_data[] <- lapply(QC_raw_table_data, function(x) as.numeric(as.character(x)))
QC_raw_table_data
sapply(QC_raw_table_data, class)


# (6) select  and order 
# (6.1)-4088. #==>Chimeric.Pairs ? BAD
QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(rownames(QC_raw_table_data), "4088")) 
QC_raw_table_data_4088[ order(rownames(QC_raw_table_data_4088)) , c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_passed_Reads_Pct", "rRNA", "rRNA.rate","Alternative.Aligments")] # 
# (6.2)
QC_raw_table_data  %>%  select ( order( colnames( QC_raw_table_data) ) )  %>% filter(str_detect(rownames(QC_raw_table_data), "4088")) %>% select("Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced")  


# (7) change column names 
colnames(dfall) = gsub("PROFESSION", "", colnames(dfall))

tool 
# (1) remove empty
install.packages("janitor") 
library(janitor) 
remove_empty(x) # remove empty line in dataframe, not in file


# (2) Selecting data frame rows based on partial string match in a column
#https://stackoverflow.com/questions/13043928/selecting-data-frame-rows-based-on-partial-string-match-in-a-column 
library(data.table)
mtcars[rownames(mtcars) %like% "Merc", ]
iris[iris$Species %like% "osa", ]

# dplyr
# select(iris,contains("Sepal"))
# See the Selection section in ?select for numerous other helpers like starts_with, ends_with, etc.

DSG_3v3_2nd_6152 = select( DSG_3v3_2nd , contains("2006152")  )



QC_raw_table_data_set1 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "iPSC"))
QC_raw_table_data_set2 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "NGN2")) 

# grep also ok
mtcars[grep("Merc", rownames(mtcars)), ]


# (3) select a column
QC_raw_table_data      %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set1 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set2 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()

Mapped_reads_lower_threshold  = quantile( ( QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ), 0.25)  * 0.5 

QC_raw_table_data[      which( (      QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set2[ which( ( QC_raw_table_data_set2$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 

QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , c( "Total.Purity.Filtered.Reads.Sequenced")] 


# (4) - order the rows of data frame to display
QC_raw_table_data[ order(rownames(QC_raw_table_data)), c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_Total.passed_Reads_percent", "rRNA", "rRNA.rate","Alternative.Aligments") ] # 


QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(Sample, "4088")) 
QC_raw_table_data_4088[ order(QC_raw_table_data_4088$Sample), c("Sample", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs", "rRNA", "rRNA.rate","Alternative.Aligments")] # QC_raw_table_data_set1[ order(QC_raw_table_data_set1$Sample), ]


# (5) change all column of Data fraome to numeric
# https://r-lang.com/how-to-convert-column-to-numeric-in-r/#:~:text=To%20convert%20all%20the%20columns,the%20columns%20were%20a%20factor.
QC_raw_table_data[] <- lapply(QC_raw_table_data, function(x) as.numeric(as.character(x)))
QC_raw_table_data
sapply(QC_raw_table_data, class)


# (6) select  and order 
# (6.1)-4088. #==>Chimeric.Pairs ? BAD
QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(rownames(QC_raw_table_data), "4088")) 
QC_raw_table_data_4088[ order(rownames(QC_raw_table_data_4088)) , c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_passed_Reads_Pct", "rRNA", "rRNA.rate","Alternative.Aligments")] # 
# (6.2)
QC_raw_table_data  %>%  select ( order( colnames( QC_raw_table_data) ) )  %>% filter(str_detect(rownames(QC_raw_table_data), "4088")) %>% select("Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced")  


# (7) change column names 
colnames(dfall) = gsub("PROFESSION", "", colnames(dfall))

# (1) remove empty
install.packages("janitor") 
library(janitor) 
remove_empty(x) # remove empty line in dataframe, not in file


# (2) Selecting data frame rows based on partial string match in a column
#https://stackoverflow.com/questions/13043928/selecting-data-frame-rows-based-on-partial-string-match-in-a-column 
library(data.table)
mtcars[rownames(mtcars) %like% "Merc", ]
iris[iris$Species %like% "osa", ]

# dplyr
# select(iris,contains("Sepal"))
# See the Selection section in ?select for numerous other helpers like starts_with, ends_with, etc.

DSG_3v3_2nd_6152 = select( DSG_3v3_2nd , contains("2006152")  )



QC_raw_table_data_set1 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "iPSC"))
QC_raw_table_data_set2 = QC_raw_table_data  %>%  filter(str_detect( rownames(QC_raw_table_data), "NGN2")) 

# grep also ok
mtcars[grep("Merc", rownames(mtcars)), ]


# (3) select a column
QC_raw_table_data      %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set1 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()
QC_raw_table_data_set2 %>% select("Total.Purity.Filtered.Reads.Sequenced") %>% summary()

Mapped_reads_lower_threshold  = quantile( ( QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ), 0.25)  * 0.5 

QC_raw_table_data[      which( (      QC_raw_table_data$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 
QC_raw_table_data_set2[ which( ( QC_raw_table_data_set2$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , ] 

QC_raw_table_data_set1[ which( ( QC_raw_table_data_set1$Total.Purity.Filtered.Reads.Sequenced ) <   Mapped_reads_lower_threshold ) , c( "Total.Purity.Filtered.Reads.Sequenced")] 


# (4) - order the rows of data frame to display
QC_raw_table_data[ order(rownames(QC_raw_table_data)), c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_Total.passed_Reads_percent", "rRNA", "rRNA.rate","Alternative.Aligments") ] # 


QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(Sample, "4088")) 
QC_raw_table_data_4088[ order(QC_raw_table_data_4088$Sample), c("Sample", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs", "rRNA", "rRNA.rate","Alternative.Aligments")] # QC_raw_table_data_set1[ order(QC_raw_table_data_set1$Sample), ]


# (5) change all column of Data fraome to numeric
# https://r-lang.com/how-to-convert-column-to-numeric-in-r/#:~:text=To%20convert%20all%20the%20columns,the%20columns%20were%20a%20factor.
QC_raw_table_data[] <- lapply(QC_raw_table_data, function(x) as.numeric(as.character(x)))
QC_raw_table_data
sapply(QC_raw_table_data, class)


# (6) select  and order 
# (6.1)-4088. #==>Chimeric.Pairs ? BAD
QC_raw_table_data_4088 = QC_raw_table_data  %>%  filter(str_detect(rownames(QC_raw_table_data), "4088")) 
QC_raw_table_data_4088[ order(rownames(QC_raw_table_data_4088)) , c( "Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced", "Chimeric.Pairs_over_passed_Reads_Pct", "rRNA", "rRNA.rate","Alternative.Aligments")] # 
# (6.2)
QC_raw_table_data  %>%  select ( order( colnames( QC_raw_table_data) ) )  %>% filter(str_detect(rownames(QC_raw_table_data), "4088")) %>% select("Chimeric.Pairs", "Total.Purity.Filtered.Reads.Sequenced")  


# (7) change column names 
colnames(dfall) = gsub("PROFESSION", "", colnames(dfall))

