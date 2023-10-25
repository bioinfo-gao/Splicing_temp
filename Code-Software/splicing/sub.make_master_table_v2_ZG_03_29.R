# this is test leafcutter and majiq in camhpc and edgehpc for discrependance 2022-03-28
.libPaths(c("/home/lhou1/R/x86_64-pc-linux-gnu-library/3.6","/camhpc/pkg/R/3.6.1/centos7/lib64/R/library" , .libPaths()))
install.packages(optparse)

library(optparse)
library(tidyverse)
library(data.table)
library(reshape2)
library(glue)

# source /etc/profile.d/modules_bash.sh; 
# module load xz/5.2.2; 
# module load hdf5/1.8.17; 
# module load R/3.6.1; 

# # prog_dir="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/splicing/"
# Rscript $prog_dir/sub.make_master_table_v2.R 
# --in_dir $in_dir 
# --out_dir $out_dir 
# --out_prefix $out_prefix 
# --majiq_cutoffs $majiq_cutoffs
# opt = list(in_dir       = "/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-12nM/", 
#            out_dir      = "/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-12nM/", 
#            out_prefix   = "TST11872_iPSC-BIO-1755497-12nM",
#            majiq_cutoffs = "0.1 0.2 0.3 0.4")
opt = list(in_dir       = "/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/", 
           out_dir      = "/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/", 
           out_prefix   = "TST11872_iPSC-BIO-2060573-3uM",
           majiq_cutoffs = "0.1 0.2 0.3 0.4")

#| qsub -N out.master_${out_prefix} -q long.q -l h_rt=272:00:00 -l h_vmem=192G -pe thread 1 -V # cpu.q

# option_list = list(
#   make_option("--in_dir", action = "store", default = NA, type = "character",
#               help = "Full path to working directory containing majiq, rmats, leafcutter results"),
#   make_option("--out_dir", action = "store", default = NA, type = "character",
#               help = "Full path to output directory."),
#   make_option("--out_prefix", action = "store", default = NA, type = "character",
#               help = "The name of the output table.  File output will look like {out_prefix}_master_table.csv"),
#   make_option("--majiq_cutoffs", action="store", default=NA, type = "character",
#               help = "The cutoffs used for majiq_v1 results.  Separate by space, comma, or semicolon.")
# )
# opt = parse_args(OptionParser(option_list=option_list))
##Debug
# if(FALSE){
#   opt = list(in_dir = "/home/mryals/splicing/HNDS_120nM_res/", 
#              out_dir = "/home/mryals/splicing/HNDS_120nM_res/", 
#              out_prefix = "HNDS_120nM_res",
#              majiq_cutoffs = "0.1 0.2 0.3 0.4")
# }

#Convert majiq cutoffs to a vector
majiq_cutoffs = unlist(strsplit(opt$majiq_cutoffs, split = "[ ,;]+"))

#Check that the folders are in the in_dir
need.dirs <- c("majiq_v1","leafcutter")# need.dirs <- c("majiq_v1","leafcutter","rmats") 
test.dirs <- file.path(opt$in_dir, need.dirs)

if(all(dir.exists(test.dirs))){
  message(paste0("all analysis directories are found in ", opt$in_dir))
} else{ stop(paste0("analysis directories are missing from ", opt$in_dir)) }

#Check that the files are in the in_dir
need.files <- c(
  # file.path(opt$in_dir, "rmats", 
  #           c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt')),
  file.path(opt$in_dir, "leafcutter", 
            c("leafcutter_ds_res_cluster_significance.txt","leafcutter_ds_res_effect_sizes.txt","leafviz.Rdata")) ####,  # leafCutter_failed________!!!!!!!!!!!!!!! 2022-02-22
  # file.path(opt$in_dir, "majiq_v1", paste0("cutoff_",majiq_cutoffs,"_reference_alternative.deltapsi.tsv") )
)
############================================
# 
# file.exists(need.files)
# [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE

if(all(file.exists(need.files))){
  message(paste0("all analysis files are found in ", opt$in_dir))
} else{ stop(paste0("analysis files are missing from ", opt$in_dir)) }

################## Prep rMATS portion of master table ###############

#rMATS helper functions

#make M
mcol <- function(df){
  #ID column is duplicated, force unique
  #names(df) <- make.unique(names(df))
  
  #Map function
  map.m <- function(x) {sum(as.numeric(unlist(strsplit(x, split = ",")))) / length(as.numeric(unlist(strsplit(x, split = ","))))}
  
  #map.m map for the collapsed columns
  df.i <- df %>% 
    select(ID, IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2) %>%
    mutate(i1 = map(IJC_SAMPLE_1, function(x) map.m(x)),
           i2 = map(IJC_SAMPLE_2, function(x) map.m(x)),
           s1 = map(SJC_SAMPLE_1, function(x) map.m(x)),
           s2 = map(SJC_SAMPLE_2, function(x) map.m(x))
    ) %>%
    select(ID, i1, i2, s1, s2) %>%
    unnest(c(i1,i2,s1,s2)) %>%
    mutate(M = pmax(i1,i2,s1,s2)) %>% 
    select(ID, M)
}

#make psi_1 / psi_2 averages
psicol <- function(df){
  #Map function
  map.m <- function(x) {mean(as.numeric(unlist(strsplit(x, split = ","))),na.rm = TRUE)}
  
  df.i <- df %>% 
    select(ID, PSI_1_persub, PSI_2_persub) %>%
    mutate(psi1 = map(PSI_1_persub, function(x) map.m(x)),
           psi2 = map(PSI_2_persub, function(x) map.m(x))
           ) %>%
    select(ID, psi1, psi2) %>%
    unnest(c(psi1, psi2)) %>%
    select(ID, PSI_1 = psi1, PSI_2 = psi2)
}

#Names of the relevant rmats files within the working directory
rmats_files = file.path(opt$in_dir, "rmats", 
                        c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
rmats_names = c('A3SS', 'A5SS', 'MXE', 'RI', 'SE')

#Read in the files
rmats_list <- lapply(rmats_files, fread)
rmats_format <- list()
#Format the files
for(i in 1:length(rmats_list)){
  
  file0 <- as.data.frame(rmats_list[[i]])
  
  #ID column is duplicated, force unique
  names(file0) <- make.unique(names(file0))
  
  #Basic reformatting
  file1 <- file0 %>% 
    mutate(Event = rmats_names[i],
           Algorithm = "rmats") %>%
    rename("dPSI" = "IncLevelDifference",
           "PSI_1_persub" = "IncLevel1",
           "PSI_2_persub" = "IncLevel2") 
  #create m
  file.m <- mcol(file0)
  #Add M
  file2 <- file1 %>% 
    left_join(file.m)
  #create PSI_1 and PSI_2 (averages)
  file.psi <- psicol(file2)
  #Add PSI_1 and PSI_2
  file3 <- file2 %>% 
    left_join(file.psi)
  #Event info (collapse all columns not in the keep.col)
  keep.col <- c('chr','strand','geneSymbol', 'GeneID', 'Event', 'Algorithm', 'FDR',
                'PValue', 'dPSI', 'M', 'PSI_1', 'PSI_2')
  file.event <- file3 %>% select(-one_of(keep.col)) %>% select(-ID.1) %>% select(ID, everything())
  file.event[] <- Map(paste, names(file.event), file.event, sep=":")
  file.event <- file.event %>% 
    mutate(ID = as.integer( gsub("ID:","",ID)) ) %>%
    unite(col = Event_info, colnames(.)[-1],sep="|") %>%
    select(ID, Event_info)
  #Join back to file2
  rmats_format[[i]] <- file3 %>% left_join(file.event) %>% select(keep.col, Event_info)
}
#Make the rMATS portion of the table
rmats_format_df <- rbindlist(rmats_format)

############################# Prep leafcutter portion of master table ##################################
#Names of the relevant leafcutter files within the working directory
leafcutter_files = file.path(opt$in_dir, "leafcutter", 
                             c("leafcutter_ds_res_cluster_significance.txt","leafcutter_ds_res_effect_sizes.txt","leafviz.Rdata"))

# read in files:
#pvalues
table_cs <- fread(leafcutter_files[1])
#drop NA from table_cs (can keep track of it with the status)
table_cs.badstatus <- table_cs %>% filter(status != "Success")
table_cs <- table_cs %>% drop_na
#psis
table_es <- fread(leafcutter_files[2])
#anno table
load(leafcutter_files[3])
table_anno <- clusters 

#add a cluster column that matches table_cs to table_es (format: chr:clu_number)
table_es$cluster <- paste0(str_split_fixed(table_es$intron, ":", 5)[,1],":",str_split_fixed(table_es$intron, ":", 5)[,4])

#make the cluster column in table_anno match as well
table_anno$cluster <- paste0(str_split_fixed(table_anno$coord,":",2)[,1],":",table_anno$clusterID)

#drop unneeded columns from table_anno, format gene
table_anno <- table_anno %>% 
  mutate(gene = gsub("(<i>|</i>)","",gene) ) %>%
  select(-clusterID, -FDR)

#merge table_es and table_cs
table_es_cs_anno <- table_es %>% full_join(table_cs, by="cluster") %>% full_join(table_anno, by="cluster")

#format the columns
leafcutter_format_df.all <- table_es_cs_anno %>%
  select(geneSymbol = gene, 
         Event = cluster,
         FDR = p.adjust,
         PValue = p,
         dPSI = deltapsi,
         M = status,
         PSI_1 = alternative,
         PSI_2 = reference,
         Annotation_leafcutter = annotation,
         everything()) %>%
  mutate(Algorithm = "leafcutter")

leafcutter_format_df.all$chr <- str_split_fixed(leafcutter_format_df.all$coord, ":", 2)[,1]
leafcutter_format_df.all$strand <- str_split_fixed(leafcutter_format_df.all$intron, ":", 5)[,5] #I don't think these are always right

leafcutter_format_df.all <- leafcutter_format_df.all %>% filter(FDR < 1)

#make the event info (use intron as an index for joining)
keep.col <- c("chr","strand","geneSymbol", "Event", "FDR", "PValue", "dPSI", "M", "PSI_1", "PSI_2", "Annotation_leafcutter", "Algorithm")
leafcutter.format.event <- leafcutter_format_df.all %>% select(-one_of(keep.col)) %>% select(intron, everything())
leafcutter.format.event[] <- Map(paste, names(leafcutter.format.event), leafcutter.format.event, sep=":")
leafcutter.format.event <- leafcutter.format.event %>% 
  mutate(intron = gsub("intron:","",intron),
         intron.info = intron) %>%
  unite(col = Event_info, colnames(.)[-1],sep="|") %>%
  select(intron, Event_info)

#make the final leafcutter file
leafcutter_format_df.final <- leafcutter_format_df.all %>% select(intron, one_of(keep.col)) %>% full_join(leafcutter.format.event, by="intron") %>%
  select(-intron)

dim(leafcutter_format_df.final )
head(leafcutter_format_df.final ,20)

write.table(leafcutter_format_df.final,  "/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/leafcutter_output_ZG.tsv", row.names=T, sep = "\t", quote=F) ###### 2022-03-29


############################# Prep majiq (v1) portion of master table #######################

#Get the files
#Names of the relevant majiq files within the working directory
majiq_files = file.path(opt$in_dir, "majiq_v1", paste0("cutoff_",majiq_cutoffs,"_reference_alternative.deltapsi.tsv"))
#Names of each cutoff (flexible)
majiq_cutoff_names <- str_split_fixed(basename(majiq_files),"_",4)[,2]

#Read in the files
majiq_list <- lapply(majiq_files, fread)
majiq_format <- list()

for(i in 1:length(majiq_list)){
  
  majiq0 <- majiq_list[[i]]
  
  #Start to format, several columns need to be un-collapsed
  majiq.format.all <- majiq0 %>% 
    select(geneSymbol = `#Gene Name`,
           GeneID = `Gene ID`,
           Event = `LSV ID`, #Use as index
           everything()
    ) %>%
    mutate(Algorithm = "majiq",
           Event_type = case_when(A5SS & !A3SS & !ES ~ "A5SS",
                                  A3SS & !A5SS & !ES ~ "A3SS",
                                  ES & !A3SS & !A5SS ~ "ES",
                                  A5SS & A3SS & !ES ~ "A5SS_A3SS",
                                  A5SS & ES & !A3SS ~ "A5SS_ES",
                                  A3SS & ES & !A5SS ~ "A3SS_ES",
                                  A3SS & A5SS & ES ~ "A5SS_A3SS_ES",
                                  !A5SS & !A3SS & !ES ~ "MXE/RI") ) #look at any of these cases, anything with >4 exons too
  
  #These columns need to be un-collapsed
  majiq.expand <- majiq.format.all %>% select(Event, 
                                              dPSI = `E(dPSI) per LSV junction`,
                                              PSI_1 = `alternative E(PSI)`,
                                              PSI_2 = `reference E(PSI)`,
                                              Junc_coords = `Junctions coords`, 
                                              PdPSI = matches("^P\\(.*>") #regex flexible for different cutoffs
                                              ) 
  
  expand.cols <- c("dPSI", "PSI_1", "PSI_2", 
                   "Junc_coords", 
                   "PdPSI")
  expand.list <- list()
  for(j in 1:length(expand.cols)){
    expand.list[[j]] <- majiq.expand %>% select(Event, !!expand.cols[[j]]) %>% separate_rows(!!expand.cols[[j]], sep=";")
  }
  
  #Make the uncollapsed dataframe and remove Event column duplicates
  majiq.expand <- do.call(cbind, expand.list)
  majiq.expand <- majiq.expand[!duplicated(as.list(majiq.expand))]
  
  #Join expand to the table
  majiq.format.df1 <- majiq.format.all %>% dplyr::select(chr, strand, geneSymbol, GeneID, Event, Event_type, Algorithm) %>%
    full_join(majiq.expand)
  
  #Some columns to keep
  keep.col = c("chr","strand",'geneSymbol', "GeneID", 
               'Event_type', 'Algorithm', 
               'PdPSI', 'dPSI', 'PSI_1', 'PSI_2', 'Junc_coords', 'Event_info')
  
  #Make the event info (Event column as index)
  majiq.format.event <- majiq.format.all %>% dplyr::select(Event, `LSV Type`, `Exons coords`, `IR coords`)
  
  majiq.format.event[] <- Map(paste, names(majiq.format.event), majiq.format.event, sep=":")
  majiq.format.event <- majiq.format.event %>% 
    mutate(Event = gsub("Event:","",Event)) %>%
    unite(col = Event_info, colnames(.)[-1],sep="|") %>%
    select(Event, Event_info)
  
  majiq_format[[i]] <- majiq.format.df1 %>% full_join(majiq.format.event) 
}
names(majiq_format) <- majiq_cutoff_names
majiq_format_df.final <- rbindlist(majiq_format, idcol = "Cutoff_dPSI")

############ Merge all tables and write out ##################
master_table <- rbindlist(list(rmats_format_df, leafcutter_format_df.final, majiq_format_df.final), fill = TRUE )
fwrite(master_table, file.path(opt$out_dir, glue("{opt$out_prefix}_master_table.csv")), row.names = FALSE)
