# install.packages("misc3d") #error: Tcl/Tk support is not available on this system
# #if (!require("BiocManager", quietly = TRUE))     install.packages("BiocManager") # BiocManager::install("edgeR")
# tmp=TST11872_DSG[, c(1), drop = F ] #  TST11872_DSG[, 1, drop = F ]   # data[is.na(data$ColWtCL_6),] # https://stackoverflow.com/questions/17013227/select-only-rows-if-its-value-in-a-particular-column-is-na-in-r
# head(tmp)
# (tmp)
# na.omit(tmp)
#i=1

require("devtools")
library(stringr)
require(limma)
require(multtest)
require(survival)
require(plyr)
library(dplyr)
require(scales)
require("edgeR")
require(DESeq2)
require(RColorBrewer)
require(grid)
library(gplots)
require(ggplot2)
library(ggrepel)
library(plot3D)
require("pheatmap")



TST11872_DSG = read.delim2("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST11872/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" , header = T, sep = "," )
Out_dir      = paste0("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST11872/TST11872_individual_compound")
#head(TST11872_DSG )
# dim(TST11872_DSG )
df_dPSI      = TST11872_DSG 

TST11955_DSG = read.delim2("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST11955/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" , header = T, sep = "," )
Out_dir      = paste0("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST11955/TST11955_individual_compound")
df_dPSI      = TST11955_DSG 

TST12086_DSG = read.delim2("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST12086/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" , header = T, sep = "," )
Out_dir      = paste0("/Users/zgao1/OneDrive - Biogen/Meetings/Project_11872_11955_12086_summary/TST12086/TST12086_individual_compound")
df_dPSI      = TST12086_DSG 

########
dir.create(Out_dir, showWarnings = F, recursive = T)
Out_dir
setwd(Out_dir)


row.names(df_dPSI ) = df_dPSI$X
df_dPSI$X = NULL

head(df_dPSI)

df_dPSI = df_dPSI[ , -grep("minpadj_", colnames(df_dPSI))]
colnames(df_dPSI) =  gsub(".BIO.", "_", colnames(df_dPSI))

df_dPSI
dim(df_dPSI)
head(df_dPSI)

#i
for(i in 1:ncol(df_dPSI)){
  temp_df = na.omit(df_dPSI[, i , drop = F ])
  temp_df[, 1] = abs( as.numeric(temp_df[, 1]) )
  
  temp_df =  temp_df[ order( temp_df[, 1],  decreasing = TRUE) , , drop=FALSE] # https://stackoverflow.com/questions/31433362/use-order-on-a-single-column-data-frame
  
  column = colnames(temp_df)[1]
  
  write.csv(temp_df,  paste0(column, ".csv"), sep = ",", row.names = TRUE) # NOT use table, miss the 1st black cell
}



##############################################=======================================================================================================================================##############################################
library(tidyverse)
require("gdata") # #  package ‘list_c’ is not available for this version of R # res <- list.cbind(ldf)


#####================== 11872
TST11872_DSG_3x_DIR = "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST11872/TST11872_individual_compound/3x" 

#https://stackoverflow.com/questions/9564489/read-all-files-in-a-folder-and-apply-a-function-to-each-data-frame
#merge all data frames together # https://www.statology.org/merge-multiple-data-frames-in-r/         list.cbind(ldf) #res <- Reduce( function(x, y) cbind(x, y), ldf) # Reduce(function(x, y) merge(x, y, all=TRUE), df_list)   
filenames <- list.files(TST11872_DSG_3x_DIR, pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
TST11872_DSG_summary = ldf %>% reduce(cbindX) 
TST11872_DSG_summary
dim(TST11872_DSG_summary)
colnames(TST11872_DSG_summary)
# change 1000uM to 3x 
colnames(TST11872_DSG_summary) =  c("x" , "dPSI_iPSC_1755497_3x" , "x" , "dPSI_iPSC_1949634_3x" , "x" , "dPSI_iPSC_2006152_3x" , 
                                    "x" , "dPSI_iPSC_2059811_3x" , "x" , "dPSI_iPSC_2060573_3x" , "x" , "dPSI_iPSC_2060884_3x" , 
                                    "x" , "dPSI_iPSC_2070692_3x" , "x" , "dPSI_iPSC_2135644_3x" , "x" , "dPSI_iPSC_2136770_3x" , 
                                    "x" , "dPSI_NGN2_1755497_3x" , "x" , "dPSI_NGN2_1949634_3x" , "x" , "dPSI_NGN2_2006152_3x" , 
                                    "x" , "dPSI_NGN2_2059811_3xM", "x" , "dPSI_NGN2_2060573_3x" , "x" , "dPSI_NGN2_2060884_3x" , 
                                    "x" , "dPSI_NGN2_2070692_3x" , "x" , "dPSI_NGN2_2135644_3x" , "x" , "dPSI_NGN2_2136770_3x" )
colnames(TST11872_DSG_summary) 
head(TST11872_DSG_summary) 

write.csv(TST11872_DSG_summary,  "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST11872/TST11872_DSG_summary.csv", row.names = TRUE) 


#####================== 11955
TST11955_DSG_3x_DIR = "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST11955/TST11955_individual_compound/3x" 

#https://stackoverflow.com/questions/9564489/read-all-files-in-a-folder-and-apply-a-function-to-each-data-frame
#merge all data frames together # https://www.statology.org/merge-multiple-data-frames-in-r/         list.cbind(ldf) #res <- Reduce( function(x, y) cbind(x, y), ldf) # Reduce(function(x, y) merge(x, y, all=TRUE), df_list)   
filenames <- list.files(TST11955_DSG_3x_DIR, pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
TST11955_DSG_summary = ldf %>% reduce(cbindX) 
TST11955_DSG_summary
dim(TST11955_DSG_summary)
colnames(TST11955_DSG_summary)
# change 4088 to 2184088
colnames(TST11955_DSG_summary) =  c( "X" ,    "dPSI_iPSC_2184088_3X" , "X" ,    "dPSI_iPSC_2184090_3X" , "X" ,    "dPSI_iPSC_2174714_3X" , "X" ,    "dPSI_iPSC_2184741_3X" , "X" ,    "dPSI_iPSC_2174748_3X" , "X" ,    "dPSI_iPSC_2175420_3X" ,
                                     "X" ,    "dPSI_iPSC_2006152_3X" , "X" ,    "dPSI_iPSC_2186527_3X" , "X" ,    "dPSI_iPSC_2186866_3X" , "X" ,    "dPSI_iPSC_2186960_3X" , "X" ,    "dPSI_iPSC_2139701_3X" , 
                                     "X" ,    "dPSI_NGN2_2184088_3X" , "X" ,    "dPSI_NGN2_2184090_3X" , "X" ,    "dPSI_NGN2_2174714_3X" , "X" ,    "dPSI_NGN2_2184741_3X" , "X" ,    "dPSI_NGN2_2174748_3X" , "X" ,    "dPSI_NGN2_2175420_3X" , 
                                     "X" ,    "dPSI_NGN2_2006152_3X" , "X" ,    "dPSI_NGN2_2186527_3X" , "X" ,    "dPSI_NGN2_2186866_3X" , "X" ,    "dPSI_NGN2_2186960_3X" , "X" ,    "dPSI_NGN2_2139701_3X" )
colnames(TST11955_DSG_summary) 
grep("2006152" , colnames(TST11955_DSG_summary) )
grep("2139701" , colnames(TST11955_DSG_summary) )
grep("2174714" , colnames(TST11955_DSG_summary) )
grep("2174748" , colnames(TST11955_DSG_summary) )
grep("2175420" , colnames(TST11955_DSG_summary) )
grep("2184088" , colnames(TST11955_DSG_summary) )
grep("2184090" , colnames(TST11955_DSG_summary) )
grep("2184741" , colnames(TST11955_DSG_summary) )
grep("2186527" , colnames(TST11955_DSG_summary) )
grep("2186866" , colnames(TST11955_DSG_summary) )
grep("2186960" , colnames(TST11955_DSG_summary) )

TST11955_DSG_summary = TST11955_DSG_summary[ ,c(13, 14,  21, 22, 5, 6, 9, 10 ,11, 12 , 1, 2 , 3, 4  ,7, 8 ,15, 16, 17, 18, 19, 20,        35, 36, 43, 44, 27, 28 ,31, 32, 33, 34, 23, 24, 25, 26, 29, 30, 37, 38 , 39,  40, 41, 42)] 
head(TST11955_DSG_summary )
write.csv(TST11955_DSG_summary,  "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST11955/TST11955_DSG_summary.csv", row.names = TRUE) 

#####================== 12086
TST12086_DSG_3x_DIR = "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST12086/TST12086_individual_compound/3X"  

#https://stackoverflow.com/questions/9564489/read-all-files-in-a-folder-and-apply-a-function-to-each-data-frame
#merge all data frames together # https://www.statology.org/merge-multiple-data-frames-in-r/         list.cbind(ldf) #res <- Reduce( function(x, y) cbind(x, y), ldf) # Reduce(function(x, y) merge(x, y, all=TRUE), df_list)   
filenames <- list.files(TST12086_DSG_3x_DIR, pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
TST12086_DSG_summary = ldf %>% reduce(cbindX) 
TST12086_DSG_summary
dim(TST12086_DSG_summary)
colnames(TST12086_DSG_summary)

colnames(TST12086_DSG_summary) 

write.csv(TST12086_DSG_summary,  "/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Meetings/TST12086/plots/Magnus_requested/2023-05-25-requested/TST12086/TST12086_DSG_summary.csv", row.names = TRUE) 


# #record_to_answer_request-magnus-2023-05-25-based-on-2022-11-07

1) TST11872
/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv # reproduced 2023-march-27

#### NOT filtered /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/ZG_analysis.03.splice_analysis_1/summarytable_dPSI_passed_padj0.05_PdPSI0.9_M50.csv

2) TST11955
/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold_3v3_to_show_10-13/backup-good/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv 


### XXX===== > NOT correct ===>> /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All-code_without_R4/analysis.02.splicing_offtargets/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv

3) TST11872
/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results/analysis.03.splice_analysis/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv





Code from:

/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_1/run.04.DSG_counts_plot_ranking_batch123_DSG_count_Magnus_Dann_final_5.R

# ##############################   111111----------------------------------------- iPSC Note ======PTC518 compound is BIO-2197294
# ##############################   1-Batch 1 ############################## ############################## ##############################  
din1 = "/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.03.splice_analysis/"
fin1 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
# param = gsub(".csv","",gsub("binary_","",fin1))
dPSI_bi = read.csv( paste0(din1, fin1)  , stringsAsFactors = F , row.names = 1 , check.names = F) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp            =  dPSI_bi[     , grep('iPSC-',  colnames(dPSI_bi) ) ] 
head(dPSI_bi_tmp)
colnames(dPSI_bi_tmp)  = gsub("BIO-", "", colnames(dPSI_bi_tmp ) )
colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c("1755497_iPSC_3x" , "1755497_iPSC_10x" , "1949634_iPSC_3x" , "1949634_iPSC_10x" , "2006152_iPSC_3x" , "2006152_iPSC_10x" , 
                           "2059811_iPSC_3x" , "2059811_iPSC_10x" , "2060573_iPSC_3x" , "2060573_iPSC_10x" , "2060884_iPSC_3x" , "2060884_iPSC_10x" , 
                           "2070692_iPSC_3x" , "2070692_iPSC_10x" , "2135644_iPSC_3x" , "2135644_iPSC_10x" , "2136770_iPSC_3x" , "2136770_iPSC_10x" ) 

# 1755497 is Branaplam, 1949634 Risdiplam, 206088s Interal Lead
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('_3x'  ,  colnames( dPSI_bi_tmp) )] #
head(dPSI_bi_tmp)
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('1949634',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('2060884',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
head(dPSI_bi_tmp)

dPSI_bi = dPSI_bi_tmp
y_1st        = colSums(dPSI_bi) ##################### get the totall DS events
y_1st        = y_1st[order(names(y_1st) )]
#names(y_1st) = gsub("dPSI", "", names(y_1st))
y_1st # class(y_1st) #y_1st_rank = rank(y_1st ) # y_rank = rank(-y )

# ##############################   1-Batch 2 ############################## ############################## ##############################  
din2 = "/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold_3v3_to_show_10-13/backup-good/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin2 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
param = gsub(".csv","",gsub("binary_","",fin1))
dPSI_bi = read.csv( paste0(din2, fin2)  , stringsAsFactors = F , row.names = 1 , check.names = F) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp            =  dPSI_bi[     , grep('iPSC-',  colnames(dPSI_bi) ) ] 
dPSI_bi_tmp            =  dPSI_bi_tmp[ , grep('-3x'  ,  colnames( dPSI_bi_tmp) )] #
head(dPSI_bi_tmp) # colnames(dPSI_bi_tmp)  = gsub("BIO-", "", colnames(dPSI_bi_tmp ) )
colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c("2184088_iPSC_3x", "2184090_iPSC_3x", "2174714_iPSC_3x", "2184741_iPSC_3x", "2174748_iPSC_3x", "2175420_iPSC_3x", 
                           "2006152_iPSC_3x", "2186527_iPSC_3x", "2176866_iPSC_3x", "2186960_iPSC_3x", "2139701_iPSC_3x" ) 
head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)

dPSI_bi = dPSI_bi_tmp
y_2nd        = colSums(dPSI_bi) ##################### get the totall DS events
y_2nd        = y_2nd[order(names(y_2nd) )]
#names(y_2nd) = gsub("dPSI", "", names(y_2nd))
y_2nd # class(y_2nd) #y_2nd_rank = rank(y_2nd ) # y_rank = rank(-y ) #y_2nd_rank


# ##############################   1-Batch 3 ############################## ############################## ############################## 
din3 = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_1/analysis.03.splice_analysis/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin3 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#param = gsub(".csv","",gsub("binary_","",fin3))
