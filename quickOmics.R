# RNA_Seq_raw2quickomics_12188_5-plot3_christina.R
# # #================>>>>>>>>>  check "/edgehpc/dept/compbio/users/zgao1/project_RNAseq/TST12188/Differential_RNAseq_12086_3vs3_based_quickOmics1.Rmd" 
#rstudioapi::writeRStudioPreference("console_max_lines", 300)
require(tidyverse)
require(reshape2)
require(DESeq2)
require(dplyr)
require(UpSetR)
require(ggrepel)

rm(list=ls())

# ### (1)4vs4 GM
# DEG_full_out_dir =              "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/GM_Full_DEG_Result_reanalysis_09_29_4vs4/" # _09_28_4vs4
# DEG_gene_out_dir =              paste0(DEG_full_out_dir, "DEG")
# load(                           "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/EA20230829_0/TST12188_GM_Ris.RData")
# sample_sheet     = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/samplesheet.tsv",  header = T, sep ="\t") #"
# comparison_sheet = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/EA20230829_0/data/compareInfo.csv",  header = T, sep =",") #"
# TST12188.tpm     = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/genes.tpm.tsv",  header = T, sep ="\t") #"genes.tpm_table.txt"
# TST12188.count   = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/genes.expected_count.tsv",  header = T, sep ="\t") #"
# 
### (2)4vs4 SH
# DEG_full_out_dir      =              "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/Sh_Full_DEG_Result_reanalysis_09_29_4vs4/"
# DEG_gene_out_dir      =              paste0(DEG_full_out_dir, "DEG")
# load(                                "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230828_0/TST12188_Sh_Ris.RData") #load( "12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230828_0/TST12188_ShSy5Y.RData")
# sample_sheet          = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/samplesheet.tsv",  header = T, sep ="\t") #"
# comparison_sheet      = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230828_0/data/compareInfo.csv",  header = T, sep =",") #"
# TST12188.tpm          = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/genes.tpm.tsv",  header = T, sep ="\t") #"genes.tpm_table.txt"
# TST12188.count        = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/genes.expected_count.tsv",  header = T, sep ="\t") #"
# sample_sheet$Index_ID = sample_sheet$SampleID
# 
# ### (3)3vs3 GM
DEG_full_out_dir =              "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/GM_Full_DEG_Result_reanalysis_09_29_3vs3/"
DEG_gene_out_dir =              paste0(DEG_full_out_dir, "DEG_plots")
load(                           "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/EA20230929_0/TST12188_GM09677c_3vs3.RData")
sample_sheet     = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/samplesheet.tsv",  header = T, sep ="\t") #"
comparison_sheet = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/EA20230829_0/data/compareInfo.csv",  header = T, sep =",") #"
TST12188.tpm     = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/genes.tpm.tsv",  header = T, sep ="\t") #"genes.tpm_table.txt"
TST12188.count   = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/genes.expected_count.tsv",  header = T, sep ="\t") #"
# 
### (4)3vs3 SH
# DEG_full_out_dir      =              "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/Sh_Full_DEG_Result_reanalysis_09_29_3vs3/"
# DEG_gene_out_dir      =              paste0(DEG_full_out_dir, "DEG")
# load(                                "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230929_0/TST12188_Sh_GR38_3vs3.RData") #load( "12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230828_0/TST12188_ShSy5Y.RData")
# sample_sheet          = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/samplesheet.tsv",  header = T, sep ="\t") #"
# comparison_sheet      = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/EA20230828_0/data/compareInfo.csv",  header = T, sep =",") #"
# TST12188.tpm          = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/genes.tpm.tsv",  header = T, sep ="\t") #"genes.tpm_table.txt"
# TST12188.count        = read.delim2( "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/genes.expected_count.tsv",  header = T, sep ="\t") #"
# sample_sheet$Index_ID = sample_sheet$SampleID

dir.create(DEG_gene_out_dir , recursive = T)
setwd(DEG_gene_out_dir) # DEG_gene_out_dir setwd(DEG_full_out_dir) 
getwd()                          

ls() 

dim(MetaData)         # 92  44 for GM === 68  44 for ShSy5Y 
dim(comp_info)        # 11*2 (2 dose vs DMSO) + 10*2 (2 dose vs Risdplam) =  42 comparison for GM cells  ================ NULL for  /TST12188.RData
dim(data_long)        # 3579996   for GM < =======  >2681602 rows for ShSy5Y
dim(data_results)     # 38913 X 60 for GM : UniqueID Gene.Name  id Intensity ====> SH_BIO_2196772_3x_Mean  && SH_BIO_2196772_3x_sd ** SH_BIO_1949634_3x_Mean  $$ SH_BIO_1949634_3x_sd  ===> SH_BIO_2196772_3x_vs_SH_DMSO_DESeq2.log2FoldChange ** SH_BIO_2196772_3x_vs_SH_DMSO_DESeq2.pvalue 
dim(data_wide)        # Sample-093 Sample-094 Sample-095  ....... Sample-159 Sample-160
dim(ProteinGeneName)  # 60760 proteins
dim(results_long)   

head(data_long)    # 
head(data_wide)    #  Sample-093 Sample-094 Sample-095  ....... Sample-159 Sample-160
head(MetaData)     #  92 samples for GM, 68  44 for ShSy5Y 

#############
comp_info          # 
head(comp_info)    # 

ProteinGeneName   
head(ProteinGeneName) 
tail(ProteinGeneName) 

data_results
data_results[1:6, 1:12]
dim(data_results) 
# head(data_results %>% select( matches("GM_BIO_2196772_3x")  ) ) # https://www.statology.org/r-select-columns-containing-string/

head(data_results) #  UniqueID Gene.Name  id Intensity ====> SH_BIO_2196772_3x_Mean  && SH_BIO_2196772_3x_sd ** SH_BIO_1949634_3x_Mean  $$ SH_BIO_1949634_3x_sd ===> SH_BIO_2196772_3x_vs_SH_DMSO_DESeq2.log2FoldChange ** SH_BIO_2196772_3x_vs_SH_DMSO_DESeq2.pvalue 

results_long   
head(results_long)   
tail(results_long)   
dim(results_long)   

########## 
sample_sheet
head(sample_sheet)
tail(sample_sheet)

comparison_sheet
head(comparison_sheet)
#data_results_summary  = data_results 

########## 
name_vector = unique(results_long$test)
name_vector
results_long_list  = list()
results_short_list = list()
DEG_with_gn        = list()

########## 

colnames(data_results)

for(i in 1:length(name_vector)){    #for(i in 1:21){   #for(i in 22:42){   
    #     i=1
    print(i)
    print( name_vector[[i]])
    
    Compare_name = as.character(name_vector[[i]]) 
    
    results_long_list[[i]] = results_long[ results_long$test == name_vector[i], ]   

    results_short_list[[i]] = results_long_list[[i]][, c("UniqueID" , "test" , "logFC" , "P.Value" , "Adj.P.Value")]
    head(results_short_list[[i]]) 
    
    results_short_list[[i]] =  results_short_list[[i]][ which( results_short_list[[i]]$Adj.P.Value < 0.05 & abs(results_short_list[[i]]$logFC) > 0.585 ),  c("UniqueID", "logFC"), drop = FALSE  ] #     DEG_with_count_TPM  =  DEG_with_count_TPM[ sort(abs(DEG_with_count_TPM$logFC), decreasing=T) ,  ]

    colnames(results_short_list[[i]]) = c("UniqueID",    Compare_name)
        
    dim(results_short_list[[i]]) 
    head(results_short_list[[i]]) 
}

(results_short_list )
num_DEG = Reduce(function(x, y) merge(x, y, all=T), results_short_list )
num_DEG

num_DEG_with_gn = merge(  ProteinGeneName[ , c( "UniqueID", "Gene.Name" ) ] , num_DEG , sort = F, all.y = T) 
num_DEG_with_gn
dim(num_DEG_with_gn)
num_DEG_with_gn[1:6, 1:6]

row.names(num_DEG_with_gn) = num_DEG_with_gn$Gene.Name
num_DEG_with_gn$Gene.Name = NULL
num_DEG_with_gn$UniqueID = NULL

getwd()

fout = paste0( DEG_gene_out_dir,  "/", "GM","_DEG_3vs3", ".csv")


fout
 
write.csv( num_DEG_with_gn, file=fout, quote=F, row.names=T)
num_DEG_with_gn_SMN_HTT = num_DEG_with_gn[which( row.names( num_DEG_with_gn) %in% c("SMN2", "SMN1",  "HTT") ), ]
num_DEG_with_gn_SMN_HTT 
# 
fout = paste0( DEG_gene_out_dir,  "/", "SMN_HTT_DEG", ".csv")
write.csv( num_DEG_with_gn_SMN_HTT, file=fout, quote=F, row.names=F)

####################################################################################################################################################
dim(num_DEG_with_gn)  # 340  42
head(num_DEG_with_gn)
num_DEG_with_gn[1:6, 1:10]

dDEG_bi = num_DEG_with_gn

row.names(dDEG_bi ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL
colnames(dDEG_bi ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL

dDEG_bi[ is.na(dDEG_bi)       ] = 0
dDEG_bi[   abs(dDEG_bi) >0    ] = 1
dDEG_bi

colnames(dDEG_bi) 
#colnames(dDEG_bi) = gsub("dPSI_" , "", colnames(dDEG_bi) ) # colnames(dDEG_bi) = gsub("dPSI_","", gsub("_max.abs.dPSI..","",colnames(dDEG_bi)))
colnames(dDEG_bi) = gsub("BIO_" , "", colnames(dDEG_bi) )  # y %>% select(sub('_ln$', '', filter_vector))

head(dDEG_bi)
dim(dDEG_bi)

# getwd()
# fout = paste0("4vs4_binary_","GM",".csv") 
# fout = paste0("4vs4_binary_","SH",".csv") 
# fout = paste0("3vs3_binary_","GM",".csv") 
# fout = paste0("3vs3_binary_","SH",".csv") 
# 
# write.csv(dDEG_bi,file = fout,quote = F)


#####################################################################################################################################################
getwd()

####==========================================####
head(dDEG_bi)
dDEG_bi[1:6, 1:6]
dim(dDEG_bi)
colnames(dDEG_bi)
colnames(dDEG_bi) = gsub( "_3x" , "_03x", colnames(dDEG_bi) )  # y %>% select(sub('_ln$', '', filter_vector))
unique( sapply(strsplit( colnames(dDEG_bi),  split="_\\d{1,2}x_vs_") , `[`, 1) )  #numDSG_DMSO$cellline_dose = sapply(strsplit(numDSG_DMSO$cellline_DMSO, "_vs_") , `[`, 1) 
# [Risdiplam] "GM_BIO_1949634"
# [1] "GM_2196772" "GM_2207180" "GM_2178782" "GM_2196895" "GM_2197306" "GM_2206678" "GM_2201042" "GM_2204984" "GM_2199562" "GM_2186827" "GM_1949634"
# [1] "SH_2196772" "SH_2207180" "SH_2178782" "SH_2196895" "SH_2197306" "SH_2204984" "SH_2199562" "SH_1949634"

dDEG_bi_DMSO = dDEG_bi[ ,  grepl( "vs_.._DMSO"      , colnames(dDEG_bi)  ) ]  # grepl NOT grep dDEG_bi_Ris  = dDEG_bi[ ,  grepl( "vs_.._Risdiplam" , colnames(dDEG_bi)  ) ]  # grepl NOT grep dDEG_bi_Ris  = dDEG_bi[ ,  grepl( "vs_.._1949634" , colnames(dDEG_bi)  ) ]  # grepl NOT grep
dDEG_bi_DMSO

# # 1
dDEG_bi_working = dDEG_bi_DMSO[ ,  grepl( "_10x_vs_"     , colnames(dDEG_bi_DMSO)  ) ]  # grepl NOT grep
# # 2
#dDEG_bi_working = dDEG_bi_DMSO[ ,  grepl( "_03x_vs_"     , colnames(dDEG_bi_DMSO)  ) ]  # grepl NOT grep 

# # 3 dDEG_bi_working = dDEG_bi_Ris[ ,  grepl( "_10x_vs_"     , colnames(dDEG_bi_Ris)  ) ]  # grepl NOT grep # # # 4 #dDEG_bi_working = dDEG_bi_Ris[ ,  grepl( "_03x_vs_"     , colnames(dDEG_bi_Ris)  ) ]  # grepl NOT grep

###### 
dDEG_bi_working = dDEG_bi_working[!(rowSums(dDEG_bi_working) == 0), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na
dDEG_bi_working
dDEG_bi_working[1:2,]
dim(dDEG_bi_working)
colSums(dDEG_bi_working)

colnames_set = colnames(dDEG_bi_working) # 15 types
colnames_set

# fout = paste0("1_UpSet_All_compound_10x_DMSO", ".png")
fout = paste0("2_UpSet_All_compound_03x_DMSO", ".png")  # fout = paste0("3_UpSet_All_compound_10x_vs_Ris", ".png") #fout = paste0("4_UpSet_All_compound_03x_vs_Ris", ".png")

getwd()
png(fout, height = 1000, width = 680+ncol(dDEG_bi_working)*50)            #png(fout,height = 500, width = length(y)*50)
upset(dDEG_bi_working , sets = colnames_set,  mainbar.y.label = "DEG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()

####==========================================####
########### ==== 
# WT_CHX (WT_100 and WT_500) --vs-- DMSO_CHX
###################################################

####===================== GM =====================####
head(dDEG_bi_DMSO)
dim(dDEG_bi_DMSO)
colnames(dDEG_bi_DMSO)
unique( sapply(strsplit( colnames(dDEG_bi_DMSO),  split="_\\d{1,2}x_vs_") , `[`, 1) )

dDEG_bi_working = dDEG_bi_DMSO[ ,  grepl( "^GM_1949634"    , colnames(dDEG_bi_DMSO)  ) ]  # grepl NOT grep
fout = paste0("5_UpSet_Ris", ".png")

dDEG_bi_working
dDEG_bi_working = dDEG_bi_working[!(rowSums(dDEG_bi_working) == 0), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na
head(dDEG_bi_working)
dim(dDEG_bi_working)

getwd()
colnames_set = colnames(dDEG_bi_working) # 15 types
colnames_set

png(fout, height = 1000, width = 680+ncol(dDEG_bi_working)*50)            #png(fout,height = 500, width = length(y)*50)
upset(dDEG_bi_working , sets = colnames_set,  mainbar.y.label = "DEG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()

# # 
####==========================================####
### 10X GM 
####==========================================####
dDEG_bi_DMSO_working = dDEG_bi_DMSO[ ,  grep( "_10x_vs_GM"    , colnames(dDEG_bi_DMSO)  ) ]  # grepl NOT grep  # grepl( "^SH_2186827"    , colnames(dPSI_bi)  )
dDEG_bi_DMSO_working                ##                                                      dPSI_bi_working =  dPSI_bi_DMSO_working[,  grep( "^GM_2178782_..x_vs_|^GM_1949634_..x_vs_"    , colnames( dPSI_bi_DMSO_working)  ) ]  # grepl NOT grep

# dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2178782_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2186827_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2196772_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2196895_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2197306_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
# # 
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2199562_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2201042_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2204984_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2206678_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
  dDEG_bi_working = dDEG_bi_DMSO_working[ ,  grepl( "^GM_2207180_..x_vs_|^GM_1949634_..x_vs_"    , colnames(dDEG_bi_DMSO_working)  ) ]  # grepl NOT grep
# 
# fout = paste0("6A_UpSet_BIO_2178782_vs_Ris_10X", ".png")
  fout = paste0("7A_UpSet_BIO_2186827_vs_Ris_10X", ".png")
  fout = paste0("8A_UpSet_BIO_2196772_vs_Ris_10X", ".png")
  fout = paste0("9A_UpSet_BIO_2196895_vs_Ris_10X", ".png")
  fout = paste0("10A_UpSet_BIO_2197306_vs_Ris_10X", ".png")
  fout = paste0("11A_UpSet_BIO_2199562_vs_Ris_10X", ".png")
  fout = paste0("12A_UpSet_BIO_2201042_vs_Ris_10X", ".png")
  fout = paste0("13A_UpSet_BIO_2204984_vs_Ris_10X", ".png")
  fout = paste0("14A_UpSet_BIO_2206678_vs_Ris_10X", ".png")
  fout = paste0("15A_UpSet_BIO_207180_vs_Ris_10X", ".png")
# 
dDEG_bi_working
dDEG_bi_working = dDEG_bi_working[!(rowSums(dDEG_bi_working) == 0), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na
head(dDEG_bi_working) # dDEG_bi_working[1:2,]
dim(dDEG_bi_working)

getwd()
colnames_set = colnames(dDEG_bi_working) # 15 types
colnames_set
png(fout, height = 1000, width = 680+ncol(dDEG_bi_working)*50)            #png(fout,height = 500, width = length(y)*50)
upset(dDEG_bi_working , sets = colnames_set,  mainbar.y.label = "DEG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()







############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments 
# https://stackoverflow.com/questions/21537782/how-to-set-fixed-continuous-colour-values-in-ggplot2
library("RColorBrewer")
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(rev(brewer.pal(7, "Dark2"))) 
sc <- scale_colour_gradientn(colours = myPalette(100) ) # sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-4, 4))  # ggplot(subset(mtcars, am==0), aes(x=wt, y=mpg, colour=carb)) +      geom_point(size=6) + sc
############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments ############################ 10x experiments 

dDEG = num_DEG_with_gn

row.names(dDEG ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL
colnames(dDEG ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL


dDEG[ ,  grepl( "_10x_vs_GM_DMSO"    , colnames(dDEG)  ) ]  # grepl NOT grep  #dDEG_Ris  = dDEG[ ,  grepl( "vs_.._1949634" , colnames(dDEG)  ) ]  # grepl NOT grep

Compounds_10x_vs_GM_DMSO = dDEG[ ,  grepl( "_10x_vs_GM_DMSO"    , colnames(dDEG)  ) ]  # grepl NOT grep  #dDEG_Ris  = dDEG[ ,  grepl( "vs_.._1949634" , colnames(dDEG)  ) ]  # grepl NOT grep
dim(Compounds_10x_vs_GM_DMSO)
head(Compounds_10x_vs_GM_DMSO)
Compounds_10x_vs_GM_DMSO[1:2, ]

# Compound1          2178782
dPSI_working = BIO_2178782_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2178782_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2178782_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

dim(dPSI_working)
head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2178782_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2178782_10x_vs_GM_DMSO)] <- 0

dPSI_working

getwd()
fout = "Compound1_BIO_2178782_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2178782_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO, color = GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()


# Compound2       2186827
dPSI_working = BIO_2186827_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2186827_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2186827_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2186827_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2186827_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound2_BIO_2186827_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2186827_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,  color = GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound3        2196772
dPSI_working = BIO_2196772_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2196772_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2196772_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2196772_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2196772_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound3_BIO_2196772_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2196772_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO, color = GM_BIO_1949634_10x_vs_GM_DMSO  ) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound4 2196895
dPSI_working = BIO_2196895_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2196895_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] 
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2196895_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2196895_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound4_BIO_2196895_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2196895_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO, color =  GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound5 2197306
dPSI_working = BIO_2197306_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2197306_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2197306_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2197306_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2197306_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound5_BIO_2197306_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2197306_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,   color = GM_BIO_1949634_10x_vs_GM_DMSO ) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound6 2199562
dPSI_working = BIO_2199562_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2199562_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2199562_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2199562_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2199562_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound6_BIO_2199562_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2199562_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,  color =  GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound7 2201042
dPSI_working = BIO_2201042_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2201042_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2201042_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2201042_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2201042_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound7_BIO_2201042_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2201042_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,   color =  GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound8    2204984
dPSI_working = BIO_2204984_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2204984_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2204984_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2204984_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2204984_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound8_BIO_2204984_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2204984_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,   color = GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound9  2206678
dPSI_working = BIO_2206678_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2206678_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2206678_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2206678_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2206678_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound9_BIO_2206678_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2206678_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO,   color = GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound10 2207180
dPSI_working = BIO_2207180_risdiplam_columns =  Compounds_10x_vs_GM_DMSO[ , grep( "2207180_..x_vs_|1949634_..x_vs_"    , colnames(Compounds_10x_vs_GM_DMSO)  ) ] # head(BIO_2207180_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_10x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2207180_10x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2207180_10x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound10_BIO_2207180_10x_VS_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2207180_10x_vs_GM_DMSO, y= GM_BIO_1949634_10x_vs_GM_DMSO, color = GM_BIO_1949634_10x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()


################# 3x   =================================================================================================================================================================================================================================================================================

dDEG = num_DEG_with_gn

row.names(dDEG ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL
colnames(dDEG ) # row.names(dDEG_bi ) =  make.names(dDEG_bi$genesymbol, unique=TRUE) # dDEG_bi $genesymbol = NULL


dDEG[ ,  grepl( "_3x_vs_GM_DMSO"    , colnames(dDEG)  ) ]  # grepl NOT grep  #dDEG_Ris  = dDEG[ ,  grepl( "vs_.._1949634" , colnames(dDEG)  ) ]  # grepl NOT grep

Compounds_3x_vs_GM_DMSO = dDEG[ ,  grepl( "_3x_vs_GM_DMSO"    , colnames(dDEG)  ) ]  # grepl NOT grep  #dDEG_Ris  = dDEG[ ,  grepl( "vs_.._1949634" , colnames(dDEG)  ) ]  # grepl NOT grep
dim(Compounds_3x_vs_GM_DMSO)
head(Compounds_3x_vs_GM_DMSO)
Compounds_3x_vs_GM_DMSO[1:2, ]

# Compound1          2178782
dPSI_working = BIO_2178782_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2178782_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2178782_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

dim(dPSI_working)
head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO
dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2178782_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2178782_3x_vs_GM_DMSO)] <- 0

dPSI_working

getwd()
fout = "Compound1_BIO_2178782_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2178782_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO, color = GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()


# Compound2       2186827
dPSI_working = BIO_2186827_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2186827_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2186827_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2186827_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2186827_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound2_BIO_2186827_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2186827_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,  color = GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound3        2196772
dPSI_working = BIO_2196772_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2196772_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2196772_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2196772_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2196772_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound3_BIO_2196772_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2196772_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO, color = GM_BIO_1949634_3x_vs_GM_DMSO  ) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound4 2196895
dPSI_working = BIO_2196895_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2196895_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] 
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2196895_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2196895_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound4_BIO_2196895_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2196895_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO, color =  GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound5 2197306
dPSI_working = BIO_2197306_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2197306_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2197306_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2197306_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2197306_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound5_BIO_2197306_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2197306_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,   color = GM_BIO_1949634_3x_vs_GM_DMSO ) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound6 2199562
dPSI_working = BIO_2199562_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2199562_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2199562_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2199562_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2199562_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound6_BIO_2199562_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2199562_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,  color =  GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound7 2201042
dPSI_working = BIO_2201042_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2201042_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2201042_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2201042_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2201042_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound7_BIO_2201042_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2201042_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,   color =  GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound8    2204984
dPSI_working = BIO_2204984_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2204984_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2204984_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2204984_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2204984_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound8_BIO_2204984_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2204984_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,   color = GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound9  2206678
dPSI_working = BIO_2206678_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2206678_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2206678_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2206678_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2206678_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound9_BIO_2206678_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2206678_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO,   color = GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Compound10 2207180
dPSI_working = BIO_2207180_risdiplam_columns =  Compounds_3x_vs_GM_DMSO[ , grep( "2207180_.*x_vs_|1949634_.*x_vs_"    , colnames(Compounds_3x_vs_GM_DMSO)  ) ] # head(BIO_2207180_risdiplam_columns)
dPSI_working = dPSI_working[ (rowSums(is.na(dPSI_working) ) != ncol(dPSI_working)), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na

head(dPSI_working, 2)
dPSI_working[1:2,]

dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_1949634_3x_vs_GM_DMSO)] <- 0
dPSI_working$GM_BIO_2207180_3x_vs_GM_DMSO[is.na(dPSI_working$GM_BIO_2207180_3x_vs_GM_DMSO)] <- 0

getwd()
fout = "Compound10_BIO_2207180_3x_vs_Risdiplam_in_DMSO.png"
png(fout, height = 1200, width = 1200)            #png(fout,height = 500, width = length(y)*50)
ggplot(dPSI_working, aes( x = GM_BIO_2207180_3x_vs_GM_DMSO, y= GM_BIO_1949634_3x_vs_GM_DMSO, color = GM_BIO_1949634_3x_vs_GM_DMSO) )  +
    geom_point(size=5) + 
    geom_text_repel(label=rownames(dPSI_working)) + 
    labs(title="DEG : internal compounds VS Risdiplam",
         y ="Risdiplam 10x dPSI vs DMSO") + 
    theme(
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
    ) + sc
dev.off()

# Ah… This is not QuickOmics, but the RANsequest analysis (all EA...) pipeline to produce the Quickomics object.
# The front end of RNAsequest is not design for clusters, since it is published (open source) and has ability to run without computational clusters. The parameter ‘parallel: slurm' in config file enable user to use computational clusters for DEG step (the most time consuming part).
# You can setup any command line tool run in the background () as:
# nohup  EArun …/config.yml &> .log &
# This way even your connection break, the program is still running. And you can check .log file for progression.
# Hope this helps. 
# O'Young

#
                 1) ==========================
PATH=/edgehpc/dept/compbio/edge_tools/RNAsequest
should be PATH=/edgehpc/dept/compbio/edge_tools/RNASequest  # upcase "S""




2) ==========================
# The file /camhpc/ngs/projects/
The file /mnt/depts/dept04/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/samplesheet.tsv 

2.1) the column "Treatment"  allow alphabetical , "_" and "." only. 
My suggetion is use  alphabetical and "_" only , such as iPSC_2060573_3uM

Hence, I replaced all "-" to "_" 

iPSC_DMSO
iPSC_2060573_3uM
iPSC_2060573_10uM
iPSC_2059811_261nM

NGN2_BIO_2006152_10x_1
NGN2_BIO_2006152_10x_2
NGN2_BIO_2006152_10x_3
regex _\d$ remove trailing for treatment


bridge_a
bridge_b
bridge_c

regex bridge_.

2.2) the same column, check name exist or not "iPSC_DMSO1" "iPSC_DMSO2" "iPSC_DMSO3"
should be iPSC_DMSO 


2.3) NO need to delete samples here, as Metadata only contains useful data


EAinit /mnt/depts/dept04/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/
EAqc /mnt/depts/dept04/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/

3) ==========================
compaireInfo.csv should uploaded to /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA20230126_0/data/  and replace the empty "compareInfo.csv"

CompareName	               Subsetting_group	Model	Covariate_levels	Group_name	Group_test	Group_ctrl	Analysis_method	Shrink_logFC	LFC_cutoff
iPSC-2060573-3uM_vs_iPSC-DMSO		        Treatment	                	Treatment	iPSC_2060573_3uM	iPSC_DMSO	DESeq2	Yes	0.2
iPSC-2060573-10uM_vs_iPSC-DMSO		      Treatment		                Treatment	iPSC_2060573_10uM	iPSC_DMSO	DESeq2	Yes	0.2


3.1) "Model"	and	"Group_name" column cannot be blank, just fill "Treatment"


3.2) "Group_test"	Group_ctrl column allow alphabetical , "_" and "." only.
My suggetion is use  alphabetical and "_" only , such as iPSC_2060573_3uM
such as iPSC_2060573_3uM and keep exact same as the column "Treatment" of "samplesheet.tsv"





https://biib.sharepoint.com/:f:/r/sites/CompBio/Shared%20Documents/TST11872%20RNA-seq%20profiling%20of%20HTT%20benchmark%20and%20Biogen%20compounds%20for%20off-target%20analysis%20OTL_91266%20-%20HTT%20Splice%20Modulator?csf=1&web=1&e=riFyqQ


https://biib.sharepoint.com/:f:/r/sites/CompBio/Shared%20Documents/TST11955%20RNA-seq%20profiling%20of%20Biogen%20HTT%20splice%20modulator%20leads%20for%20off-target%20analysis%20OTL_91266%20-%20HTT%20Splice%20Modulator?csf=1&web=1&e=hgf1He

compbio main entry
https://biib.sharepoint.com/:f:/r/sites/CompBio
https://biib.sharepoint.com/:f:/r/sites/CompBio/TST11955

TST12145
Status:
  Data Analysis - Pending
Previous Status:
  Data Generation
Archived:
  No
Requested Date:
  2023-02-23
Requested By:
  Seung Joon Lee
Approved By:
  Thomas Carlile
Last Updated By:
  Thomas Carlile
Last Updated Date:
  2023-05-22
Details
Study Title:
  Pilot study to assess splicing changes promoted by modified U1 snRNA
Protocol:
  KAPA mRNA HyperPrep
Study Participants / Point of Contact:
  Seung Joon Lee
Oracle Corporate Project Code:
  91009 Splice Modulator Platform
############################## ==================================================== >>>>>>>>>>>>>>>>    >>>>>>>>>>>>>>>>        Objectives & Background:
  Majority of the ongoing splice modulating small molecule projects aim to target 5’ splice site – which basepairing with U1 snRNA constitutes the initial step of splicing. 
  The A nucleotide bulge formed during this 5’ splice site/U1 snRNA duplex provides a druggable space to promote splicing. 
  A modified U1 snRNA has been shown to mimic splice modulator and promote cryptic splicing. 
  This will allow us to 
  
  1) discover new cryptic splice sites which can be targeted to achieve desirable outcomes (target mRNA down/up regulation), 
  2) provide a list of potential off-targets when the particular 5’ splice site is strengthened and 
  3) provide mechanistic clue of how splice modulators could influence splicing in only the subset of DSGs modulated by U1 snRNA.

############################## ==================================================== >>>>>>>>>>>>>>>>    >>>>>>>>>>>>>>>>        Experiment Design:

  293 cells will be transiently transfected with U1 snRNA plasmids. 
  48 hrs later, total RNA will be purified for bulk RNA-seq. Optionally, we could include cycloheximide, a known nonsense-mediated decay inhibitor, to protect the cryptic NMD transcripts from degradation. 
  Also, another control is to include benchmark splice modulators (BIO-2197294 and internal leads).
  
Four biological replicates will be prepared to support differential splicing analyses.

Test three different versions of U1 snRNA – 
1) Wild type, 2) A bulge repaired (target : HTT), 3) UU mismatch repair (target : MSH3).

Splice modulator treatment : need to finalize treatment duration and concentration.

############################## ==================================================== >>>>>>>>>>>>>>>>    >>>>>>>>>>>>>>>>
ZhenGao added 2023-07-16: from Joon, CHX is very toxic will kill most cell if time is long, hence treated 6 hours (maybe after 18 hour drug treamtment) , then collcect cell samples. 

############################## ==================================================== >>>>>>>>>>>>>>>>    >>>>>>>>>>>>>>>>       Data Analysis Deliverables:

  DSG and DEG analyses of U1 snRNA/splice modulator
  
Group the DSGs into different categories 

– DSGs modulated only by U1 snRNA, DSGs modulated only by splice modulators and DSGs modulated by both.

Motif enrichment analyses could be carried out to discover any feature/motif which is necessary for splice modulation by either U1 snRNA or splice modulator.

Timeline Information
Requested Completion Date:
  2023-06-30
Sample Availability Date:
  2023-03-31
Other Details
Related Existing Studies:
  TST11872, TST11955, TST12086 (HTT splice modulator RNA-seq studies)
Partial Project Support:
  -
  Requested Priority:
  2
Therapeutic Area:
  Cross-RU Exploratory Research
R&D Pipeline Stage:
  Exploratory Target
Sample Type:
  RNA
Estimated Sample Number:
  Need to be determined after scoping meeting
Summary of Findings and Follow-up:
  -

# # (Step 0) Launch FASTR splicing pipeline (Event-level analysis depends on DNAnexus splicing outputs) 

# (Step 1) Make a working dictionary (WD)
# WD=/edgehpc/dept/compbio/projects/TSTID*/splicing_events/

# (Step 2) Preparing the data:
# mkdir and cd /data folder under working dictionary (WD):

# (Step 3) download data by using a script: example please see 
# /home/ychen12/splicing/TST12086/link_to_projects/splicing_events/data/data_download.sh
# /home/zgao1
# != /edgehpc/dept/compbio/users/ychen12/splicing/TST12086/link_to_projects/splicing_events/data/data_download.sh
# != /edgehpc/dept/compbio/users/zgao1/

#
WD=/edgehpc/dept/compbio/projects/TST12188/splicing_events/
$WD/data

cd  /edgehpc/dept/compbio/projects/TST12188/splicing_events/data/
    
cp /home/ychen12/splicing/TST12086/link_to_projects/splicing_events/data/data_download.sh .    



To continue, please double check the completion of data preparation. 

(Step 3) bash main_all.sh: example please see /home/ychen12/splicing/TST12121/link_to_projects/splicing_events2/main_all.sh

cp /home/ychen12/splicing/TST12121/link_to_projects/splicing_events2/main_all.sh /edgehpc/dept/compbio/projects/TST12188/splicing_events/SH_data 
cp /home/ychen12/splicing/TST12121/link_to_projects/splicing_events2/main_all.sh /edgehpc/dept/compbio/projects/TST12188/splicing_events/GM_data 
cp /home/ychen12/splicing_test/splicing_harmonization/DIRS/create_scripts.py /edgehpc/dept/compbio/projects/TST12188/splicing_events/GM_data 
cp /home/ychen12/splicing_test/splicing_harmonization/DIRS/create_scripts.py /edgehpc/dept/compbio/projects/TST12188/splicing_events/SH_data 


python /home/ychen12/splicing_test/splicing_harmonization/DIRS/create_scripts.py -indir $WD -ref $reference -comparisonRef_name $ref_in_comparision
Traceback (most recent call last):
    File "/home/ychen12/splicing_test/splicing_harmonization/DIRS/create_scripts.py", line 15, in <module>
    for i, file in enumerate(os.listdir( main_folder + "/data/stringtie")):
    OSError: [Errno 2] No such file or directory: '/mnt/depts/dept04/compbio/projects/TST12188/splicing_events/SH_data/data/stringtie'

mv '/mnt/depts/dept04/compbio/projects/TST12188/splicing_events/SH_data/data/stringtie' to '/mnt/depts/dept04/compbio/projects/TST12188/'

mv '/mnt/depts/dept04/compbio/projects/TST12188/splicing_events/SH_data/' "/mnt/depts/dept04/compbio/projects/TST12188/SH_splicing_events/" 

# nohup wget -i DNAnexus_export_urls-20230905-123819.txt & or curl -O individual
