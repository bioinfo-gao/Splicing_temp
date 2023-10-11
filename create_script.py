#setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_downstream/")
R.Version()
#rm(list = ls()) #install.packages("ggVennDiagram") # 15min
#library("ggVennDiagram")
library(ggplot2)
require("UpSetR") # movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header = T, sep = ";")
library(dplyr)

getwd()
out_dir = "/edgehpc/dept/compbio/projects/TST12188/dnanexus/2023-09-24-4vs4-3vs3-difference/DSG_Result_difference_09-25/"
dir.create(out_dir, recursive = T)               #dir.create(din, recursive = T)
setwd(out_dir )


######
GM_4vs4_dPSI_bi =  read.csv("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/DSG_Result/splicing_offtargets_rerun_4vs4_09-25/binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", stringsAsFactors=F)   
row.names(GM_4vs4_dPSI_bi ) # row.names(dPSI_bi ) =  make.names(dPSI_bi$genesymbol, unique=TRUE) # dPSI_bi $genesymbol = NULL
row.names(GM_4vs4_dPSI_bi )  = GM_4vs4_dPSI_bi $X
GM_4vs4_dPSI_bi$X = NULL

######
SH_4vs4_dPSI_bi =  read.csv("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/DSG_Result/splicing_offtargets_rerun_4vs4_09-25/binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", stringsAsFactors=F)   
row.names(SH_4vs4_dPSI_bi  ) # row.names(SH_4vs4_dPSI_bi  ) =  make.names(SH_4vs4_dPSI_bi $genesymbol, unique=TRUE) # SH_4vs4_dPSI_bi  $genesymbol = NULL
row.names(SH_4vs4_dPSI_bi  )  = SH_4vs4_dPSI_bi $X
SH_4vs4_dPSI_bi$X = NULL
#####


######
GM_3vs3_dPSI_bi =  read.csv("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/splicing_offtargets_3vs3_09-25//binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", stringsAsFactors=F)   
row.names(GM_3vs3_dPSI_bi  ) # row.names(GM_3vs3_dPSI_bi  ) =  make.names(GM_3vs3_dPSI_bi $genesymbol, unique=TRUE) # GM_3vs3_dPSI_bi  $genesymbol = NULL
row.names(GM_3vs3_dPSI_bi  )  = GM_3vs3_dPSI_bi $X
GM_3vs3_dPSI_bi$X = NULL

######
SH_3vs3_dPSI_bi =  read.csv("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/splicing_offtargets_3vs3_09-25/binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", stringsAsFactors=F)   
row.names(SH_3vs3_dPSI_bi  ) # row.names(SH_3vs3_dPSI_bi  ) =  make.names(SH_3vs3_dPSI_bi $genesymbol, unique=TRUE) # SH_3vs3_dPSI_bi  $genesymbol = NULL
row.names(SH_3vs3_dPSI_bi  )  = SH_3vs3_dPSI_bi $X
SH_3vs3_dPSI_bi$X = NULL
#####


dim(GM_4vs4_dPSI_bi)
dim(SH_4vs4_dPSI_bi)
dim(GM_3vs3_dPSI_bi)
dim(SH_3vs3_dPSI_bi)
GM_4vs4_dPSI_bi[1:6,1:6]
SH_4vs4_dPSI_bi[1:6,1:6]
GM_3vs3_dPSI_bi[1:6,1:6]
SH_3vs3_dPSI_bi[1:6,1:6]

colnames(GM_4vs4_dPSI_bi) =  paste0( "4vs4_", colnames(GM_4vs4_dPSI_bi) )
colnames(SH_4vs4_dPSI_bi) =  paste0( "4vs4_", colnames(SH_4vs4_dPSI_bi) )
colnames(GM_3vs3_dPSI_bi) =  paste0( "3vs3_", colnames(GM_3vs3_dPSI_bi) )
colnames(SH_3vs3_dPSI_bi) =  paste0( "3vs3_", colnames(SH_3vs3_dPSI_bi) )


GM_4vs4_3vs3_dPSI_bi = merge(GM_4vs4_dPSI_bi, GM_3vs3_dPSI_bi, by = 0, all = T)
dim(GM_4vs4_3vs3_dPSI_bi)
dPSI_bi = GM_4vs4_3vs3_dPSI_bi      # 4vs4_GM_BIO_1949634_10x

SH_4vs4_3vs3_dPSI_bi = merge(SH_4vs4_dPSI_bi, SH_3vs3_dPSI_bi, by = 0, all = T)
dim(SH_4vs4_3vs3_dPSI_bi)
dPSI_bi = SH_4vs4_3vs3_dPSI_bi

head(dPSI_bi)
row.names(dPSI_bi) = dPSI_bi$Row.names
dPSI_bi$Row.names = NULL

dPSI_bi[is.na(dPSI_bi)] <- 0 

head(dPSI_bi)
dim(dPSI_bi)
colnames(dPSI_bi)


# "GM_BIO_2178782" "GM_BIO_2186827" "GM_BIO_2196772" "GM_BIO_2196895" "GM_BIO_2197306" "GM_BIO_2199562" "GM_BIO_2201042" "GM_BIO_2204984" "GM_BIO_2206678" "GM_BIO_2207180" ==== #############
dPSI_bi_working = dPSI_bi[ ,  grepl( "^...._.._..._1949634_10x"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep 4vs4_GM_BIO_1949634_10x 3vs3_GM_BIO_1949634_10x_
dPSI_bi_working = dPSI_bi[ ,  grepl( "^...._.._..._1949634_3x"     , colnames(dPSI_bi)  ) ]  # grepl NOT grep 4vs4_GM_BIO_1949634_3x 3vs3_GM_BIO_1949634_3x_


dPSI_bi_working[1:6, ] # only 2 column

# 
#  dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2178782"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2186827"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2196772"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2196895"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2197306"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
# # # 
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2199562"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2201042"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2204984"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2206678"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
#   dPSI_bi_working = dPSI_bi[ ,  grepl( "^GM_2207180"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep

###### 
dPSI_bi_working 
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na
head(dPSI_bi_working)
dPSI_bi_working[1:2,]
dim(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) # 15 types
colnames_set

getwd()
# fout = paste0("GM_5_UpSet_Ris_10x", ".png")
# fout = paste0("GM_5_UpSet_Ris_3x", ".png")

fout = paste0("SH_5_UpSet_Ris_10x", ".png")
#fout = paste0("SH_5_UpSet_Ris_3x", ".png")

# fout = paste0("5_UpSet_Ris", ".png")
 # fout = paste0("6_UpSet_BIO_2178782_vs_DMSO_and_Ris", ".png")
 # fout = paste0("7_UpSet_BIO_2186827_vs_DMSO_and_Ris", ".png")
 # fout = paste0("8_UpSet_BIO_2196772_vs_DMSO_and_Ris", ".png")
 # fout = paste0("9_UpSet_BIO_2196895_vs_DMSO_and_Ris", ".png")
 # fout = paste0("10_UpSet_BIO_2197306_vs_DMSO_and_Ris", ".png")
 # fout = paste0("11_UpSet_BIO_2199562_vs_DMSO_and_Ris", ".png")
 # fout = paste0("12_UpSet_BIO_2201042_vs_DMSO_and_Ris", ".png")
 # fout = paste0("13_UpSet_BIO_2204984_vs_DMSO_and_Ris", ".png")
 # fout = paste0("14_UpSet_BIO_2206678_vs_DMSO_and_Ris", ".png")
 # fout = paste0("15_UpSet_BIO_2207180_vs_DMSO_and_Ris", ".png")

png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50)            #png(fout,height = 500, width = length(y)*50)
upset(dPSI_bi_working , sets = colnames_set,  mainbar.y.label = "Differential Splicing Gene Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()


### "SH_1949634" "SH_2178782" "SH_2196772" "SH_2196895" "SH_2197306" "SH_2199562" "SH_2204984" "SH_2207180"
dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_1949634"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep

 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2178782"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep  # dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2186827"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2196772"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2196895"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2197306"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep

 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2199562"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep  # dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2201042"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2204984"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep  # dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2206678"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep
 dPSI_bi_working = dPSI_bi[ ,  grepl( "^SH_2207180"    , colnames(dPSI_bi)  ) ]  # grepl NOT grep

###### 
dPSI_bi_working 
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #U1_WT_EtOH_vs_DMSO_EtOH_dPSI # m[rowSums(is.na(m)) != ncol(m), ] # https://stackoverflow.com/questions/6471689/remove-rows-in-r-matrix-where-all-data-is-na
head(dPSI_bi_working)
dPSI_bi_working[1:2,]
dim(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) # 15 types
colnames_set

 fout = paste0("SH_5_UpSet_Ris", ".png")
 fout = paste0("6_UpSet_BIO_2178782_vs_DMSO_and_Ris", ".png") #  fout = paste0("7_UpSet_BIO_2186827_vs_DMSO_and_Ris", ".png")
 fout = paste0("8_UpSet_BIO_2196772_vs_DMSO_and_Ris", ".png") 
 fout = paste0("9_UpSet_BIO_2196895_vs_DMSO_and_Ris", ".png")
 fout = paste0("10_UpSet_BIO_2197306_vs_DMSO_and_Ris", ".png")
 fout = paste0("11_UpSet_BIO_2199562_vs_DMSO_and_Ris", ".png") #  fout = paste0("12_UpSet_BIO_2201042_vs_DMSO_and_Ris", ".png")
 fout = paste0("13_UpSet_BIO_2204984_vs_DMSO_and_Ris", ".png") #  fout = paste0("14_UpSet_BIO_2206678_vs_DMSO_and_Ris", ".png")
 fout = paste0("15_UpSet_BIO_207180_vs_DMSO_and_Ris", ".png")

 getwd()
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50)            #png(fout,height = 500, width = length(y)*50)
upset(dPSI_bi_working , sets = colnames_set,  mainbar.y.label = "Differential Splicing Gene Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()

