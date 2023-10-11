# 3rd_catherine4.R
#din="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG
#setwd("/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/ZG_analysis.03.splice_analysis/")
# din="/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG
#DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/ZG_analysis.03.splice_analysis/summarytable_minpadj_maxPdPSI_dPSI_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T) # no row-name

#require(gdata)
# require(xlsx)
# require(readxl)
#install.packages("openxlsx")

# xx <- pheatmap(test) # https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

library(pheatmap)
library("openxlsx")
#setwd("~/2T_Disk/Download/DATA")
getwd()



### 3rd exp
#DSG_3v3_3rd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
#dim(DSG_3v3_3rd)
DSG_3v3_3rd = read.csv("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results/analysis.03.splice_analysis_3vs3/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_3rd)
DSG_3v3_3rd[1:10, 1:10]
DSG_3v3_3rd[1:10, ]
colnames(DSG_3v3_3rd)
dim(DSG_3v3_3rd)

row.names(DSG_3v3_3rd) = DSG_3v3_3rd$X
DSG_3v3_3rd$X = NULL
DSG_3v3_3rd

getwd()
################# (1) iPSC ==========================================================================
###== (1.1) dPSI === ####
DSG_3v3_3rd_dPSI  = DSG_3v3_3rd[ , grep("^dPSI_",   colnames(DSG_3v3_3rd  ) ) ] 
DSG_3v3_3rd_dPSI[1:10,]
colnames(DSG_3v3_3rd_dPSI)

## === (1.1.1) 3x iPSC
DSG_3v3_3rd_dPSI_3x      = DSG_3v3_3rd_dPSI[   grepl("3x",   colnames(DSG_3v3_3rd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_dPSI_3x_iPSC = DSG_3v3_3rd_dPSI_3x[grepl("iPSC", colnames(DSG_3v3_3rd_dPSI_3x) ) ]
DSG_3v3_3rd_dPSI_3x_iPSC
 
DSG_3v3_3rd_dPSI_3x_iPSC = DSG_3v3_3rd_dPSI_3x_iPSC[ ,  c( "dPSI_iPSC_1755497_3x" , "dPSI_iPSC_2006152_3x" , "dPSI_iPSC_2189972_3x" , "dPSI_iPSC_2194943_3x" , "dPSI_iPSC_2195127_3x" , "dPSI_iPSC_2195327_3x" , "dPSI_iPSC_2197294_3x" )]

#dPSI_iPSC_1755497_3x_bridge
#dPSI_iPSC_2006152_3x_bridge 
# dPSI_iPSC_1755497_3x 
# dPSI_iPSC_2006152_3x 
# dPSI_iPSC_2189972_3x 
# dPSI_iPSC_2194943_3x 
# dPSI_iPSC_2195127_3x 
# dPSI_iPSC_2195327_3x 
# dPSI_iPSC_2197294_3x

DSG_3v3_3rd_dPSI_3x_iPSC$dPSI_iPSC_1755497_3x
DSG_3v3_3rd_dPSI_3x_iPSC

#DSG_3v3_3rd_dPSI_3x_iPSC_1755497_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC$dPSI_iPSC_1755497_3x )

DSG_3v3_3rd_dPSI_3x_iPSC_1755497_abs = abs( DSG_3v3_3rd_dPSI_3x_iPSC[ ,"dPSI_iPSC_1755497_3x", drop=FALSE  ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2006152_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2006152_3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2189972_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2189972_3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2194943_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2194943_3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2195127_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2195127_3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2195327_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2195327_3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs = abs(DSG_3v3_3rd_dPSI_3x_iPSC[, "dPSI_iPSC_2197294_3x", drop=FALSE ] )


DSG_3v3_3rd_dPSI_3x_iPSC_abs = cbind(DSG_3v3_3rd_dPSI_3x_iPSC_1755497_abs , DSG_3v3_3rd_dPSI_3x_iPSC_2006152_abs , DSG_3v3_3rd_dPSI_3x_iPSC_2189972_abs , DSG_3v3_3rd_dPSI_3x_iPSC_2194943_abs ,
                                DSG_3v3_3rd_dPSI_3x_iPSC_2195127_abs , DSG_3v3_3rd_dPSI_3x_iPSC_2195327_abs , DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs )

dir.create("/mnt/depts/dept04/compbio/users/zgao1/project_splicing/TST12806/All_1st_vs_2nd_vs_3rd/DSG" , recursive =T) 
setwd("/mnt/depts/dept04/compbio/users/zgao1/project_splicing/TST12806/All_1st_vs_2nd_vs_3rd/DSG")
getwd()

fout = "DSG_3v3_3rd_dPSI_3x_iPSC_abs_with_name.csv"
write.csv(DSG_3v3_3rd_dPSI_3x_iPSC_abs,file = fout, quote = F, row.names = T)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pheatmap(DSG_3v3_3rd_dPSI_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_3x_iPSC_abs.pdf")

#### (1.1.2) 10x iPSC
DSG_3v3_3rd_dPSI_10x      = DSG_3v3_3rd_dPSI[   grepl("10x",   colnames(DSG_3v3_3rd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_dPSI_10x_iPSC = DSG_3v3_3rd_dPSI_10x[grepl("iPSC", colnames(DSG_3v3_3rd_dPSI_10x) ) ]
DSG_3v3_3rd_dPSI_10x_iPSC

DSG_3v3_3rd_dPSI_10x_iPSC_1755497_abs = abs( DSG_3v3_3rd_dPSI_10x_iPSC[ ,"dPSI_iPSC_1755497_10x", drop=FALSE  ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2006152_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2006152_10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2189972_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2189972_10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2194943_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2194943_10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2195127_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2195127_10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2195327_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2195327_10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[, "dPSI_iPSC_2197294_10x", drop=FALSE ] )

DSG_3v3_3rd_dPSI_10x_iPSC_abs = cbind(DSG_3v3_3rd_dPSI_10x_iPSC_1755497_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2006152_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2189972_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2194943_abs ,
                                DSG_3v3_3rd_dPSI_10x_iPSC_2195127_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2195327_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs )

fout = "DSG_3v3_3rd_dPSI_10x_iPSC_abs.csv"
write.csv(DSG_3v3_3rd_dPSI_10x_iPSC_abs,file = fout, quote = F, row.names = T)

pheatmap(DSG_3v3_3rd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_10x_iPSC_abs.pdf")

###== (1.2) significance === ####
#DSG_3v3_3rd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_3rd)
DSG_3v3_3rd[1:10, 1:10]
DSG_3v3_3rd[1:10, ]
colnames(DSG_3v3_3rd)
dim(DSG_3v3_3rd)
row.names(DSG_3v3_3rd) = DSG_3v3_3rd$X
DSG_3v3_3rd$X = NULL

DSG_3v3_3rd
colnames(DSG_3v3_3rd)
head(DSG_3v3_3rd)

### 1.2 Q and Prob
DSG_3v3_3rd_Q_P  = DSG_3v3_3rd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_3rd  ) ) ] 
DSG_3v3_3rd_Q_P[1:10,]
colnames(DSG_3v3_3rd_Q_P)
head(DSG_3v3_3rd_Q_P)

## === (1.2.1) 3x iPSC
DSG_3v3_3rd_Q_P_3x      = DSG_3v3_3rd_Q_P[   grepl("3x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_3x_iPSC = DSG_3v3_3rd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_3x) ) ]
DSG_3v3_3rd_Q_P_3x_iPSC
colnames(DSG_3v3_3rd_Q_P_3x_iPSC)
head(DSG_3v3_3rd_Q_P_3x_iPSC)


DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs ,
                                DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs , 
                                DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_3rd_Q_P_3x_iPSC)
fout = "DSG_3v3_3rd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_3x_iPSC.pdf")

#### (1.2.2) 10x iPSC
DSG_3v3_3rd_Q_P_10x      = DSG_3v3_3rd_Q_P[   grepl("10x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_10x_iPSC = DSG_3v3_3rd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_10x) ) ]
DSG_3v3_3rd_Q_P_10x_iPSC

DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs ,
                                 DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs , 
                                 DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_3rd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_10x_iPSC.pdf")

################# (2) NGN2 ==========================================================================
###== (2.1) dPSI === ####
# colnames(DSG_3v3_3rd  )
# head(DSG_3v3_3rd, 2  )
# DSG_3v3_3rd_dPSI  = DSG_3v3_3rd[ , grep("^dPSI_",   colnames(DSG_3v3_3rd  ) ) ] 
# DSG_3v3_3rd_dPSI[1:10,]
# colnames(DSG_3v3_3rd_dPSI)

## === (2.1.1) 3x NGN2
DSG_3v3_3rd_dPSI_3x      = DSG_3v3_3rd_dPSI[   grepl("3x",   colnames(DSG_3v3_3rd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_dPSI_3x_NGN2 = DSG_3v3_3rd_dPSI_3x[grepl("NGN2", colnames(DSG_3v3_3rd_dPSI_3x) ) ]
DSG_3v3_3rd_dPSI_3x_NGN2

DSG_3v3_3rd_dPSI_3x_NGN2_1755497_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`1755497`, "dPSI_NGN2.1755497.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2006152_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2006152`, "dPSI_NGN2.2006152.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2189972_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2189972`, "dPSI_NGN2.2189972.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2194943_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2194943`, "dPSI_NGN2.2194943.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2195127_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2195127`, "dPSI_NGN2.2195127.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2195327_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2195327`, "dPSI_NGN2.2195327.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_2197294_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`2197294`, "dPSI_NGN2.2197294.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_6527_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`6527`, "dPSI_NGN2.6527.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_6866_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`6866`, "dPSI_NGN2.6866.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_6960_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`6960`, "dPSI_NGN2.6960.3x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_NGN2_9701_abs = abs(DSG_3v3_3rd_dPSI_3x_NGN2[NGN2_data$`9701`, "dPSI_NGN2.9701.3x", drop=FALSE ] )

DSG_3v3_3rd_dPSI_3x_NGN2_abs = cbind(DSG_3v3_3rd_dPSI_3x_NGN2_1755497_abs , DSG_3v3_3rd_dPSI_3x_NGN2_2006152_abs , DSG_3v3_3rd_dPSI_3x_NGN2_2189972_abs , DSG_3v3_3rd_dPSI_3x_NGN2_2194943_abs ,
                                     DSG_3v3_3rd_dPSI_3x_NGN2_2195127_abs , DSG_3v3_3rd_dPSI_3x_NGN2_2195327_abs , DSG_3v3_3rd_dPSI_3x_NGN2_2197294_abs , DSG_3v3_3rd_dPSI_3x_NGN2_6527_abs , 
                                     DSG_3v3_3rd_dPSI_3x_NGN2_6866_abs , DSG_3v3_3rd_dPSI_3x_NGN2_6960_abs , DSG_3v3_3rd_dPSI_3x_NGN2_9701_abs )
fout = "DSG_3v3_3rd_dPSI_3x_NGN2_abs.csv"
write.csv(DSG_3v3_3rd_dPSI_3x_NGN2_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_dPSI_3x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_3x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_3x_NGN2_abs.pdf")

#### (2.1.2) 10x NGN2
DSG_3v3_3rd_dPSI_10x      = DSG_3v3_3rd_dPSI[   grepl("10x",   colnames(DSG_3v3_3rd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_dPSI_10x_NGN2 = DSG_3v3_3rd_dPSI_10x[grepl("NGN2", colnames(DSG_3v3_3rd_dPSI_10x) ) ]
DSG_3v3_3rd_dPSI_10x_NGN2

DSG_3v3_3rd_dPSI_10x_NGN2_1755497_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`1755497`, "dPSI_NGN2.1755497.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2006152_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2006152`, "dPSI_NGN2.2006152.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2189972_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2189972`, "dPSI_NGN2.2189972.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2194943_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2194943`, "dPSI_NGN2.2194943.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2195127_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2195127`, "dPSI_NGN2.2195127.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2195327_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2195327`, "dPSI_NGN2.2195327.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_2197294_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`2197294`, "dPSI_NGN2.2197294.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_6527_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`6527`, "dPSI_NGN2.6527.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_6866_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`6866`, "dPSI_NGN2.6866.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_6960_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`6960`, "dPSI_NGN2.6960.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_NGN2_9701_abs = abs(DSG_3v3_3rd_dPSI_10x_NGN2[NGN2_data$`9701`, "dPSI_NGN2.9701.10x", drop=FALSE ] )

DSG_3v3_3rd_dPSI_10x_NGN2_abs = cbind(DSG_3v3_3rd_dPSI_10x_NGN2_1755497_abs , DSG_3v3_3rd_dPSI_10x_NGN2_2006152_abs , DSG_3v3_3rd_dPSI_10x_NGN2_2189972_abs , DSG_3v3_3rd_dPSI_10x_NGN2_2194943_abs ,
                                      DSG_3v3_3rd_dPSI_10x_NGN2_2195127_abs , DSG_3v3_3rd_dPSI_10x_NGN2_2195327_abs , DSG_3v3_3rd_dPSI_10x_NGN2_2197294_abs , DSG_3v3_3rd_dPSI_10x_NGN2_6527_abs , 
                                      DSG_3v3_3rd_dPSI_10x_NGN2_6866_abs , DSG_3v3_3rd_dPSI_10x_NGN2_6960_abs , DSG_3v3_3rd_dPSI_10x_NGN2_9701_abs )
fout = "DSG_3v3_3rd_dPSI_10x_NGN2_abs.csv"
write.csv(DSG_3v3_3rd_dPSI_10x_NGN2_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_dPSI_10x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_10x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_10x_NGN2_abs.pdf")

###== (2.2) significance === ####
### 2.2 Q and Prob
DSG_3v3_3rd_Q_P  = DSG_3v3_3rd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_3rd  ) ) ] 
DSG_3v3_3rd_Q_P[1:10,]
colnames(DSG_3v3_3rd_Q_P)
head(DSG_3v3_3rd_Q_P)

## === (2.2.1) 3x iPSC
DSG_3v3_3rd_Q_P_3x      = DSG_3v3_3rd_Q_P[   grepl("3x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_3x_iPSC = DSG_3v3_3rd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_3x) ) ]
DSG_3v3_3rd_Q_P_3x_iPSC
colnames(DSG_3v3_3rd_Q_P_3x_iPSC)
head(DSG_3v3_3rd_Q_P_3x_iPSC)
head(NGN2_data)

DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[NGN2_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs ,
                                    DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs , 
                                    DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_3rd_Q_P_3x_iPSC)
fout = "DSG_3v3_3rd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_3x_iPSC.pdf")

#### (2.2.2) 10x iPSC
DSG_3v3_3rd_Q_P_10x      = DSG_3v3_3rd_Q_P[   grepl("10x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_10x_iPSC = DSG_3v3_3rd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_10x) ) ]
DSG_3v3_3rd_Q_P_10x_iPSC

DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[NGN2_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs ,
                                     DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs , 
                                     DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_3rd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_10x_iPSC.pdf")


################# (3) 2197294 ==========================================================================
dim(Data_2197294)
head(Data_2197294)
tail(Data_2197294)
colnames(Data_2197294) # [1] "iPSC_3x_1st"                  "iPSC_3x_3rd"                  "iPSC_3x_supported_by_10x_1st" "iPSC_3x_supported_by_10x_3rd" "NGN2_3x_1st"                 
#[6] "NGN2_3x_3rd"                  "NGN2_3x_supported_by_10x_1st" "NGN2_3x_supported_by_10x_3rd"
Data_2197294[,1] # 41] " "          NA     

###== (3.1) dPSI === ####
DSG_3v3_1st_2197294  = DSG_3v3_1st[ , grep("2197294",   colnames(DSG_3v3_1st  ) ) ] 
DSG_3v3_3rd_2197294  = DSG_3v3_3rd[ , grep("2197294",   colnames(DSG_3v3_3rd  ) ) ] 
DSG_3v3_1st_2197294
colnames(DSG_3v3_1st_2197294)
DSG_3v3_3rd_2197294
colnames(DSG_3v3_3rd_2197294)

## === (3.1.1) 3x 2197294


DSG_3v3_1st_dPSI_3x_iPSC_2197294_abs = abs(DSG_3v3_1st_2197294[Data_2197294$iPSC_3x_1st, "dPSI_TST11872_iPSC.BIO.2006152.45nM", drop=FALSE ] )
DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs = abs(DSG_3v3_3rd_2197294[Data_2197294$iPSC_3x_3rd, "dPSI_iPSC.2197294.3x", drop=FALSE ] )
DSG_3v3_1st_dPSI_3x_iPSC_2197294_abs
DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs

DSG_3v3_1st_dPSI_10x_iPSC_2197294_abs = abs(DSG_3v3_1st_2197294[Data_2197294$iPSC_10x_1st, "dPSI_TST11872_iPSC.BIO.2006152.150nM", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_2197294[Data_2197294$iPSC_10x_3rd, "dPSI_iPSC.2197294.10x", drop=FALSE ] )

DSG_3v3_1st_dPSI_10x_iPSC_2197294_abs
DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs


DSG_3v3_dPSI_iPSC_2197294_abs = cbind( )

fout = "DSG_3v3_3rd_dPSI_2197294_abs.csv"
write.csv(DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_3x_iPSC_2197294_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_2197294_abs.pdf")

#### (3.1.2) 10x iPSC
DSG_3v3_3rd_dPSI_10x      = DSG_3v3_3rd_dPSI[   grepl("10x",   colnames(DSG_3v3_3rd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_dPSI_10x_iPSC = DSG_3v3_3rd_dPSI_10x[grepl("iPSC", colnames(DSG_3v3_3rd_dPSI_10x) ) ]
DSG_3v3_3rd_dPSI_10x_iPSC

DSG_3v3_3rd_dPSI_10x_iPSC_1755497_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`1755497`, "dPSI_iPSC.1755497.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2006152_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2006152`, "dPSI_iPSC.2006152.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2189972_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2189972`, "dPSI_iPSC.2189972.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2194943_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2194943`, "dPSI_iPSC.2194943.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2195127_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2195127`, "dPSI_iPSC.2195127.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2195327_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2195327`, "dPSI_iPSC.2195327.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`2197294`, "dPSI_iPSC.2197294.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_6527_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`6527`, "dPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_6866_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`6866`, "dPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_6960_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`6960`, "dPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_3rd_dPSI_10x_iPSC_9701_abs = abs(DSG_3v3_3rd_dPSI_10x_iPSC[iPSC_data$`9701`, "dPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_3rd_dPSI_10x_iPSC_abs = cbind(DSG_3v3_3rd_dPSI_10x_iPSC_1755497_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2006152_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2189972_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2194943_abs ,
                                      DSG_3v3_3rd_dPSI_10x_iPSC_2195127_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2195327_abs , DSG_3v3_3rd_dPSI_10x_iPSC_2197294_abs , DSG_3v3_3rd_dPSI_10x_iPSC_6527_abs , 
                                      DSG_3v3_3rd_dPSI_10x_iPSC_6866_abs , DSG_3v3_3rd_dPSI_10x_iPSC_6960_abs , DSG_3v3_3rd_dPSI_10x_iPSC_9701_abs )

fout = "DSG_3v3_3rd_dPSI_10x_iPSC_abs.csv"
write.csv(DSG_3v3_3rd_dPSI_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_dPSI_10x_iPSC_abs.pdf")

###== (3.2) significance === ####
DSG_3v3_3rd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_3rd)
DSG_3v3_3rd[1:10, 1:10]
DSG_3v3_3rd[1:10, ]
colnames(DSG_3v3_3rd)
dim(DSG_3v3_3rd)
row.names(DSG_3v3_3rd) = DSG_3v3_3rd$X
DSG_3v3_3rd$X = NULL

DSG_3v3_3rd
colnames(DSG_3v3_3rd)
head(DSG_3v3_3rd)

### 1.2 Q and Prob
DSG_3v3_3rd_Q_P  = DSG_3v3_3rd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_3rd  ) ) ] 
DSG_3v3_3rd_Q_P[1:10,]
colnames(DSG_3v3_3rd_Q_P)
head(DSG_3v3_3rd_Q_P)

## === (3.2.1) 3x iPSC
DSG_3v3_3rd_Q_P_3x      = DSG_3v3_3rd_Q_P[   grepl("3x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_3x_iPSC = DSG_3v3_3rd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_3x) ) ]
DSG_3v3_3rd_Q_P_3x_iPSC
colnames(DSG_3v3_3rd_Q_P_3x_iPSC)
head(DSG_3v3_3rd_Q_P_3x_iPSC)
head(iPSC_data)

DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_3x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_3x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2194943_abs ,
                                    DSG_3v3_3rd_Q_P_3x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_3x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6527_abs , 
                                    DSG_3v3_3rd_Q_P_3x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_3x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_3rd_Q_P_3x_iPSC)
fout = "DSG_3v3_3rd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_3x_iPSC.pdf")

#### (3.2.2) 10x iPSC
DSG_3v3_3rd_Q_P_10x      = DSG_3v3_3rd_Q_P[   grepl("10x",   colnames(DSG_3v3_3rd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_3rd_Q_P_10x_iPSC = DSG_3v3_3rd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_3rd_Q_P_10x) ) ]
DSG_3v3_3rd_Q_P_10x_iPSC

DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`1755497`, "minpadj_maxPdPSI_iPSC.1755497.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2006152`, "minpadj_maxPdPSI_iPSC.2006152.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2189972`, "minpadj_maxPdPSI_iPSC.2189972.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2194943`, "minpadj_maxPdPSI_iPSC.2194943.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2195127`, "minpadj_maxPdPSI_iPSC.2195127.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2195327`, "minpadj_maxPdPSI_iPSC.2195327.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`2197294`, "minpadj_maxPdPSI_iPSC.2197294.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_3rd_Q_P_10x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_3rd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_3rd_Q_P_10x_iPSC_1755497_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2006152_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2189972_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2194943_abs ,
                                     DSG_3v3_3rd_Q_P_10x_iPSC_2195127_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2195327_abs , DSG_3v3_3rd_Q_P_10x_iPSC_2197294_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6527_abs , 
                                     DSG_3v3_3rd_Q_P_10x_iPSC_6866_abs , DSG_3v3_3rd_Q_P_10x_iPSC_6960_abs , DSG_3v3_3rd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_3rd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_3rd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_3rd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_3rd_Q_P_10x_iPSC.pdf")



#final[rowSums(is.na(final[ , 5:6])) == 0, ]
# (base) /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code zgao1@camhpcps01 $ 
#   $  tree | grep heatmap
# │   ├── run.02.DSG_counts_minmaxtable_07-26-merged-df-for-heatmap.R
# │   ├── run.02.DSG_counts_minmaxtable_plot_2022-07-25-for-correlation-heatmap.R
# │   ├── pheatmap_comparions_1.pdf
# ├── Miyoung_Shin_heatmap.R
# ├── run.02-1.DSG_counts_plot_ZG_2022-07-26-similarity-heatmap.R
# ├── run.02-1.DSG_counts_plot_ZG_2022-07-27-similarity-heatmap-CAM.R

### 2197294

DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_3rd/summarytable_with_gene-1st_experiment.csv", header = T)
dim(DSG_3v3_1st)
DSG_3v3_1st[1:10, 1:10]
DSG_3v3_1st[1:10, ]
# rownames(minmaxall) = minmaxall$genesymbol;   minmaxall$genesymbol = NULL 
# df <- df %>% dplyr:: select(starts_with("ABC"))
# df[, grep("ABC", names(df)), with = FALSE]
row.names(DSG_3v3_1st) = DSG_3v3_1st$X
DSG_3v3_1st$X = NULL

DSG_3v3_1st  = DSG_3v3_1st[ , grep("^dPSI_",   colnames(DSG_3v3_1st  ) ) ] 
colnames(DSG_3v3_1st) = gsub("dPSI_TST11872","1st_exp", colnames(DSG_3v3_1st)) 
colnames(DSG_3v3_1st) = gsub(".BIO",         "", colnames(DSG_3v3_1st)) 
DSG_3v3_1st[1:10, 1:10]

DSG_3v3_1st_2197294           = DSG_3v3_1st[      , grep("2006152", colnames(DSG_3v3_1st      ) ) ] 
colnames(DSG_3v3_1st_2197294) = gsub("150nM",   "10x" , colnames(DSG_3v3_1st_2197294)) 
colnames(DSG_3v3_1st_2197294) = gsub("45nM",    "3x"  , colnames(DSG_3v3_1st_2197294)) 
colnames(DSG_3v3_1st_2197294) = gsub("200",     ""    , colnames(DSG_3v3_1st_2197294)) 
dim(DSG_3v3_1st_2197294)
head(DSG_3v3_1st_2197294)

### 3rd exp
DSG_3v3_3rd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_3rd)
DSG_3v3_3rd[1:10, 1:10]
DSG_3v3_3rd[1:10, ]

row.names(DSG_3v3_3rd) = DSG_3v3_3rd$X
DSG_3v3_3rd$X = NULL

DSG_3v3_3rd  = DSG_3v3_3rd[ , grep("^dPSI_",   colnames(DSG_3v3_3rd  ) ) ] 
colnames(DSG_3v3_3rd) = gsub("dPSI","3rd_exp", colnames(DSG_3v3_3rd)) #colnames(DSG_3v3_3rd) = gsub(".BIO",         "", colnames(DSG_3v3_3rd)) 
DSG_3v3_3rd[1:10, 1:10]

DSG_3v3_3rd_2197294           = DSG_3v3_3rd[      , grep("2197294", colnames(DSG_3v3_3rd      ) ) ] 
dim(DSG_3v3_3rd_2197294)
head(DSG_3v3_3rd_2197294)

DSG_3v3_2197294 = merge(DSG_3v3_1st_2197294, DSG_3v3_3rd_2197294, by = 0, all = T)
head(DSG_3v3_2197294)
tail(DSG_3v3_2197294)
DSG_3v3_2197294[ 500:550, ]
dim(DSG_3v3_2197294) # 3927
class(DSG_3v3_2197294)

row.names(DSG_3v3_2197294) =  DSG_3v3_2197294$Row.names
DSG_3v3_2197294$Row.names  = NULL

dim(DSG_3v3_2197294)

DSG_3v3_2197294[ is.na(DSG_3v3_2197294) ]  = 0 # df[rowSums(is.na(df)) == 0, ]
DSG_3v3_2197294        =  DSG_3v3_2197294[ , c("1st_exp_iPSC.2197294.10x", "3rd_exp_iPSC.2197294.10x", "1st_exp_iPSC.2197294.3x", "3rd_exp_iPSC.2197294.3x", "1st_exp_NGN2.2197294.10x", "3rd_exp_NGN2.2197294.10x", "1st_exp_NGN2.2197294.3x", "3rd_exp_NGN2.2197294.3x")]
DSG_3v3_2197294.abs.mx =  abs(as.matrix(DSG_3v3_2197294))
DSG_3v3_2197294        =  as.data.frame(DSG_3v3_2197294.abs.mx)
DSG_3v3_2197294        = DSG_3v3_2197294[order(DSG_3v3_2197294$`1st_exp_iPSC.2197294.10x`,decreasing = T) ,]
DSG_3v3_2197294 = DSG_3v3_2197294[ rowSums(DSG_3v3_2197294[]) > 0.1,]

head(DSG_3v3_2197294)
tail(DSG_3v3_2197294)
dim(DSG_3v3_2197294)

getwd()
fout = "compound_2197294_induced_DSG_in_1st_and_3rd_experiments1.csv"
write.csv(DSG_3v3_2197294,file = fout,quote = F, row.names = T)

library(pheatmap)
#pheatmap(DSG_3v3_2197294 ,fontsize=9, fontsize_row=6, cluster_rows = T,         cluster_cols = F) 
pheatmap(DSG_3v3_2197294 ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
#pheatmap(DSG_3v3_2197294 ,fontsize=9, fontsize_row=6) 


# overlap_1st_v_2nd_2_for_catherine3.R
#din="/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG
#setwd("/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/code/ZG_analysis.03.splice_analysis/")
# din="/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG
#DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/ZG_analysis.03.splice_analysis/summarytable_minpadj_maxPdPSI_dPSI_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T) # no row-name

#require(gdata)
# require(xlsx)
# require(readxl)
#install.packages("openxlsx")

# xx <- pheatmap(test) # https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }

library(pheatmap)
library("openxlsx")
#setwd("~/2T_Disk/Download/DATA")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/")
#csv_files <- list.files(pattern = "\\.csv$")
#iPSC_data <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_ comparison_short.xlsx", sheet = 1, colNames = TRUE) 
iPSC_data <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_comparison_short1.xlsx", sheet = 1, colNames = TRUE) 
dim(iPSC_data)
head(iPSC_data)
tail(iPSC_data)
iPSC_data[,1]


#NGN2_data <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_ comparison_short.xlsx", sheet = 2, colNames = TRUE) 
NGN2_data <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_comparison_short1.xlsx", sheet = 2, colNames = TRUE) 
dim(NGN2_data)
head(NGN2_data)
tail(NGN2_data)
NGN2_data[,1]


#Data_6152 <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_ comparison_short.xlsx", sheet = 3, colNames = TRUE) 
Data_6152 <- openxlsx::read.xlsx("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/RNAseq2_3v3 DSG compiled lis_comparison_short1.xlsx", sheet = 3, colNames = TRUE)[, 1:8] 
dim(Data_6152)
head(Data_6152)
tail(Data_6152)
Data_6152[,1] # 41] " "          NA     

### 1st exp
DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/summarytable_with_gene-1st_experiment.csv", header = T)
dim(DSG_3v3_1st)
DSG_3v3_1st[1:10, 1:10]
DSG_3v3_1st[1:10, ]

row.names(DSG_3v3_1st) = DSG_3v3_1st$X
DSG_3v3_1st$X = NULL
# 
#          DSG_3v3_1st  = DSG_3v3_1st[ , grep("^dPSI_",   colnames(DSG_3v3_1st  ) ) ] 
# colnames(DSG_3v3_1st) = gsub("dPSI_TST11872","1st_exp", colnames(DSG_3v3_1st)) 
# colnames(DSG_3v3_1st) = gsub(".BIO",         "", colnames(DSG_3v3_1st)) 
# colnames(DSG_3v3_1st)
# 
# # "1st_exp_iPSC.1755497.12nM"   "1st_exp_iPSC.1755497.40nM"   "1st_exp_iPSC.1949634.10uM"   "1st_exp_iPSC.1949634.3uM"    "1st_exp_iPSC.2006152.150nM"  "1st_exp_iPSC.2006152.45nM"   "1st_exp_iPSC.2059811.261nM"  "1st_exp_iPSC.2059811.870nM" 
# # "1st_exp_iPSC.2060573.10uM"   "1st_exp_iPSC.2060573.3uM"    "1st_exp_iPSC.2060884.1000nM" "1st_exp_iPSC.2060884.3700nM" "1st_exp_iPSC.2070692.1680nM" "1st_exp_iPSC.2070692.5600nM" "1st_exp_iPSC.2135644.1uM"    "1st_exp_iPSC.2135644.300nM" 
# # "1st_exp_iPSC.2136770.1000nM" "1st_exp_iPSC.2136770.300nM"  "1st_exp_NGN2.1755497.12nM"   "1st_exp_NGN2.1755497.40nM"   "1st_exp_NGN2.1949634.10uM"   "1st_exp_NGN2.1949634.3uM"    "1st_exp_NGN2.2006152.150nM"  "1st_exp_NGN2.2006152.45nM"  
# # "1st_exp_NGN2.2059811.261nM"  "1st_exp_NGN2.2059811.870nM"  "1st_exp_NGN2.2060573.10uM"   "1st_exp_NGN2.2060573.3uM"    "1st_exp_NGN2.2060884.1000nM" "1st_exp_NGN2.2060884.3700nM" "1st_exp_NGN2.2070692.1680nM" "1st_exp_NGN2.2070692.5600nM"
# # "1st_exp_NGN2.2135644.1uM"    "1st_exp_NGN2.2135644.300nM"  "1st_exp_NGN2.2136770.1000nM" "1st_exp_NGN2.2136770.300nM" 
# 
# colnames(DSG_3v3_1st) = c( "1st_exp_iPSC.5497.3x"   , "1st_exp_iPSC.5497.10x" , "1st_exp_iPSC.9634.10x"  , "1st_exp_iPSC.9634.3x"  , "1st_exp_iPSC.6152.10x"  , "1st_exp_iPSC.6152.3x"  , "1st_exp_iPSC.9811.3x"   , "1st_exp_iPSC.9811.10x" , 
#                            "1st_exp_iPSC.0573.10x"  , "1st_exp_iPSC.0573.3x"  , "1st_exp_iPSC.0884.3x"   , "1st_exp_iPSC.0884.10x" , "1st_exp_iPSC.0692.3x"   , "1st_exp_iPSC.0692.10x" , "1st_exp_iPSC.5644.10x"  , "1st_exp_iPSC.5644.3x"  , 
#                            "1st_exp_iPSC.6770.10x"  , "1st_exp_iPSC.6770.3x"  , "1st_exp_NGN2.5497.3x"   , "1st_exp_NGN2.5497.10x" , "1st_exp_NGN2.9634.10x"  , "1st_exp_NGN2.9634.3x"  , "1st_exp_NGN2.6152.10x"  , "1st_exp_NGN2.6152.3x"  , 
#                            "1st_exp_NGN2.9811.3x"   , "1st_exp_NGN2.9811.10x" , "1st_exp_NGN2.0573.10x"  , "1st_exp_NGN2.0573.3x"  , "1st_exp_NGN2.0884.3x"   , "1st_exp_NGN2.0884.10x" , "1st_exp_NGN2.0692.3x"   , "1st_exp_NGN2.0692.10x" ,
#                            "1st_exp_NGN2.5644.10x"  , "1st_exp_NGN2.5644.3x"  , "1st_exp_NGN2.6770.10x"  , "1st_exp_NGN2.6770.3x") 
# 
# head(DSG_3v3_1st)
# 
# DSG_3v3_1st = DSG_3v3_1st[ ,  c("1st_exp_iPSC.0573.10x"  , "1st_exp_iPSC.0573.3x"  , "1st_exp_iPSC.0692.10x"  , "1st_exp_iPSC.0692.3x"  , 
#                                 "1st_exp_iPSC.0884.10x"  , "1st_exp_iPSC.0884.3x"  , "1st_exp_iPSC.5497.10x"  , "1st_exp_iPSC.5497.3x"  , 
#                                 "1st_exp_iPSC.5644.10x"  , "1st_exp_iPSC.5644.3x"  , "1st_exp_iPSC.6152.10x"  , "1st_exp_iPSC.6152.3x"  ,
#                                 "1st_exp_iPSC.6770.10x"  , "1st_exp_iPSC.6770.3x"  , "1st_exp_iPSC.9634.10x"  , "1st_exp_iPSC.9634.3x"  ,
#                                 "1st_exp_iPSC.9811.10x"  , "1st_exp_iPSC.9811.3x"  , "1st_exp_NGN2.0573.10x"  , "1st_exp_NGN2.0573.3x"  ,
#                                 "1st_exp_NGN2.0692.10x"  , "1st_exp_NGN2.0692.3x"  , "1st_exp_NGN2.0884.10x"  , "1st_exp_NGN2.0884.3x"  , 
#                                 "1st_exp_NGN2.5497.10x"  , "1st_exp_NGN2.5497.3x"  , "1st_exp_NGN2.5644.10x"  , "1st_exp_NGN2.5644.3x"  , 
#                                 "1st_exp_NGN2.6152.10x"  , "1st_exp_NGN2.6152.3x"  , "1st_exp_NGN2.6770.10x"  , "1st_exp_NGN2.6770.3x"  , 
#                                 "1st_exp_NGN2.9634.10x"  , "1st_exp_NGN2.9634.3x"  , "1st_exp_NGN2.9811.10x"  , "1st_exp_NGN2.9811.3x"  )    ]
# 
# head(DSG_3v3_1st)
# dim(DSG_3v3_1st)


### 2nd exp
DSG_3v3_2nd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_2nd)
DSG_3v3_2nd[1:10, 1:10]
DSG_3v3_2nd[1:10, ]
colnames(DSG_3v3_2nd)
dim(DSG_3v3_2nd)

row.names(DSG_3v3_2nd) = DSG_3v3_2nd$X
DSG_3v3_2nd$X = NULL

getwd()
################# (1) iPSC ==========================================================================
###== (1.1) dPSI === ####
DSG_3v3_2nd_dPSI  = DSG_3v3_2nd[ , grep("^dPSI_",   colnames(DSG_3v3_2nd  ) ) ] 
DSG_3v3_2nd_dPSI[1:10,]
colnames(DSG_3v3_2nd_dPSI)

## === (1.1.1) 3x iPSC
DSG_3v3_2nd_dPSI_3x      = DSG_3v3_2nd_dPSI[   grepl("3x",   colnames(DSG_3v3_2nd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_dPSI_3x_iPSC = DSG_3v3_2nd_dPSI_3x[grepl("iPSC", colnames(DSG_3v3_2nd_dPSI_3x) ) ]
DSG_3v3_2nd_dPSI_3x_iPSC

DSG_3v3_2nd_dPSI_3x_iPSC_4088_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`4088`, "dPSI_iPSC.4088.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_4090_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`4090`, "dPSI_iPSC.4090.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_4714_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`4714`, "dPSI_iPSC.4714.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_4741_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`4741`, "dPSI_iPSC.4741.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_4748_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`4748`, "dPSI_iPSC.4748.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_5420_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`5420`, "dPSI_iPSC.5420.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`6152`, "dPSI_iPSC.6152.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_6527_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`6527`, "dPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_6866_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`6866`, "dPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_6960_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`6960`, "dPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_9701_abs = abs(DSG_3v3_2nd_dPSI_3x_iPSC[iPSC_data$`9701`, "dPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_2nd_dPSI_3x_iPSC_abs = cbind(DSG_3v3_2nd_dPSI_3x_iPSC_4088_abs , DSG_3v3_2nd_dPSI_3x_iPSC_4090_abs , DSG_3v3_2nd_dPSI_3x_iPSC_4714_abs , DSG_3v3_2nd_dPSI_3x_iPSC_4741_abs ,
                                DSG_3v3_2nd_dPSI_3x_iPSC_4748_abs , DSG_3v3_2nd_dPSI_3x_iPSC_5420_abs , DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs , DSG_3v3_2nd_dPSI_3x_iPSC_6527_abs , 
                                DSG_3v3_2nd_dPSI_3x_iPSC_6866_abs , DSG_3v3_2nd_dPSI_3x_iPSC_6960_abs , DSG_3v3_2nd_dPSI_3x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_dPSI_3x_iPSC_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_3x_iPSC_abs.pdf")

#### (1.1.2) 10x iPSC
DSG_3v3_2nd_dPSI_10x      = DSG_3v3_2nd_dPSI[   grepl("10x",   colnames(DSG_3v3_2nd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_dPSI_10x_iPSC = DSG_3v3_2nd_dPSI_10x[grepl("iPSC", colnames(DSG_3v3_2nd_dPSI_10x) ) ]
DSG_3v3_2nd_dPSI_10x_iPSC

DSG_3v3_2nd_dPSI_10x_iPSC_4088_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4088`, "dPSI_iPSC.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4090_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4090`, "dPSI_iPSC.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4714_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4714`, "dPSI_iPSC.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4741_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4741`, "dPSI_iPSC.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4748_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4748`, "dPSI_iPSC.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_5420_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`5420`, "dPSI_iPSC.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6152`, "dPSI_iPSC.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6527_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6527`, "dPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6866_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6866`, "dPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6960_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6960`, "dPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_9701_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`9701`, "dPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_dPSI_10x_iPSC_abs = cbind(DSG_3v3_2nd_dPSI_10x_iPSC_4088_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4090_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4714_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4741_abs ,
                                DSG_3v3_2nd_dPSI_10x_iPSC_4748_abs , DSG_3v3_2nd_dPSI_10x_iPSC_5420_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6527_abs , 
                                DSG_3v3_2nd_dPSI_10x_iPSC_6866_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6960_abs , DSG_3v3_2nd_dPSI_10x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_dPSI_10x_iPSC_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_10x_iPSC_abs.pdf")

###== (1.2) significance === ####
DSG_3v3_2nd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_2nd)
DSG_3v3_2nd[1:10, 1:10]
DSG_3v3_2nd[1:10, ]
colnames(DSG_3v3_2nd)
dim(DSG_3v3_2nd)
row.names(DSG_3v3_2nd) = DSG_3v3_2nd$X
DSG_3v3_2nd$X = NULL

DSG_3v3_2nd
colnames(DSG_3v3_2nd)
head(DSG_3v3_2nd)

### 1.2 Q and Prob
DSG_3v3_2nd_Q_P  = DSG_3v3_2nd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_2nd  ) ) ] 
DSG_3v3_2nd_Q_P[1:10,]
colnames(DSG_3v3_2nd_Q_P)
head(DSG_3v3_2nd_Q_P)

## === (1.2.1) 3x iPSC
DSG_3v3_2nd_Q_P_3x      = DSG_3v3_2nd_Q_P[   grepl("3x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_3x_iPSC = DSG_3v3_2nd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_3x) ) ]
DSG_3v3_2nd_Q_P_3x_iPSC
colnames(DSG_3v3_2nd_Q_P_3x_iPSC)
head(DSG_3v3_2nd_Q_P_3x_iPSC)
head(iPSC_data)

DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs ,
                                DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs , 
                                DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_2nd_Q_P_3x_iPSC)
fout = "DSG_3v3_2nd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_3x_iPSC.pdf")

#### (1.2.2) 10x iPSC
DSG_3v3_2nd_Q_P_10x      = DSG_3v3_2nd_Q_P[   grepl("10x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_10x_iPSC = DSG_3v3_2nd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_10x) ) ]
DSG_3v3_2nd_Q_P_10x_iPSC

DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs ,
                                 DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs , 
                                 DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_10x_iPSC.pdf")

################# (2) NGN2 ==========================================================================
###== (2.1) dPSI === ####
# colnames(DSG_3v3_2nd  )
# head(DSG_3v3_2nd, 2  )
# DSG_3v3_2nd_dPSI  = DSG_3v3_2nd[ , grep("^dPSI_",   colnames(DSG_3v3_2nd  ) ) ] 
# DSG_3v3_2nd_dPSI[1:10,]
# colnames(DSG_3v3_2nd_dPSI)

## === (2.1.1) 3x NGN2
DSG_3v3_2nd_dPSI_3x      = DSG_3v3_2nd_dPSI[   grepl("3x",   colnames(DSG_3v3_2nd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_dPSI_3x_NGN2 = DSG_3v3_2nd_dPSI_3x[grepl("NGN2", colnames(DSG_3v3_2nd_dPSI_3x) ) ]
DSG_3v3_2nd_dPSI_3x_NGN2

DSG_3v3_2nd_dPSI_3x_NGN2_4088_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`4088`, "dPSI_NGN2.4088.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_4090_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`4090`, "dPSI_NGN2.4090.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_4714_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`4714`, "dPSI_NGN2.4714.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_4741_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`4741`, "dPSI_NGN2.4741.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_4748_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`4748`, "dPSI_NGN2.4748.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_5420_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`5420`, "dPSI_NGN2.5420.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_6152_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`6152`, "dPSI_NGN2.6152.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_6527_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`6527`, "dPSI_NGN2.6527.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_6866_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`6866`, "dPSI_NGN2.6866.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_6960_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`6960`, "dPSI_NGN2.6960.3x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_NGN2_9701_abs = abs(DSG_3v3_2nd_dPSI_3x_NGN2[NGN2_data$`9701`, "dPSI_NGN2.9701.3x", drop=FALSE ] )

DSG_3v3_2nd_dPSI_3x_NGN2_abs = cbind(DSG_3v3_2nd_dPSI_3x_NGN2_4088_abs , DSG_3v3_2nd_dPSI_3x_NGN2_4090_abs , DSG_3v3_2nd_dPSI_3x_NGN2_4714_abs , DSG_3v3_2nd_dPSI_3x_NGN2_4741_abs ,
                                     DSG_3v3_2nd_dPSI_3x_NGN2_4748_abs , DSG_3v3_2nd_dPSI_3x_NGN2_5420_abs , DSG_3v3_2nd_dPSI_3x_NGN2_6152_abs , DSG_3v3_2nd_dPSI_3x_NGN2_6527_abs , 
                                     DSG_3v3_2nd_dPSI_3x_NGN2_6866_abs , DSG_3v3_2nd_dPSI_3x_NGN2_6960_abs , DSG_3v3_2nd_dPSI_3x_NGN2_9701_abs )
fout = "DSG_3v3_2nd_dPSI_3x_NGN2_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_3x_NGN2_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_3x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_3x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_3x_NGN2_abs.pdf")

#### (2.1.2) 10x NGN2
DSG_3v3_2nd_dPSI_10x      = DSG_3v3_2nd_dPSI[   grepl("10x",   colnames(DSG_3v3_2nd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_dPSI_10x_NGN2 = DSG_3v3_2nd_dPSI_10x[grepl("NGN2", colnames(DSG_3v3_2nd_dPSI_10x) ) ]
DSG_3v3_2nd_dPSI_10x_NGN2

DSG_3v3_2nd_dPSI_10x_NGN2_4088_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`4088`, "dPSI_NGN2.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_4090_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`4090`, "dPSI_NGN2.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_4714_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`4714`, "dPSI_NGN2.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_4741_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`4741`, "dPSI_NGN2.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_4748_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`4748`, "dPSI_NGN2.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_5420_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`5420`, "dPSI_NGN2.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_6152_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`6152`, "dPSI_NGN2.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_6527_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`6527`, "dPSI_NGN2.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_6866_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`6866`, "dPSI_NGN2.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_6960_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`6960`, "dPSI_NGN2.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_NGN2_9701_abs = abs(DSG_3v3_2nd_dPSI_10x_NGN2[NGN2_data$`9701`, "dPSI_NGN2.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_dPSI_10x_NGN2_abs = cbind(DSG_3v3_2nd_dPSI_10x_NGN2_4088_abs , DSG_3v3_2nd_dPSI_10x_NGN2_4090_abs , DSG_3v3_2nd_dPSI_10x_NGN2_4714_abs , DSG_3v3_2nd_dPSI_10x_NGN2_4741_abs ,
                                      DSG_3v3_2nd_dPSI_10x_NGN2_4748_abs , DSG_3v3_2nd_dPSI_10x_NGN2_5420_abs , DSG_3v3_2nd_dPSI_10x_NGN2_6152_abs , DSG_3v3_2nd_dPSI_10x_NGN2_6527_abs , 
                                      DSG_3v3_2nd_dPSI_10x_NGN2_6866_abs , DSG_3v3_2nd_dPSI_10x_NGN2_6960_abs , DSG_3v3_2nd_dPSI_10x_NGN2_9701_abs )
fout = "DSG_3v3_2nd_dPSI_10x_NGN2_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_10x_NGN2_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_10x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_10x_NGN2_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_10x_NGN2_abs.pdf")

###== (2.2) significance === ####
### 2.2 Q and Prob
DSG_3v3_2nd_Q_P  = DSG_3v3_2nd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_2nd  ) ) ] 
DSG_3v3_2nd_Q_P[1:10,]
colnames(DSG_3v3_2nd_Q_P)
head(DSG_3v3_2nd_Q_P)

## === (2.2.1) 3x iPSC
DSG_3v3_2nd_Q_P_3x      = DSG_3v3_2nd_Q_P[   grepl("3x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_3x_iPSC = DSG_3v3_2nd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_3x) ) ]
DSG_3v3_2nd_Q_P_3x_iPSC
colnames(DSG_3v3_2nd_Q_P_3x_iPSC)
head(DSG_3v3_2nd_Q_P_3x_iPSC)
head(NGN2_data)

DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[NGN2_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs ,
                                    DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs , 
                                    DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_2nd_Q_P_3x_iPSC)
fout = "DSG_3v3_2nd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_3x_iPSC.pdf")

#### (2.2.2) 10x iPSC
DSG_3v3_2nd_Q_P_10x      = DSG_3v3_2nd_Q_P[   grepl("10x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_10x_iPSC = DSG_3v3_2nd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_10x) ) ]
DSG_3v3_2nd_Q_P_10x_iPSC

DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[NGN2_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs ,
                                     DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs , 
                                     DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_10x_iPSC.pdf")


################# (3) 6152 ==========================================================================
dim(Data_6152)
head(Data_6152)
tail(Data_6152)
colnames(Data_6152) # [1] "iPSC_3x_1st"                  "iPSC_3x_2nd"                  "iPSC_3x_supported_by_10x_1st" "iPSC_3x_supported_by_10x_2nd" "NGN2_3x_1st"                 
#[6] "NGN2_3x_2nd"                  "NGN2_3x_supported_by_10x_1st" "NGN2_3x_supported_by_10x_2nd"
Data_6152[,1] # 41] " "          NA     

###== (3.1) dPSI === ####
DSG_3v3_1st_6152  = DSG_3v3_1st[ , grep("6152",   colnames(DSG_3v3_1st  ) ) ] 
DSG_3v3_2nd_6152  = DSG_3v3_2nd[ , grep("6152",   colnames(DSG_3v3_2nd  ) ) ] 
DSG_3v3_1st_6152
colnames(DSG_3v3_1st_6152)
DSG_3v3_2nd_6152
colnames(DSG_3v3_2nd_6152)

## === (3.1.1) 3x 6152


DSG_3v3_1st_dPSI_3x_iPSC_6152_abs = abs(DSG_3v3_1st_6152[Data_6152$iPSC_3x_1st, "dPSI_TST11872_iPSC.BIO.2006152.45nM", drop=FALSE ] )
DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs = abs(DSG_3v3_2nd_6152[Data_6152$iPSC_3x_2nd, "dPSI_iPSC.6152.3x", drop=FALSE ] )
DSG_3v3_1st_dPSI_3x_iPSC_6152_abs
DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs

DSG_3v3_1st_dPSI_10x_iPSC_6152_abs = abs(DSG_3v3_1st_6152[Data_6152$iPSC_10x_1st, "dPSI_TST11872_iPSC.BIO.2006152.150nM", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_6152[Data_6152$iPSC_10x_2nd, "dPSI_iPSC.6152.10x", drop=FALSE ] )

DSG_3v3_1st_dPSI_10x_iPSC_6152_abs
DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs


DSG_3v3_dPSI_iPSC_6152_abs = cbind( )

fout = "DSG_3v3_2nd_dPSI_6152_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_3x_iPSC_6152_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_6152_abs.pdf")

#### (3.1.2) 10x iPSC
DSG_3v3_2nd_dPSI_10x      = DSG_3v3_2nd_dPSI[   grepl("10x",   colnames(DSG_3v3_2nd_dPSI) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_dPSI_10x_iPSC = DSG_3v3_2nd_dPSI_10x[grepl("iPSC", colnames(DSG_3v3_2nd_dPSI_10x) ) ]
DSG_3v3_2nd_dPSI_10x_iPSC

DSG_3v3_2nd_dPSI_10x_iPSC_4088_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4088`, "dPSI_iPSC.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4090_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4090`, "dPSI_iPSC.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4714_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4714`, "dPSI_iPSC.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4741_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4741`, "dPSI_iPSC.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_4748_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`4748`, "dPSI_iPSC.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_5420_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`5420`, "dPSI_iPSC.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6152`, "dPSI_iPSC.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6527_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6527`, "dPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6866_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6866`, "dPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_6960_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`6960`, "dPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_dPSI_10x_iPSC_9701_abs = abs(DSG_3v3_2nd_dPSI_10x_iPSC[iPSC_data$`9701`, "dPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_dPSI_10x_iPSC_abs = cbind(DSG_3v3_2nd_dPSI_10x_iPSC_4088_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4090_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4714_abs , DSG_3v3_2nd_dPSI_10x_iPSC_4741_abs ,
                                      DSG_3v3_2nd_dPSI_10x_iPSC_4748_abs , DSG_3v3_2nd_dPSI_10x_iPSC_5420_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6152_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6527_abs , 
                                      DSG_3v3_2nd_dPSI_10x_iPSC_6866_abs , DSG_3v3_2nd_dPSI_10x_iPSC_6960_abs , DSG_3v3_2nd_dPSI_10x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_dPSI_10x_iPSC_abs.csv"
write.csv(DSG_3v3_2nd_dPSI_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_dPSI_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_dPSI_10x_iPSC_abs.pdf")

###== (3.2) significance === ####
DSG_3v3_2nd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_2nd)
DSG_3v3_2nd[1:10, 1:10]
DSG_3v3_2nd[1:10, ]
colnames(DSG_3v3_2nd)
dim(DSG_3v3_2nd)
row.names(DSG_3v3_2nd) = DSG_3v3_2nd$X
DSG_3v3_2nd$X = NULL

DSG_3v3_2nd
colnames(DSG_3v3_2nd)
head(DSG_3v3_2nd)

### 1.2 Q and Prob
DSG_3v3_2nd_Q_P  = DSG_3v3_2nd[ , grep("^minpadj_maxP",   colnames(DSG_3v3_2nd  ) ) ] 
DSG_3v3_2nd_Q_P[1:10,]
colnames(DSG_3v3_2nd_Q_P)
head(DSG_3v3_2nd_Q_P)

## === (3.2.1) 3x iPSC
DSG_3v3_2nd_Q_P_3x      = DSG_3v3_2nd_Q_P[   grepl("3x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_3x_iPSC = DSG_3v3_2nd_Q_P_3x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_3x) ) ]
DSG_3v3_2nd_Q_P_3x_iPSC
colnames(DSG_3v3_2nd_Q_P_3x_iPSC)
head(DSG_3v3_2nd_Q_P_3x_iPSC)
head(iPSC_data)

DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.3x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_3x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.3x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_3x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_3x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_3x_iPSC_4741_abs ,
                                    DSG_3v3_2nd_Q_P_3x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_3x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6527_abs , 
                                    DSG_3v3_2nd_Q_P_3x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_3x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_3x_iPSC_9701_abs )
head(DSG_3v3_2nd_Q_P_3x_iPSC)
fout = "DSG_3v3_2nd_Q_P_3x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_3x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_3x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_3x_iPSC.pdf")

#### (3.2.2) 10x iPSC
DSG_3v3_2nd_Q_P_10x      = DSG_3v3_2nd_Q_P[   grepl("10x",   colnames(DSG_3v3_2nd_Q_P) ) ] # https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
DSG_3v3_2nd_Q_P_10x_iPSC = DSG_3v3_2nd_Q_P_10x[grepl("iPSC", colnames(DSG_3v3_2nd_Q_P_10x) ) ]
DSG_3v3_2nd_Q_P_10x_iPSC

DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4088`, "minpadj_maxPdPSI_iPSC.4088.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4090`, "minpadj_maxPdPSI_iPSC.4090.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4714`, "minpadj_maxPdPSI_iPSC.4714.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4741`, "minpadj_maxPdPSI_iPSC.4741.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`4748`, "minpadj_maxPdPSI_iPSC.4748.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`5420`, "minpadj_maxPdPSI_iPSC.5420.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6152`, "minpadj_maxPdPSI_iPSC.6152.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6527`, "minpadj_maxPdPSI_iPSC.6527.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6866`, "minpadj_maxPdPSI_iPSC.6866.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`6960`, "minpadj_maxPdPSI_iPSC.6960.10x", drop=FALSE ] )
DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs = abs(DSG_3v3_2nd_Q_P_10x_iPSC[iPSC_data$`9701`, "minpadj_maxPdPSI_iPSC.9701.10x", drop=FALSE ] )

DSG_3v3_2nd_Q_P_10x_iPSC_abs = cbind(DSG_3v3_2nd_Q_P_10x_iPSC_4088_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4090_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4714_abs , DSG_3v3_2nd_Q_P_10x_iPSC_4741_abs ,
                                     DSG_3v3_2nd_Q_P_10x_iPSC_4748_abs , DSG_3v3_2nd_Q_P_10x_iPSC_5420_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6152_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6527_abs , 
                                     DSG_3v3_2nd_Q_P_10x_iPSC_6866_abs , DSG_3v3_2nd_Q_P_10x_iPSC_6960_abs , DSG_3v3_2nd_Q_P_10x_iPSC_9701_abs )

fout = "DSG_3v3_2nd_Q_P_10x_iPSC.csv"
write.csv(DSG_3v3_2nd_Q_P_10x_iPSC_abs,file = fout, quote = F, row.names = F)

pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
xx = pheatmap(DSG_3v3_2nd_Q_P_10x_iPSC_abs ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
save_pheatmap_pdf(xx, "DSG_3v3_2nd_Q_P_10x_iPSC.pdf")



#final[rowSums(is.na(final[ , 5:6])) == 0, ]
# (base) /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code zgao1@camhpcps01 $ 
#   $  tree | grep heatmap
# │   ├── run.02.DSG_counts_minmaxtable_07-26-merged-df-for-heatmap.R
# │   ├── run.02.DSG_counts_minmaxtable_plot_2022-07-25-for-correlation-heatmap.R
# │   ├── pheatmap_comparions_1.pdf
# ├── Miyoung_Shin_heatmap.R
# ├── run.02-1.DSG_counts_plot_ZG_2022-07-26-similarity-heatmap.R
# ├── run.02-1.DSG_counts_plot_ZG_2022-07-27-similarity-heatmap-CAM.R

### 6152

DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/summarytable_with_gene-1st_experiment.csv", header = T)
dim(DSG_3v3_1st)
DSG_3v3_1st[1:10, 1:10]
DSG_3v3_1st[1:10, ]
# rownames(minmaxall) = minmaxall$genesymbol;   minmaxall$genesymbol = NULL 
# df <- df %>% dplyr:: select(starts_with("ABC"))
# df[, grep("ABC", names(df)), with = FALSE]
row.names(DSG_3v3_1st) = DSG_3v3_1st$X
DSG_3v3_1st$X = NULL

DSG_3v3_1st  = DSG_3v3_1st[ , grep("^dPSI_",   colnames(DSG_3v3_1st  ) ) ] 
colnames(DSG_3v3_1st) = gsub("dPSI_TST11872","1st_exp", colnames(DSG_3v3_1st)) 
colnames(DSG_3v3_1st) = gsub(".BIO",         "", colnames(DSG_3v3_1st)) 
DSG_3v3_1st[1:10, 1:10]

DSG_3v3_1st_6152           = DSG_3v3_1st[      , grep("2006152", colnames(DSG_3v3_1st      ) ) ] 
colnames(DSG_3v3_1st_6152) = gsub("150nM",   "10x" , colnames(DSG_3v3_1st_6152)) 
colnames(DSG_3v3_1st_6152) = gsub("45nM",    "3x"  , colnames(DSG_3v3_1st_6152)) 
colnames(DSG_3v3_1st_6152) = gsub("200",     ""    , colnames(DSG_3v3_1st_6152)) 
dim(DSG_3v3_1st_6152)
head(DSG_3v3_1st_6152)

### 2nd exp
DSG_3v3_2nd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
dim(DSG_3v3_2nd)
DSG_3v3_2nd[1:10, 1:10]
DSG_3v3_2nd[1:10, ]

row.names(DSG_3v3_2nd) = DSG_3v3_2nd$X
DSG_3v3_2nd$X = NULL

DSG_3v3_2nd  = DSG_3v3_2nd[ , grep("^dPSI_",   colnames(DSG_3v3_2nd  ) ) ] 
colnames(DSG_3v3_2nd) = gsub("dPSI","2nd_exp", colnames(DSG_3v3_2nd)) #colnames(DSG_3v3_2nd) = gsub(".BIO",         "", colnames(DSG_3v3_2nd)) 
DSG_3v3_2nd[1:10, 1:10]

DSG_3v3_2nd_6152           = DSG_3v3_2nd[      , grep("6152", colnames(DSG_3v3_2nd      ) ) ] 
dim(DSG_3v3_2nd_6152)
head(DSG_3v3_2nd_6152)

DSG_3v3_6152 = merge(DSG_3v3_1st_6152, DSG_3v3_2nd_6152, by = 0, all = T)
head(DSG_3v3_6152)
tail(DSG_3v3_6152)
DSG_3v3_6152[ 500:550, ]
dim(DSG_3v3_6152) # 3927
class(DSG_3v3_6152)

row.names(DSG_3v3_6152) =  DSG_3v3_6152$Row.names
DSG_3v3_6152$Row.names  = NULL

dim(DSG_3v3_6152)

DSG_3v3_6152[ is.na(DSG_3v3_6152) ]  = 0 # df[rowSums(is.na(df)) == 0, ]
DSG_3v3_6152        =  DSG_3v3_6152[ , c("1st_exp_iPSC.6152.10x", "2nd_exp_iPSC.6152.10x", "1st_exp_iPSC.6152.3x", "2nd_exp_iPSC.6152.3x", "1st_exp_NGN2.6152.10x", "2nd_exp_NGN2.6152.10x", "1st_exp_NGN2.6152.3x", "2nd_exp_NGN2.6152.3x")]
DSG_3v3_6152.abs.mx =  abs(as.matrix(DSG_3v3_6152))
DSG_3v3_6152        =  as.data.frame(DSG_3v3_6152.abs.mx)
DSG_3v3_6152        = DSG_3v3_6152[order(DSG_3v3_6152$`1st_exp_iPSC.6152.10x`,decreasing = T) ,]
DSG_3v3_6152 = DSG_3v3_6152[ rowSums(DSG_3v3_6152[]) > 0.1,]

head(DSG_3v3_6152)
tail(DSG_3v3_6152)
dim(DSG_3v3_6152)

getwd()
fout = "compound_6152_induced_DSG_in_1st_and_2nd_experiments1.csv"
write.csv(DSG_3v3_6152,file = fout,quote = F, row.names = T)

library(pheatmap)
#pheatmap(DSG_3v3_6152 ,fontsize=9, fontsize_row=6, cluster_rows = T,         cluster_cols = F) 
pheatmap(DSG_3v3_6152 ,fontsize=9, fontsize_row=6, cluster_rows = F,         cluster_cols = F) 
#pheatmap(DSG_3v3_6152 ,fontsize=9, fontsize_row=6) 

# Data For Magnus
# run.04.DEG_counts_plot_ranking_batch123_edge7-3x-validated2-with_UnigquID.R
#install.packages("patchwork") # require(devtools) # install_version("rowr", version = "1.1.3", repos = "http://cran.us.r-project.org")
R.Version()
rm(list = ls())        #library(gplots)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(withr);          # in the Packages withr, select the 2.4.2 version  shown as library(withr, lib.loc = "/opt/R/4.0.3/lib/R/library")
library(rowr)
library(patchwork)         #https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors
Date           = Sys.Date()
Date           = gsub("-", "_", Date) # to avoid java commandline parse error
Date

dout       = paste0("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/" ,  Date, "/DEG/")
dir.create(dout, recursive = T)
dir.create(paste0(dout, "/Magnus"))
dir.create(paste0(dout, "/Jessica"))
setwd(dout); 

getwd() #setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/2023_07_27/DEG/"); 

# ##############################   Experiment 111111----------------------------------------- iPSC Note ======PTC518 compound is BIO-2197294 
# ##############################   Experiment1  iPSC 3x  ############################## ############################## ##############################  
din1 = "/edgehpc/dept/compbio/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din1)
fin1s = list.files( path = din1 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin1s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursivelyfin1s
fin1s = fin1s[! fin1s %in% c("iPSC-1949634-3xvsiPSC-DMSO.tsv", "iPSC-2060884-3xvsiPSC-DMSO.tsv") ] # Remove 2 samples habe more than 1 thousand reduced genes
fin1s
str(fin1s)

Names = gsub( "vsiPSC-DMSO.tsv" ,  "Batch1" , fin1s)
Names
Names
str(Names)

############################### (1.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin1s)     

datalist <- setNames(
  lapply(fin1s, function(i){
    read.table(i, header=TRUE) }) ,
  paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
  ) 

str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (1.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (1.3) ==> merge the DataSet to A Single Table Magnus ############################### 
# i=1
# i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "UniqueID", "Gene.Name", "log2FoldChange" )]  
  # tmp 
  # colnames(tmp)  = c( paste0( "Index", i, "_",  Names[i], "_",  "Gene.Name"),   "log2FoldChange"  )
  #tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "UniqueID",  "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" )  )
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = rowr::cbind.fill(temp_count_table,  tmp, fill=NA)  
  }
  print(head(temp_count_table ))
  
  
}

getwd()
fout = paste0( dout,"/Magnus" , "/TST11872_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)

###############################   ========================  Experiment 22222222222222 ############################## ############################## ##############################  
din2 = "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din2)
fin2s = list.files( "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" , pattern = "3xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin2s
str(fin2s)

Names = gsub( "vsiPSC_DMSO.tsv" ,  "Batch2" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Names = gsub( "vsiPSC"-"DMSO.tsv" ,  "" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Names
Names
str(Names)
#i=1
i
Names[i]
############################### (2.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin2s)                       #  paste0("iPSC_10x", s))
datalist <- setNames(
  lapply(fin2s, function(i){
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
) 
str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (2.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (2.3) ==> merge the DataSet to A Single Table Magnus ############################### 
# i=1
# i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "UniqueID", "Gene.Name", "log2FoldChange" )]  
  # tmp 
  # colnames(tmp)  = c( paste0( "Index", i, "_",  Names[i], "_",  "Gene.Name"),   "log2FoldChange"  )
  #tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "UniqueID",  "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" )  )
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = rowr::cbind.fill(temp_count_table,  tmp, fill=NA)  
  }
  print(head(temp_count_table ))
  
}  

getwd() #fout = paste0("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/2023_07_27/DEG/TST11872_reduced_EXpression_Jessica",".csv") 
fout = paste0( dout,"/Magnus" , "/TST11955_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)


###############################   Experiment3 ############################## ############################## ############################## 
din3 = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din3)

#fin3s = list.files( path = din1 ) # list.dirs all sub_dir recursively
fin3s = list.files( path = din3 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin3s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin3s
str(fin3s)

Names = gsub( "vsiPSC-DMSO.tsv" ,  "Batch3" , fin3s)
Names
Names
str(Names)
# i=1
#i
Names[i]
############################### (3.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin3s)                       #  paste0("iPSC_10x", s))
datalist <- setNames(
  lapply(fin3s, function(i) { ####### 1s 2s 3s
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
) 

str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (3.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (3.3) ==> merge the DataSet to A Single Table Magnus ############################### 
i=1
i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "UniqueID", "Gene.Name", "log2FoldChange" )]  
  # tmp 
  # colnames(tmp)  = c( paste0( "Index", i, "_",  Names[i], "_",  "Gene.Name"),   "log2FoldChange"  )
  #tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "UniqueID",  "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" )  )
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = rowr::cbind.fill(temp_count_table,  tmp, fill=NA)  
  }
  print(head(temp_count_table ))
  
}

getwd()
fout = paste0( dout,"/Magnus" , "/TST12086_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)

# for Jessica
# run.04.DEG_counts_plot_ranking_batch123_1edge-3x-good-and_merge1.R

#install.packages("patchwork") # require(devtools) # install_version("rowr", version = "1.1.3", repos = "http://cran.us.r-project.org")
R.Version()
rm(list = ls())        #library(gplots)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(withr);          # in the Packages withr, select the 2.4.2 version  shown as library(withr, lib.loc = "/opt/R/4.0.3/lib/R/library")
library(rowr)
library(patchwork)         #https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors
Date           = Sys.Date()
Date           = gsub("-", "_", Date) # to avoid java commandline parse error
Date

dout       = paste0("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/" ,  Date, "/DEG/")
dir.create(dout, recursive = T)
dir.create(paste0(dout, "/Magnus"))
dir.create(paste0(dout, "/Jessica"))
setwd(dout); 

getwd() #setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/2023_07_27/DEG/"); 

# ##############################   Experiment 111111----------------------------------------- iPSC Note ======PTC518 compound is BIO-2197294 
# ##############################   Experiment1  iPSC 3x  ############################## ############################## ##############################  
din1 = "/edgehpc/dept/compbio/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din1)
fin1s = list.files( path = din1 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin1s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursivelyfin1s
fin1s = fin1s[! fin1s %in% c("iPSC-1949634-3xvsiPSC-DMSO.tsv", "iPSC-2060884-3xvsiPSC-DMSO.tsv") ] # Remove 2 samples habe more than 1 thousand reduced genes
fin1s
str(fin1s)

Names = gsub( "vsiPSC-DMSO.tsv" ,  "Batch1" , fin1s)
Names
Names
str(Names)

############################### (1.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin1s)     

datalist <- setNames(
  lapply(fin1s, function(i){
    read.table(i, header=TRUE) }) ,
  paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
  ) 

str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (1.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (1.3) ==> merge the DataSet to A Single Table Magnus ############################### 
# i=1
# i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "Gene.Name", "log2FoldChange", "padj" )]  
  tmp 
  #row.names(tmp) = make.names( tmp$Gene.Name  , unique=TRUE)
  tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" ) ,  paste0( "Index", i, "_",  "padj"  ) )

  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = merge(temp_count_table,  tmp, by = "Gene.Name", all = T)  
  }
  print( head(temp_count_table ) )
}

getwd()
fout = paste0( dout,"/Jessica" , "/TST11872_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)

###############################   ========================  Experiment 22222222222222 ############################## ############################## ##############################  
din2 = "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din2)
fin2s = list.files( "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" , pattern = "3xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin2s
str(fin2s)

Names = gsub( "vsiPSC_DMSO.tsv" ,  "Batch2" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Names = gsub( "vsiPSC"-"DMSO.tsv" ,  "" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Names
Names
str(Names)
#i=1
i
Names[i]
############################### (2.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin2s)                       #  paste0("iPSC_10x", s))
datalist <- setNames(
  lapply(fin2s, function(i){
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
) 
str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (2.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (2.3) ==> merge the DataSet to A Single Table Magnus ############################### 
# i=1
# i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "Gene.Name", "log2FoldChange", "padj" )]  
  tmp 
  #row.names(tmp) = make.names( tmp$Gene.Name  , unique=TRUE)
  tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" ) ,  paste0( "Index", i, "_",  "padj"  ) )
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = merge(temp_count_table,  tmp, by = "Gene.Name", all = T)  
  }
  
  print(head(temp_count_table ))
  (temp_count_table )
}

getwd() #fout = paste0("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/2023_07_27/DEG/TST11872_reduced_EXpression_Jessica",".csv") 
fout = paste0( dout,"/Jessica" , "/TST11955_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)


###############################   Experiment3 ############################## ############################## ############################## 
din3 = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
setwd(din3)

#fin3s = list.files( path = din1 ) # list.dirs all sub_dir recursively
fin3s = list.files( path = din3 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin3s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin3s
str(fin3s)

Names = gsub( "vsiPSC-DMSO.tsv" ,  "Batch3" , fin3s)
Names
Names
str(Names)
# i=1
i
Names[i]
############################### (3.1) ==> Prepare The DataSet ###############################
s <- 1:length(fin3s)                       #  paste0("iPSC_10x", s))
datalist <- setNames(
  lapply(fin3s, function(i) { ####### 1s 2s 3s
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s, "_",  Names[s] ) #    ) , paste0("iPSC_3x", s) =================
) 

str(datalist)
head(datalist[[1]])
head(datalist[[3]])
dim(datalist[[3]])

############################### (3.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist)) {
  temp = datalist[[i]] 
  print( names(datalist)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (3.3) ==> merge the DataSet to A Single Table Magnus ############################### 
i=1
i=i+1
i

for(i in 1:length(datalist)) {
  tmp = datalist[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c( "Gene.Name", "log2FoldChange", "padj" )]  
  tmp 
  #row.names(tmp) = make.names( tmp$Gene.Name  , unique=TRUE)
  tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  colnames(tmp)  = c( "Gene.Name", paste0( "Index", i, "_", Names[i]  ,"_", "log2FoldChange" ) ,  paste0( "Index", i, "_",  "padj"  ) )
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = merge(temp_count_table,  tmp, by = "Gene.Name", all = T)  
  }
  
  print(head(temp_count_table ))
  (temp_count_table )
  
}

getwd()
fout = paste0( dout,"/Jessica" , "/TST12086_reduced_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)


############################### (4)  DataSet combined from DataSet 1, 2, 3############################### ############################### (4.1) ==> Prepare The DataSet ###############################
############################### (4.1) ==> Prepare The DataSet ###############################
din1 = "/edgehpc/dept/compbio/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#setwd(din1)
fin1s = list.files( path = din1 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin1s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursivelyfin1s
fin1s = fin1s[! fin1s %in% c("iPSC-1949634-3xvsiPSC-DMSO.tsv", "iPSC-2060884-3xvsiPSC-DMSO.tsv") ] # Remove 2 samples habe more than 1 thousand reduced genes
fin1s
str(fin1s)

Names_Batch1 = gsub( "vsiPSC-DMSO.tsv" ,  "Batch1" , fin1s)
Names_Batch1
Names_Batch1


din2 = "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#setwd(din2)
fin2s = list.files( "/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" , pattern = "3xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin2s
str(fin2s)

Names_Batch2 = gsub( "vsiPSC_DMSO.tsv" ,  "Batch2" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Names = gsub( "vsiPSC"-"DMSO.tsv" ,  "" , fin2s) ########## ====================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Names_Batch2
Names_Batch2
str(Names_Batch2)
#i=1
Names_Batch2[1]

din3 = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/EA_DEG_automatic_files/Omics_DEG/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#setwd(din3)
#fin3s = list.files( path = din1 ) # list.dirs all sub_dir recursively
fin3s = list.files( path = din3 , pattern = "3xvsiPSC-DMSO.tsv$") # list.dirs all sub_dir recursively #fin3s = list.files( path=din1 , pattern = "10xvsiPSC_DMSO.tsv$") # list.dirs all sub_dir recursively
fin3s
str(fin3s)

Names_Batch3 = gsub( "vsiPSC-DMSO.tsv" ,  "Batch3" , fin3s)
Names_Batch3
Names_Batch3
str(Names_Batch3)

s1 <- 1:length(fin1s)                       #  paste0("iPSC_10x", s))
s2 <- 1:length(fin2s)                       #  paste0("iPSC_10x", s))
s3 <- 1:length(fin3s)                       #  paste0("iPSC_10x", s))
s1 
s2 
s3 

datalist_Batch1 <- setNames(
  lapply( paste0(din1, fin1s) , function(i) { ####### 1s 2s 3s
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s1, "_",  Names_Batch1[s1] ) #    ) , paste0("iPSC_3x", s) =================
) 

datalist_Batch2 <- setNames(
  lapply(paste0(din2, fin2s), function(i) { ####### 1s 2s 3s
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s2, "_",  Names_Batch2[s2] ) #    ) , paste0("iPSC_3x", s) =================
) 

datalist_Batch3 <- setNames(
  lapply( paste0(din3,  fin3s) , function(i) { ####### 1s 2s 3s
    read.table(i, header=TRUE) }
  ) , paste0( "Index", s3, "_",  Names_Batch3[s3] ) #    ) , paste0("iPSC_3x", s) =================
) 



datalist_Batch123  = c(datalist_Batch1, datalist_Batch2, datalist_Batch3) 

str(datalist_Batch1)
str(datalist_Batch2)
str(datalist_Batch3)
str(datalist_Batch123)

head(datalist_Batch1[[1]])
head(datalist_Batch2[[3]])
head(datalist_Batch3[[3]])
dim(datalist_Batch3[[3]])

length(datalist_Batch123)

############################### (4.2) ==> show the Dimension of the DataSet ###############################
for(i in 1:length(datalist_Batch123)) {
  temp = datalist_Batch123[[i]] 
  print( names(datalist_Batch123)[i]  ) # names(my_list)[2]  # https://statisticsglobe.com/extract-names-of-list-elements-in-r
  #print(dim(temp)) 
  temp1 = temp[which(temp$padj < 0.05 & temp$log2FoldChange < -0.585 ), ] # 0.585 => 0.5, 0.2 => 0.263 choose Reduction Only here  for reduction only in earlier ploting, HERE I choose 2023-07-01 to keep
  print(dim(temp1))
}

############################### (4.3) ==> merge the DataSet to A Single Table Magnus ############################### 
# i=1
# i=i+1
i

for(i in 1:length(datalist_Batch123)) {
  tmp = datalist_Batch123[[i]]
  tmp = tmp[which(tmp$log2FoldChange< -0.585 & tmp$padj< 0.05 ), c("UniqueID", "Gene.Name", "log2FoldChange", "padj" )]  
  tmp 
  #row.names(tmp) = make.names( tmp$Gene.Name  , unique=TRUE)
  #tmp$Gene.Name = make.names( tmp$Gene.Name  , unique=TRUE)
  
  colnames(tmp)  = c( "UniqueID", "Gene.Name", paste0( names(datalist_Batch123)[i]   ,"_", "log2FoldChange" ) ,   "padj"  ) 
  
  if (i == 1){
    temp_count_table = tmp 
  } else {
    temp_count_table = merge(temp_count_table,  tmp, by = "UniqueID", all = T)  #     temp_count_table = merge(temp_count_table,  tmp, by = "Gene.Name", all = T)  
  }
  
  print(head(temp_count_table ))
  (temp_count_table )
  
}



detected_genes = temp_count_table$UniqueID # detected_genes = temp_count_table$Gene.Name

#i=1
# i=i+1
i

#"UniqueID"
for(i in 1:length(datalist_Batch123)) {
  tmp_full = datalist_Batch123[[i]]
  tmp_full = tmp_full[ which(tmp_full$UniqueID %in% detected_genes),  c( "UniqueID", "Gene.Name", "log2FoldChange",   "padj" )  ]  #   tmp_full = tmp_full[ which(tmp_full$Gene.Name %in% detected_genes), ]  
  
  tmp_full 
  #row.names(tmp_full) = make.names( tmp_full$Gene.Name  , unique=TRUE)
  #tmp_full$Gene.Name = make.names( tmp_full$Gene.Name  , unique=TRUE)
  
  colnames(tmp_full)  = c( "UniqueID", "Gene.Name", paste0( names(datalist_Batch123)[i]   ,"_", "log2FoldChange" ) ,   "padj"  ) 
  
  if (i == 1){
    temp_count_table = tmp_full 
  } else {
    temp_count_table = merge(temp_count_table,  tmp_full, by = "UniqueID", all = T)  
  }
  
  print(head(temp_count_table ))
  (temp_count_table )
  
}

getwd()
fout = paste0( dout,"/Jessica" , "/TST11872_11955_12086_reduced_merged_EXpression",".csv")
write.csv(temp_count_table, file=fout,  row.names=F)

