R.Version()
rm(list = ls())
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(withr); # in the Packages withr, select the 2.4.2 version  shown as library(withr, lib.loc = "/opt/R/4.0.3/lib/R/library") #install.packages("patchwork")
library(patchwork) #https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r
library(dplyr) #outer_join  #library(plyr) #join # 
library(forcats)  # fct_relevel
library(stringr)
library(VennDiagram)


pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors

# dout = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis/" # din = "./res.02.DSG_counts/" 

dout = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig/Magnus_05_11_request/" # din = "./res.02.DSG_counts/" 
dir.create(dout, recursive = T)
setwd(dout); 
getwd()                  # [1] "/mnt/depts/dept04/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig"

# ##############################   111111----------------------------------------- iPSC Note ======PTC518 compound is BIO-2197294
##############################=========      2023-04-18 ===========================
# 1) Ven diagrams of BIO-2184090, BIO-2195127 BIO-2186960 à what is the percent overlap of off-targets?
# 2) Compound 6152 was run in all 3 batches, create list of off-targets for iPSC for each run. How much overlap between runs?
##############################=========      2023-04-25 ===========================
# =============================================== PTC518 comp_3_title     = "PTC-518-3x: iPSC vs NGN2 vs Sy5Y"
# =============================================== # 1755497 is Branaplam, 1949634 Risdiplam, 2186960
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

# 1755497 is ============>> Branaplam,            05-16========                         1949634 ===========>>> Risdiplam, 206088s Interal Lead
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
param = gsub(".csv","",gsub("binary_","",fin2))
dPSI_bi = read.csv( paste0(din2, fin2)  , stringsAsFactors = F , row.names = 1 , check.names = F) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp            =  dPSI_bi[     , grep('iPSC-',  colnames(dPSI_bi) ) ] 
#dPSI_bi_tmp            =  dPSI_bi[     , grep('NGN2-',  colnames(dPSI_bi) ) ] 

dPSI_bi_tmp            =  dPSI_bi_tmp[ , grep('-3x'  ,  colnames( dPSI_bi_tmp) )] #

head(dPSI_bi_tmp) # colnames(dPSI_bi_tmp)  = gsub("BIO-", "", colnames(dPSI_bi_tmp ) )
colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c("2184088_iPSC_3x", "2184090_iPSC_3x", "2174714_iPSC_3x", "2184741_iPSC_3x", "2174748_iPSC_3x", "2175420_iPSC_3x", 
                           "2006152_iPSC_3x", "2186527_iPSC_3x", "2176866_iPSC_3x", "2186960_iPSC_3x", "2139701_iPSC_3x" ) 
dPSI_bi_batch2_iPSC = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25

# colnames(dPSI_bi_tmp) =  c("2184088_NGN2_3x", "2184090_NGN2_3x", "2174714_NGN2_3x", "2184741_NGN2_3x", "2174748_NGN2_3x", "2175420_NGN2_3x", 
#                            "2006152_NGN2_3x", "2186527_NGN2_3x", "2176866_NGN2_3x", "2186960_NGN2_3x", "2139701_NGN2_3x" ) 
# dPSI_bi_batch2_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25

head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)


# ##############################   1-Batch 3 ############################## ############################## ############################## 
din3 = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_1/analysis.03.splice_analysis/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin3 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#param = gsub(".csv","",gsub("binary_","",fin3))
dPSI_bi = read.csv( paste0(din3, fin3) , stringsAsFactors = F , row.names = 1 , check.names = F)
dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp =  dPSI_bi[     , grep('iPSC_',  colnames(dPSI_bi) ) ] 
#dPSI_bi_tmp =  dPSI_bi[     , grep('NGN2_',  colnames(dPSI_bi) ) ] 

dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('_3x'  ,  colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('bridge',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
head(dPSI_bi_tmp)

colnames(dPSI_bi_tmp)

colnames(dPSI_bi_tmp) =  c(  "1755497_iPSC_3x", "2006152_iPSC_3x", "2189972_iPSC_3x", "2194943_iPSC_3x", "2195127_iPSC_3x", "2195327_iPSC_3x", "2197294_iPSC_3x" ) 
dPSI_bi_batch3_iPSC = dPSI_bi_tmp     #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425

#colnames(dPSI_bi_tmp) =  c(  "1755497_NGN2_3x", "2006152_NGN2_3x", "2189972_NGN2_3x", "2194943_NGN2_3x", "2195127_NGN2_3x", "2195327_NGN2_3x", "2197294_NGN2_3x" ) 
#dPSI_bi_batch3_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425

head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)


####### iPSC
##############################=========      2023-04-18 ===========================
# 1) Ven diagrams of BIO-2184090,  BIO-2186960 , BIO-2195127  à what is the percent overlap of off-targets? 2184090    2186960 Batch2         BIO-2195127 batch3
# 2) Compound 6152 was run in all 3 batches, create list of off-targets for iPSC for each run. How much overlap between runs?

head(dPSI_bi_batch2_iPSC)
head(dPSI_bi_batch3_iPSC)

iPSC_2184090_3x_vector =   rownames( dPSI_bi_batch2_iPSC[  (dPSI_bi_batch2_iPSC$`2184090_iPSC_3x` == 1), ] ) # 6960 from batch2
iPSC_2184090_3x_vector

iPSC_2186960_3x_vector =   rownames( dPSI_bi_batch2_iPSC[  (dPSI_bi_batch2_iPSC$`2186960_iPSC_3x` == 1), ] ) # 2186960
iPSC_2186960_3x_vector

iPSC_2195127_3x_vector =   rownames( dPSI_bi_batch3_iPSC[  (dPSI_bi_batch3_iPSC$`2195127_iPSC_3x` == 1), ] ) # PTC518
iPSC_2195127_3x_vector

library(stringr)
library(VennDiagram)
library(dplyr)


getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(iPSC_2184090_3x_vector, iPSC_2186960_3x_vector, iPSC_2195127_3x_vector  ),
  #category.names = c("2184090", "2186960", "2195127"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2184090_2186960_2195127.3x_iPSC.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2184090_2186960_2195127.3x_iPSC_no_name.png'),      ###### output in the working directory!!!
  output=TRUE
)
dev.off()


#######NGN2
##############################=========      2023-04-18 ===========================
# 1) Ven diagrams of BIO-2184090,  BIO-2184090 , BIO-2195127  à what is the percent overlap of off-targets? 4090 6960 Batch2 BIO-2195127 batch3
# 2) Compound 6152 was run in all 3 batches, create list of off-targets for iPSC for each run. How much overlap between runs?


head(dPSI_bi_batch2_NGN2) #2
head(dPSI_bi_batch3_NGN2) #3

NGN2_2184090_3x_vector =   rownames( dPSI_bi_batch2_NGN2[  (dPSI_bi_batch2_NGN2$`2184090_NGN2_3x` == 1), ] ) # 6960 from batch2
NGN2_2184090_3x_vector

NGN2_2186960_3x_vector =   rownames( dPSI_bi_batch2_NGN2[  (dPSI_bi_batch2_NGN2$`2186960_NGN2_3x` == 1), ] ) # 2186960
NGN2_2186960_3x_vector

NGN2_2195127_3x_vector =   rownames( dPSI_bi_batch3_NGN2[  (dPSI_bi_batch3_NGN2$`2195127_NGN2_3x` == 1), ] ) # PTC518
NGN2_2195127_3x_vector

getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(NGN2_2184090_3x_vector, NGN2_2186960_3x_vector, NGN2_2195127_3x_vector  ),
  #category.names = c("2184090", "2186960", "2195127"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2184090_2186960_2195127.3x_NGN2.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2184090_2186960_2195127.3x_no_name_NGN2.png'),      ###### output in the working directory!!!
  output=TRUE
)

dev.off()


########

# =============================================== PTC518
comp_3_title     = "PTC-518-3x: iPSC vs NGN2 vs Sy5Y"
comp_3_compounds = c("iPSC_2197294_3x",  "NGN2_2197294_3x" ,  "Sy5Y_2197294_3x" ) 

# comp_3_compounds_mx =   prepair_venn_plot_data( compounds = comp_3_compounds )   
# comp_3_compounds_mx =  comp_3_compounds_mx[ which( rowSums(is.na(comp_3_compounds_mx)) < ncol( comp_3_compounds_mx )  ) , ] #comp_3_compounds_mx =   comp_3_compounds_mx[ complete.cases( comp_3_compounds_mx )   , ]
# head(comp_3_compounds_mx)
# dim(comp_3_compounds_mx)
# write.table(comp_3_compounds_mx, "comp3_1st_2nd_NGN2_6152_3x_DSG_genes.csv", row.names=T, sep = ",", col.names=NA) # write.table(comp_3_compounds, "comp_3_genes.tsv", row.names=T, sep = "\t")

comp_3_compounds_df = as.data.frame(comp_3_compounds_mx)
head(comp_3_compounds_df)

iPSC_2197294_3x_vector =   rownames( dPSI_bi_tmp_PTC518[  (dPSI_bi_tmp_PTC518$iPSC_2197294_3x == 1), ] )
NGN2_2197294_3x_vector =   rownames( dPSI_bi_tmp_PTC518[  (dPSI_bi_tmp_PTC518$NGN2_2197294_3x == 1), ] )
Sy5Y_2197294_3x_vector =   rownames( dPSI_bi_tmp_PTC518[  (dPSI_bi_tmp_PTC518$Sy5Y_2197294_3x == 1), ] )


library(stringr)
library(VennDiagram)
library(dplyr)


getwd()
out_dir=getwd()
dev.off()
venn.diagram(
  x = list(iPSC_2197294_3x_vector, NGN2_2197294_3x_vector, Sy5Y_2197294_3x_vector ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  category.names = c("", "",  ""),
  imagetype = "png",
  filename = paste0(out_dir, '/venn_diagramm_compound_PTC518.3x_in_3_cells_no_title.png'),      ###### output in the working directory!!!
  output=TRUE
)
dev.off()

################################################## ======== ========== ##################################################
#### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels
y_1st
y_2nd
y_3rd
# y_1st_BK = y_1st
# y_2nd_BK = y_2nd 
# y_3rd_BK = y_3rd
y_1st_df  = as.data.frame(y_1st) 
y_2nd_df  = as.data.frame(y_2nd)
y_3rd_df  = as.data.frame(y_3rd)

y_12_df  = merge(y_1st_df, y_2nd_df, by = 0, all = T)
y_123_df = merge(y_12_df, y_3rd_df, by.x = "Row.names",  by.y = 0, all = T)

y_123_df

y_123_df$Row.names <- str_replace(y_123_df$Row.names, "_iPSC_3x", "") #library(stringr)

y_123_df_Batch1  = y_123_df[ ,c ("Row.names",  "y_1st")] #y_123_df_Batch1$y_1st[is.na(y_123_df_Batch1$y_1st )]  = 0
y_123_df_Batch2  = y_123_df[ ,c ("Row.names",  "y_2nd")] #y_123_df_Batch2$y_2nd[is.na(y_123_df_Batch2$y_2nd )]  = 0
y_123_df_Batch3  = y_123_df[ ,c ("Row.names",  "y_3rd")] #y_123_df_Batch3$y_3rd[is.na(y_123_df_Batch3$y_3rd)]  = 0 # NA will not show 0 in plots

y_123_df_Batch1 # axis_labels <- y_123_df_Batch1$Row.names

############### 1-1 raw-3 plots format
dev.off()
p1 <- ggplot(y_123_df_Batch1, aes(x=Row.names, y=y_1st ,label = y_1st ) )+ 
  ylim(0,max(y_1st)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels = y_123_df_Batch1$Row.names, guide = guide_axis(angle = 90) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  #+annotate(  "text", label = "iPSC cells",     x = 12, y = 15, size = 6, colour = "red"  )
p2 <- ggplot(y_123_df_Batch2, aes(x=Row.names, y=y_2nd, label = y_2nd ) ) + 
  ylim(0,max(y_2nd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels =  y_123_df_Batch2$Row.names) +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  
p3 <- ggplot(y_123_df_Batch3, aes(x=Row.names, y=y_3rd , label = y_3rd) ) + #p3 <- ggplot(y_123_df_Batch3, aes(x=factor(Row.names), y=y_3rd , label = y_3rd) ) + 
  ylim(0,max(y_3rd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  + ggtitle("DSG induced in iPSC cells with 3x concentration") #+   xlab("Dose (mg)") + ylab("Teeth length")

p3 / p2 / p1 +  plot_layout(guides = 'collect') # save as pdf iPSC_3_plots_alphabetical 7x7 DSG_3_panel_rankied_iPSC 7x9


############### 1.2 ranked  3 plots format
# y_1st_df  = as.data.frame(y_1st) 
# y_2nd_df  = as.data.frame(y_2nd)
# y_3rd_df  = as.data.frame(y_3rd)
y_123_df # Un-sorted_for_raw_alphabeticall range

y_1st_df
y_2nd_df  
y_3rd_df  

y_1st_df$compound  = row.names(y_1st_df) 
y_2nd_df$compound  = row.names(y_2nd_df)
y_3rd_df$compound  = row.names(y_3rd_df)

class(y_1st_df)
y_1st_df_sorted = y_1st_df[ order(y_1st_df$y_1st, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector y_1st_df_sorted = y_1st_df[ order(y_1st_df$y_1st, decreasing=TRUE ), c(1), drop = FALSE]  # add c(1) avoid degraded to vector
y_2nd_df_sorted = y_2nd_df[ order(y_2nd_df$y_2nd, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector
y_3rd_df_sorted = y_3rd_df[ order(y_3rd_df$y_3rd, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector
y_1st_df_sorted
y_2nd_df_sorted
y_3rd_df_sorted


y_12_df_sorted  = full_join(y_1st_df_sorted, y_2nd_df_sorted, all = T )  # keep the order
y_12_df_sorted
y_123_df_sorted = full_join(y_12_df_sorted,  y_3rd_df_sorted, all = T )  
y_123_df_sorted

# y_1st        compound y_2nd y_3rd
# 1    189 2060573_iPSC_3x    NA    NA
# 2    185 2070692_iPSC_3x    NA    NA
# 3    147 2059811_iPSC_3x    NA    NA
# 4    116 2135644_iPSC_3x    NA    NA
# 5    103 2136770_iPSC_3x    NA    NA
# 8     NA 2176866_iPSC_3x    59    NA
# 9     NA 2175420_iPSC_3x    57    NA
# 10    NA 2184741_iPSC_3x    55    NA
# 11    NA 2184088_iPSC_3x    47    NA
# 18    NA 2189972_iPSC_3x    NA    39
# 12    NA 2174714_iPSC_3x    36    NA
# 19    NA 2195327_iPSC_3x    NA    36
# 6     44 2006152_iPSC_3x    32    31
# 13    NA 2139701_iPSC_3x    32    NA
# 14    NA 2174748_iPSC_3x    32    NA
# 15    NA 2184090_iPSC_3x    31    NA
# 20    NA 2195127_iPSC_3x    NA    31
# 21    NA 2197294_iPSC_3x    NA    25
# 16    NA 2186527_iPSC_3x    24    NA
# 17    NA 2186960_iPSC_3x    24    NA
# 7     20 1755497_iPSC_3x    NA    23
# 22    NA 2194943_iPSC_3x    NA    20
y_123_df_sorted  = y_123_df_sorted[c(1:5, 8:11, 18,12,19, 6, 13:15, 20:21, 16:17, 7, 22), ] 

y_123_df_sorted
y_123_df_sorted$compound
y_123_df_sorted$compound <- str_replace(y_123_df_sorted$compound, "_iPSC_3x", "") #library(stringr) #print(df)
y_123_df_sorted$compound

y_123_df_sorted$compound[ which(y_123_df_sorted$compound == "1755497")] <- "Branaplam" # https://www.geeksforgeeks.org/how-to-replace-specific-values-in-column-in-r-dataframe/
y_123_df_sorted$compound[ which(y_123_df_sorted$compound == "2197294")] <- "PTC-518" # ======PTC518 compound is BIO-2197294
y_123_df_sorted # 

# library(dplyr)
# library(forcats)  # fct_relevel
y_123_df_sorted$compound
y_123_df_sorted <- y_123_df_sorted %>% # qqplot automatically sort with string alphabetically, has to change to factor, and change their level
  mutate(
    compound = fct_relevel(compound, 
                           "2060573", "2070692", "2059811", "2135644", "2136770", "2176866", "2175420", "2184741",   "2184088", "2189972", "2174714", "2195327", "2006152",  
                           "2139701", "2174748", "2184090", "2195127", "PTC-518" , "2186527", "2186960", "Branaplam", "2194943"             
                           ))
levels(y_123_df_sorted$compound)

y_123_df_Batch1  = y_123_df_sorted[ ,c ("compound",  "y_1st" )]  #, "Row.name")]  #
y_123_df_Batch2  = y_123_df_sorted[ ,c ("compound",  "y_2nd" )]  #, "Row.name")] 
y_123_df_Batch3  = y_123_df_sorted[ ,c ("compound",  "y_3rd" )]  #, "Row.name")] 

y_123_df_Batch1
y_123_df_Batch1$compound # y_123_df_Batch1$Row.name

y_123_df_Batch1
y_123_df_Batch2
y_123_df_Batch3

dev.off() #
p1 <- ggplot(y_123_df_Batch1, aes(x=compound, y=y_1st ,label = y_1st ) )+ 
  ylim(0,max(y_1st)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels = y_123_df_Batch1$compound, guide = guide_axis(angle = 90) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  #+   annotate(  "text", label = "iPSC cells with 3x concentration",     x = 10, y = 2, size = 6, colour = "red"  )


p2 <- ggplot(y_123_df_Batch2, aes(x=compound, y=y_2nd, label = y_2nd ) ) + 
  ylim(0,max(y_2nd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels =  y_123_df_Batch2$compound) +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  

#p3 <- ggplot(y_123_df_Batch3, aes(x=factor(compound), y=y_3rd , label = y_3rd) ) + 
p3 <- ggplot(y_123_df_Batch3, aes(x=compound, y=y_3rd , label = y_3rd) ) + 
  ylim(0,max(y_3rd)*1.5) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  + ggtitle("DSG induced in iPSC cells with 3x concentration") #+   xlab("Dose (mg)") + ylab("Teeth length")

p3 + p2 + p1 + plot_layout(heights = c(1, 1.6, 5.2))  # p1 + p2 + p3 + p4 +  plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
getwd() # DSG_in_iPSC_3_panel.pdf 7 x 9      pdf iPSC_3_panel_ranked_DSG

#### ==========+++++++++++========== 1-3 united plot 
y_1st_df  = as.data.frame(y_1st) 
y_2nd_df  = as.data.frame(y_2nd)
y_3rd_df  = as.data.frame(y_3rd)

y_1st_df$batch  = "Batch1" 
y_2nd_df$batch  = "Batch2" 
y_3rd_df$batch  = "Batch3" 
y_3rd_df

rownames(y_1st_df)[rownames(y_1st_df) == "1755497_iPSC_3x"] <- "Branaplam_Batch_1"
rownames(y_3rd_df)[rownames(y_3rd_df) == "1755497_iPSC_3x"] <- "Branaplam_Batch_3"

rownames(y_1st_df)[rownames(y_1st_df) == "2006152_iPSC_3x"] <- "2006152_Batch_1"
rownames(y_2nd_df)[rownames(y_2nd_df) == "2006152_iPSC_3x"] <- "2006152_Batch_2"
rownames(y_3rd_df)[rownames(y_3rd_df) == "2006152_iPSC_3x"] <- "2006152_Batch_3"
#  ======PTC518 compound is BIO-2197294
rownames(y_3rd_df)[rownames(y_3rd_df) == "2197294_iPSC_3x"] <- "PTC-518"

names(y_1st_df)[names(y_1st_df) == "y_1st"] <- "DSG_count"
names(y_2nd_df)[names(y_2nd_df) == "y_2nd"] <- "DSG_count"
names(y_3rd_df)[names(y_3rd_df) == "y_3rd"] <- "DSG_count"

y_1st_df
y_2nd_df
y_3rd_df

y_123_df = rbind(y_1st_df, y_2nd_df, y_3rd_df) #y_123_df_back = y_123_df
y_123_df

Average_of_Branaplam =  c( round( mean( y_123_df[ c('Branaplam_Batch_1', 'Branaplam_Batch_3'                ), c("DSG_count") ] ) , 1 ), "Average" )
Average_of_2006152   =  c( round( mean( y_123_df[ c('2006152_Batch_1'  ,'2006152_Batch_2', '2006152_Batch_3'), c("DSG_count") ] ) , 1),  "Average" )


y_123_df= rbind (y_123_df, Branaplam = Average_of_Branaplam,  "2006152" = Average_of_2006152) # this way add row name
y_123_df = y_123_df[ !(row.names(y_123_df) %in% c( "Branaplam_Batch_1" , "Branaplam_Batch_3", "2006152_Batch_1", "2006152_Batch_2", "2006152_Batch_3" ) ), ]
y_123_df

row.names(y_123_df) <- str_replace( row.names(y_123_df), "_iPSC_3x", "") #library(stringr) OR sub

y_123_df$DSG_count = as.numeric(y_123_df$DSG_count)
y_123_df$compounds = rownames(y_123_df)
y_123_df = y_123_df[order(y_123_df$DSG_count, decreasing=TRUE), ] 
y_123_df$compounds
y_123_df 

str(y_123_df) 

########################
y_123_df$compounds
y_123_df <- y_123_df %>% # qqplot automatically sort with string alphabetically, has to change to factor, and change their level
  mutate(
    compounds = fct_relevel(compounds, 
                           "2060573","2070692","2059811","2135644","2136770","2176866","2175420","2184741","2184088","2189972","2174714","2195327","2006152","2139701","2174748","2184090","2195127",  
                           "PTC-518","2186527","2186960","Branaplam","2194943"            
    ))
levels(y_123_df$compounds)

dev.off()

#ggplot(y_123_df1, aes(x= reorder(compounds, -DSG_count) , y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
ggplot(y_123_df, aes(x= compounds, y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
  geom_bar(stat = "identity" ) +   #geom_bar(stat = "identity", fill="white") + 
  scale_x_discrete(labels = y_123_df$compounds, guide = guide_axis(angle = 90) ) +
  theme(plot.title  = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(20, 20), "pt"))  + ggtitle("DSG induced in iPSC cells with 3x concentration") +
  geom_text(size = 3, vjust = -0.5) 
# save as iPSC_merged_DSG 7X9


# ##############################   222222----------------------------------------- NGN2 Note ======PTC518 compound is BIO-2197294===============================================================================================================================================
# ##############################   2-Batch 1 ############################## ############################## ##############################  
din1 = "/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.03.splice_analysis/"
fin1 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
# param = gsub(".csv","",gsub("binary_","",fin1))
dPSI_bi = read.csv( paste0(din1, fin1)  , stringsAsFactors = F , row.names = 1 , check.names = F) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp            =  dPSI_bi[     , grep('NGN2-',  colnames(dPSI_bi) ) ] 
head(dPSI_bi_tmp)
colnames(dPSI_bi_tmp)  = gsub("BIO-", "", colnames(dPSI_bi_tmp ) )
colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c("1755497_NGN2_3x" , "1755497_NGN2_10x" , "1949634_NGN2_3x" , "1949634_NGN2_10x" , "2006152_NGN2_3x" , "2006152_NGN2_10x" , 
                           "2059811_NGN2_3x" , "2059811_NGN2_10x" , "2060573_NGN2_3x" , "2060573_NGN2_10x" , "2060884_NGN2_3x" , "2060884_NGN2_10x" , 
                           "2070692_NGN2_3x" , "2070692_NGN2_10x" , "2135644_NGN2_3x" , "2135644_NGN2_10x" , "2136770_NGN2_3x" , "2136770_NGN2_10x" ) 

# 1755497 is Branaplam, 1949634 Risdiplam, 206088s Interal Lead
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('_3x'  ,  colnames( dPSI_bi_tmp) )] #
head(dPSI_bi_tmp)
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('1949634',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('2060884',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
head(dPSI_bi_tmp)

dPSI_bi = dPSI_bi_tmp
y_1st        = colSums(dPSI_bi) ##################### get the totall DS events
y_1st        = y_1st[order(names(y_1st) )]
y_1st # class(y_1st) #y_1st_rank = rank(y_1st ) # y_rank = rank(-y )

# ##############################   2-Batch 2 ############################## ############################## ##############################  
din2 = "/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold_3v3_to_show_10-13/backup-good/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin2 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
param = gsub(".csv","",gsub("binary_","",fin1))
dPSI_bi = read.csv( paste0(din2, fin2)  , stringsAsFactors = F , row.names = 1 , check.names = F) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp            =  dPSI_bi[     , grep('NGN2-',  colnames(dPSI_bi) ) ] 
dPSI_bi_tmp            =  dPSI_bi_tmp[ , grep('-3x'  ,  colnames( dPSI_bi_tmp) )] #
head(dPSI_bi_tmp) # colnames(dPSI_bi_tmp)  = gsub("BIO-", "", colnames(dPSI_bi_tmp ) )
colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c("2184088_NGN2_3x", "2184090_NGN2_3x", "2174714_NGN2_3x", "2184741_NGN2_3x", "2174748_NGN2_3x", "2175420_NGN2_3x", 
                           "2006152_NGN2_3x", "2186527_NGN2_3x", "2176866_NGN2_3x", "2186960_NGN2_3x", "2139701_NGN2_3x" ) 
head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)

dPSI_bi = dPSI_bi_tmp
y_2nd        = colSums(dPSI_bi) ##################### get the totall DS events
y_2nd        = y_2nd[order(names(y_2nd) )]
y_2nd # class(y_2nd) #y_2nd_rank = rank(y_2nd ) # y_rank = rank(-y ) #y_2nd_rank


# ##############################   2-Batch 3 ############################## ############################## ############################## 
din3 = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_1/analysis.03.splice_analysis/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin3 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#param = gsub(".csv","",gsub("binary_","",fin3))
dPSI_bi = read.csv( paste0(din3, fin3) , stringsAsFactors = F , row.names = 1 , check.names = F)
dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp =  dPSI_bi[     , grep('NGN2_',  colnames(dPSI_bi) ) ] 
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('_3x'  ,  colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('bridge',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
head(dPSI_bi_tmp)

colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c(  "1755497_NGN2_3x", "2006152_NGN2_3x", "2189972_NGN2_3x", "2194943_NGN2_3x", "2195127_NGN2_3x", "2195327_NGN2_3x", "2197294_NGN2_3x" ) 

dPSI_bi = dPSI_bi_tmp
y_3rd        = colSums(dPSI_bi) ##################### get the totall DS events
y_3rd        = y_3rd[order(names(y_3rd) , decreasing = T)]
y_3rd #class(y_3rd)  #y_3rd_rank = rank(y_3rd ) # y_rank = rank(-y )

#### ==========+++++++++++========== 2-1  ==> 1 plot with 3 panels
y_1st
y_2nd
y_3rd

y_1st_df  = as.data.frame(y_1st) 
y_2nd_df  = as.data.frame(y_2nd)
y_3rd_df  = as.data.frame(y_3rd)

y_12_df  = merge(y_1st_df, y_2nd_df, by = 0, all = T)
y_123_df = merge(y_12_df, y_3rd_df, by.x = "Row.names",  by.y = 0, all = T)

y_123_df

y_123_df$Row.names <- str_replace(y_123_df$Row.names, "_NGN2_3x", "") #library(stringr)

y_123_df_Batch1  = y_123_df[ ,c ("Row.names",  "y_1st")] #y_123_df_Batch1$y_1st[is.na(y_123_df_Batch1$y_1st )]  = 0
y_123_df_Batch2  = y_123_df[ ,c ("Row.names",  "y_2nd")] #y_123_df_Batch2$y_2nd[is.na(y_123_df_Batch2$y_2nd )]  = 0
y_123_df_Batch3  = y_123_df[ ,c ("Row.names",  "y_3rd")] #y_123_df_Batch3$y_3rd[is.na(y_123_df_Batch3$y_3rd)]  = 0 # NA will not show 0 in plots

y_123_df_Batch1 # axis_labels <- y_123_df_Batch1$Row.names

############### 2-1 raw-3 plots format
dev.off()
p1 <- ggplot(y_123_df_Batch1, aes(x=Row.names, y=y_1st ,label = y_1st ) )+ 
  ylim(0,max(y_1st)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels = y_123_df_Batch1$Row.names, guide = guide_axis(angle = 90) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  #+annotate(  "text", label = "NGN2 cells",     x = 12, y = 15, size = 6, colour = "red"  )
p2 <- ggplot(y_123_df_Batch2, aes(x=Row.names, y=y_2nd, label = y_2nd ) ) + 
  ylim(0,max(y_2nd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels =  y_123_df_Batch2$Row.names) +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  
p3 <- ggplot(y_123_df_Batch3, aes(x=Row.names, y=y_3rd , label = y_3rd) ) + #p3 <- ggplot(y_123_df_Batch3, aes(x=factor(Row.names), y=y_3rd , label = y_3rd) ) + 
  ylim(0,max(y_3rd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  + ggtitle("DSG induced in NGN2 cells with 3x concentration") #+   xlab("Dose (mg)") + ylab("Teeth length")

p3 / p2 / p1 +  plot_layout(guides = 'collect') # save as pdf NGN2_3_plots_alphabetical_compound 7x7 DSG_3_panel_rankied_NGN2 7x9


############### 2-2 ranked  3 plots format
y_1st_df
y_2nd_df  
y_3rd_df  

y_1st_df$compound  = row.names(y_1st_df) 
y_2nd_df$compound  = row.names(y_2nd_df)
y_3rd_df$compound  = row.names(y_3rd_df)

class(y_1st_df)
y_1st_df_sorted = y_1st_df[ order(y_1st_df$y_1st, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector y_1st_df_sorted = y_1st_df[ order(y_1st_df$y_1st, decreasing=TRUE ), c(1), drop = FALSE]  # add c(1) avoid degraded to vector
y_2nd_df_sorted = y_2nd_df[ order(y_2nd_df$y_2nd, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector
y_3rd_df_sorted = y_3rd_df[ order(y_3rd_df$y_3rd, decreasing=TRUE ), ]  # add c(1) avoid degraded to vector
y_1st_df_sorted
y_2nd_df_sorted
y_3rd_df_sorted


y_12_df_sorted  = full_join(y_1st_df_sorted, y_2nd_df_sorted, all = T )  # keep the order
y_12_df_sorted
y_123_df_sorted = full_join(y_12_df_sorted,  y_3rd_df_sorted, all = T )  
y_123_df_sorted
# y_1st        compound y_2nd y_3rd
# 1    433 2070692_NGN2_3x    NA    NA
# 2    366 2060573_NGN2_3x    NA    NA
# 3    316 2059811_NGN2_3x    NA    NA
# 4    245 2136770_NGN2_3x    NA    NA
# 5    228 2135644_NGN2_3x    NA    NA

# 8     NA 2175420_NGN2_3x   131    NA
# 9     NA 2184741_NGN2_3x   105    NA
# 10    NA 2174714_NGN2_3x    80    NA
# 11    NA 2184088_NGN2_3x    76    NA

# 6     95 2006152_NGN2_3x    66    63
# 7     73 1755497_NGN2_3x    NA    79

# 12    NA 2176866_NGN2_3x    62    NA
# 13    NA 2186960_NGN2_3x    61    NA
# 14    NA 2174748_NGN2_3x    60    NA
# 15    NA 2184090_NGN2_3x    55    NA
# 16    NA 2186527_NGN2_3x    52    NA

# 18    NA 2194943_NGN2_3x    NA    47

# 17    NA 2139701_NGN2_3x    40    NA

# 19    NA 2195127_NGN2_3x    NA    43
# 20    NA 2195327_NGN2_3x    NA    40
# 21    NA 2189972_NGN2_3x    NA    40
# 22    NA 2197294_NGN2_3x    NA    38


y_123_df_sorted  = y_123_df_sorted[c(1:5, 8:11, 6:7,12:16, 18,17, 19:22), ]  #######<<<<<<<<<<<<<<<<<<<==================== MANUALLY check!!!!!!!!!!!!!

y_123_df_sorted
y_123_df_sorted$compound
y_123_df_sorted$compound <- str_replace(y_123_df_sorted$compound, "_NGN2_3x", "") #library(stringr) #print(df)
y_123_df_sorted$compound

y_123_df_sorted$compound[ which(y_123_df_sorted$compound == "1755497")] <- "Branaplam" # https://www.geeksforgeeks.org/how-to-replace-specific-values-in-column-in-r-dataframe/
y_123_df_sorted$compound[ which(y_123_df_sorted$compound == "2197294")] <- "PTC-518" # ======PTC518 compound is BIO-2197294
y_123_df_sorted # 

# library(dplyr)
# library(forcats)  # fct_relevel
y_123_df_sorted$compound #  "2070692"   "2060573"   "2059811"   "2136770"   "2135644"   "2175420"   "2184741"   "2174714"   "2184088"   "2006152"   "Branaplam" "2176866"   "2186960"   "2174748"   "2184090"   "2186527"   "2194943"   "2139701"   "2195127"   "2195327"   "2189972"   "PTC-518"  
y_123_df_sorted <- y_123_df_sorted %>% # qqplot automatically sort with string alphabetically, has to change to factor, and change their level
  mutate(
    compound = fct_relevel(compound, 
                           "2070692", "2060573", "2059811", "2136770", "2135644", "2175420", "2184741", "2174714", "2184088", "2006152", "Branaplam",
                           "2176866", "2186960", "2174748", "2184090", "2186527", "2194943", "2139701", "2195127", "2195327", "2189972", "PTC-518"              
    ))
levels(y_123_df_sorted$compound)

y_123_df_Batch1  = y_123_df_sorted[ ,c ("compound",  "y_1st" )]  #, "Row.name")]  #
y_123_df_Batch2  = y_123_df_sorted[ ,c ("compound",  "y_2nd" )]  #, "Row.name")] 
y_123_df_Batch3  = y_123_df_sorted[ ,c ("compound",  "y_3rd" )]  #, "Row.name")] 

y_123_df_Batch1
y_123_df_Batch1$compound # y_123_df_Batch1$Row.name

y_123_df_Batch1
y_123_df_Batch2
y_123_df_Batch3

dev.off() #
p1 <- ggplot(y_123_df_Batch1, aes(x=compound, y=y_1st ,label = y_1st ) )+ 
  ylim(0,max(y_1st)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels = y_123_df_Batch1$compound, guide = guide_axis(angle = 90) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  #+   annotate(  "text", label = "NGN2 cells with 3x concentration",     x = 10, y = 2, size = 6, colour = "red"  )


p2 <- ggplot(y_123_df_Batch2, aes(x=compound, y=y_2nd, label = y_2nd ) ) + 
  ylim(0,max(y_2nd)*1.2) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_discrete(labels =  y_123_df_Batch2$compound) +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  

#p3 <- ggplot(y_123_df_Batch3, aes(x=factor(compound), y=y_3rd , label = y_3rd) ) + 
p3 <- ggplot(y_123_df_Batch3, aes(x=compound, y=y_3rd , label = y_3rd) ) + 
  ylim(0,max(y_3rd)*1.5) +
  geom_col(position = position_dodge()) +
  geom_text(size = 3, vjust = -0.5) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  + ggtitle("DSG induced in NGN2 cells with 3x concentration") #+   xlab("Dose (mg)") + ylab("Teeth length")

p3 + p2 + p1 + plot_layout(heights = c(1, 1.6, 5.2))  # p1 + p2 + p3 + p4 +  plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
getwd() #  save as NGN2_3_panel_ranked_DSG 7x9 pad

#### ==========+++++++++++========== 2-3 united plot 
y_1st_df  = as.data.frame(y_1st) 
y_2nd_df  = as.data.frame(y_2nd)
y_3rd_df  = as.data.frame(y_3rd)

y_1st_df$batch  = "Batch1" 
y_2nd_df$batch  = "Batch2" 
y_3rd_df$batch  = "Batch3" 

y_1st_df
y_2nd_df
y_3rd_df

rownames(y_1st_df)[rownames(y_1st_df) == "1755497_NGN2_3x"] <- "Branaplam_Batch_1"
rownames(y_3rd_df)[rownames(y_3rd_df) == "1755497_NGN2_3x"] <- "Branaplam_Batch_3"

rownames(y_1st_df)[rownames(y_1st_df) == "2006152_NGN2_3x"] <- "2006152_Batch_1"
rownames(y_2nd_df)[rownames(y_2nd_df) == "2006152_NGN2_3x"] <- "2006152_Batch_2"
rownames(y_3rd_df)[rownames(y_3rd_df) == "2006152_NGN2_3x"] <- "2006152_Batch_3"
#  ======PTC518 compound is BIO-2197294
rownames(y_3rd_df)[rownames(y_3rd_df) == "2197294_NGN2_3x"] <- "PTC-518"

names(y_1st_df)[names(y_1st_df) == "y_1st"] <- "DSG_count"
names(y_2nd_df)[names(y_2nd_df) == "y_2nd"] <- "DSG_count"
names(y_3rd_df)[names(y_3rd_df) == "y_3rd"] <- "DSG_count"

y_1st_df
y_2nd_df
y_3rd_df

y_123_df = rbind(y_1st_df, y_2nd_df, y_3rd_df) #y_123_df_back = y_123_df
y_123_df

Average_of_Branaplam =  c( round( mean( y_123_df[ c('Branaplam_Batch_1', 'Branaplam_Batch_3'                ), c("DSG_count") ] ) , 1 ), "Average" )
Average_of_2006152   =  c( round( mean( y_123_df[ c('2006152_Batch_1'  ,'2006152_Batch_2', '2006152_Batch_3'), c("DSG_count") ] ) , 1),  "Average" )


y_123_df= rbind (y_123_df, Branaplam = Average_of_Branaplam,  "2006152" = Average_of_2006152) # this way add row name
y_123_df = y_123_df[ !(row.names(y_123_df) %in% c( "Branaplam_Batch_1" , "Branaplam_Batch_3", "2006152_Batch_1", "2006152_Batch_2", "2006152_Batch_3" ) ), ]
y_123_df

row.names(y_123_df) <- str_replace( row.names(y_123_df), "_NGN2_3x", "") #library(stringr) OR sub

y_123_df$DSG_count = as.numeric(y_123_df$DSG_count)
y_123_df$compounds = rownames(y_123_df)
y_123_df = y_123_df[order(y_123_df$DSG_count, decreasing=TRUE), ] 
y_123_df$compounds
y_123_df 

str(y_123_df) 
y_123_df 

########################
y_123_df$compounds #[1] "2070692"   "2060573"   "2059811"   "2136770"   "2135644"   "2175420"   "2184741"   "2174714"   "2184088"   "Branaplam" "2006152"   "2176866"   "2186960"   "2174748"   "2184090"   "2186527"   "2194943"   "2195127"   "2139701"   "2195327"   "2189972"   "PTC-518"  
y_123_df <- y_123_df %>% # qqplot automatically sort with string alphabetically, has to change to factor, and change their level
  mutate(
    compounds = fct_relevel(compounds, 
                            "2070692", "2060573", "2059811", "2136770", "2135644", "2175420", "2184741", "2174714", "2184088", "Branaplam",
                            "2006152", "2176866", "2186960", "2174748", "2184090", "2186527", "2194943", "2195127", "2139701", "2195327", "2189972", "PTC-518"          
    ))
levels(y_123_df$compounds)

dev.off()

#ggplot(y_123_df1, aes(x= reorder(compounds, -DSG_count) , y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
ggplot(y_123_df, aes(x= compounds, y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
  geom_bar(stat = "identity" ) +   #geom_bar(stat = "identity", fill="white") + 
  scale_x_discrete(labels = y_123_df$compounds, guide = guide_axis(angle = 90) ) +
  theme(plot.title  = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(20, 20), "pt"))  + ggtitle("DSG induced in NGN2 cells with 3x concentration") +
  geom_text(size = 3, vjust = -0.5) 
# save as NGN2_merged_DSG_labled 7 X 9


###############################   333333----------------------------------------- Sy5Y   ====================================================================================================================================================================================
###############################   3-Batch 3 ############################## ############################## ############################## 
din3 = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_1/analysis.03.splice_analysis/" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
fin3 = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
#param = gsub(".csv","",gsub("binary_","",fin3))
dPSI_bi = read.csv( paste0(din3, fin3) , stringsAsFactors = F , row.names = 1 , check.names = F)
dim(dPSI_bi)
head(dPSI_bi[, 1:10])
head(dPSI_bi)
dim(dPSI_bi)

## total number of DSG ###
dPSI_bi_tmp =  dPSI_bi[     , grep('Sy5Y_',  colnames(dPSI_bi) ) ] 
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('_3x'  ,  colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
dPSI_bi_tmp =  dPSI_bi_tmp[ , grep('bridge',invert = TRUE, colnames( dPSI_bi_tmp) )] # grep('_(abc|zxy)',str, value = TRUE) # https://stackoverflow.com/questions/15284316/grep-on-two-strings
head(dPSI_bi_tmp)

colnames(dPSI_bi_tmp)
colnames(dPSI_bi_tmp) =  c(  "1755497_Sy5Y_3x", "2006152_Sy5Y_3x", "2189972_Sy5Y_3x", "2194943_Sy5Y_3x", "2195127_Sy5Y_3x", "2195327_Sy5Y_3x", "2197294_Sy5Y_3x" ) 

dPSI_bi = dPSI_bi_tmp
y_3rd        = colSums(dPSI_bi) ##################### get the totall DS events
y_3rd        = y_3rd[order(names(y_3rd) , decreasing = T)]
#names(y_3rd) = gsub("dPSI", "", names(y_3rd))
y_3rd #class(y_3rd)  #y_3rd_rank = rank(y_3rd ) # y_rank = rank(-y )

#### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels
#y_123_df$Row.names <- str_replace(y_123_df$Row.names, "_Sy5Y_3x", "") #library(stringr)
#### ==========+++++++++++========== 3-3 united plot 
y_3rd_df  = as.data.frame(y_3rd)

########################


########## 

y_3rd_df$batch  = "Batch3" 

rownames(y_3rd_df)[rownames(y_3rd_df) == "1755497_Sy5Y_3x"] <- "Branaplam"
rownames(y_3rd_df)[rownames(y_3rd_df) == "2006152_Sy5Y_3x"] <- "2006152"
rownames(y_3rd_df)[rownames(y_3rd_df) == "2197294_Sy5Y_3x"] <- "PTC-518"
names(y_3rd_df)[names(y_3rd_df) == "y_3rd"] <- "DSG_count"
y_3rd_df

y_123_df = y_3rd_df


y_123_df

y_123_df$compounds = rownames(y_123_df)

row.names(y_123_df) <- str_replace( row.names(y_123_df), "_Sy5Y_3x", "") #library(stringr) OR sub
y_123_df$compounds <- str_replace( y_123_df$compounds, "_Sy5Y_3x", "") #library(stringr) OR sub

y_123_df$DSG_count = as.numeric(y_123_df$DSG_count)
y_123_df = y_123_df[order(y_123_df$DSG_count, decreasing=TRUE), ] 
y_123_df$compounds
y_123_df 
str(y_123_df) 

y_123_df$compounds
y_123_df <- y_123_df %>% # qqplot automatically sort with string alphabetically, has to change to factor, and change their level
  mutate(
    compounds = fct_relevel(compounds, 
                            "2006152","2195127","2195327","2189972","Branaplam" ,"PTC-518","2194943"  
    ))
levels(y_123_df$compounds)

dev.off()

#ggplot(y_123_df1, aes(x= reorder(compounds, -DSG_count) , y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
ggplot(y_123_df, aes(x= compounds, y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
  geom_bar(stat = "identity" ) +   #geom_bar(stat = "identity", fill="white") + 
  scale_x_discrete(labels = y_123_df$compounds, guide = guide_axis(angle = 90) ) +
  theme(plot.title  = element_text(color="red", size=14, face="bold.italic", hjust = 0.5) ) +
  theme(plot.margin = unit(rep(20, 20), "pt"))  + ggtitle("DSG induced in Sy5Y cells with 3x concentration") +
  geom_text(size = 3, vjust = -0.5) 
# save as Sy5Y_merged_DSG_labeled 6 X 7

ggplot(y_123_df1, aes(x= reorder(compounds, -DSG_count) , y = DSG_count, label = DSG_count, color=batch, fill=batch) ) +
  geom_bar(stat = "identity" ) +   #geom_bar(stat = "identity", fill="white") + 
  scale_x_discrete(labels = y_123_df1$compounds, guide = guide_axis(angle = 90) ) +
