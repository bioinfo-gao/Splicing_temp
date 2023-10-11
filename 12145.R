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


# RNA_Seq_raw2quickomics_12145_4-good2.2.R


require(tidyverse)
require(reshape2)
require(DESeq2)
require("UpSetR")
rm(list=ls())

Joon_out_dir     = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DEG_Results_Joon_1/"
DEG_gene_out_dir = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DEG_Results_Joon_1/Joon_HTT_MSH3_FoxM1/"
DEG_plot_out_dir = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DEG_Results_Joon_1/Bionary_Plot/"
dir.create(Joon_out_dir , recursive = T)
dir.create(DEG_gene_out_dir , recursive = T)
dir.create(DEG_plot_out_dir , recursive = T)

setwd(Joon_out_dir) # Joon getwd()

load("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/EA20230815_0/TST12145.RData") 
ls() # "comp_info"       "data_long"       "data_results"    "data_wide"       "MetaData"        "ProteinGeneName" "results_long"   
dim(ProteinGeneName) ###########
head(ProteinGeneName) ###########
tail(ProteinGeneName) ###########
results_long   
dim(results_long)   
head(results_long)   
tail(results_long)   

name_vector_Group1  = c( "DMSO-CHX__vs_DMSO" ,  
                         "U1-WT100-CHX_vs_DMSO" , "U1-WT100-CHX_vs_DMSO-CHX" , "U1-WT100-CHX_vs_WT100" , "U1-WT100-EtOH_vs_DMSO" , 
                         "U1-WT500-CHX_vs_DMSO" , "U1-WT500-CHX_vs_DMSO-CHX" , "U1-WT500-CHX_vs_WT500" , "U1-WT500-EtOH_vs_DMSO" )  

name_vector_Group2  = c( "U1-mut1-100-CHX_vs_DMSO"    , "U1-mut1-100-CHX_vs_DMSO-CHX" , "U1-mut1-100-CHX_vs_mut1-100" , "U1-mut1-100-EtOH_vs_DMSO" ,
                         "U1-mut1-100-CHX_vs_WT_EtOH" , "U1-mut1-100-CHX_vs_WT_CHX"   , "U1-mut1-100-EtOH_vs_WT_EtOH" ,
                         "U1-mut1-500-CHX_vs_DMSO"    , "U1-mut1-500-CHX_vs_DMSO-CHX" , "U1-mut1-500-CHX_vs_mut1-500" , "U1-mut1-500-EtOH_vs_DMSO" , 
                         "U1-mut1-500-CHX_vs_WT_EtOH" , "U1-mut1-500-CHX_vs_WT_CHX"   , "U1-mut1-500-EtOH_vs_WT_EtOH" )

name_vector_Group3  = c( "U1-mut2-100-CHX_vs_DMSO"    , "U1-mut2-100-CHX_vs_DMSO-CHX" , "U1-mut2-100-CHX_vs_mut2-100" , "U1-mut2-100-EtOH_vs_DMSO" ,
                         "U1-mut2-100-CHX_vs_WT_EtOH" , "U1-mut2-100-CHX_vs_WT_CHX"   , "U1-mut2-100-EtOH_vs_WT_EtOH" ,
                         "U1-mut2-500-CHX_vs_DMSO"    , "U1-mut2-500-CHX_vs_DMSO-CHX" , "U1-mut2-500-CHX_vs_mut5-500" , "U1-mut2-500-EtOH_vs_DMSO" ,
                         "U1-mut2-500-CHX_vs_WT_EtOH" , "U1-mut2-500-CHX_vs_WT_CHX"   , "U1-mut2-500-EtOH_vs_WT_EtOH" )
                              
name_vector_Group4  = c( "BIO2195127-CHX_vs_DMSO"                 , "BIO2195127-CHX_vs_DMSO-CHX"             , "BIO2195127-CHX_vs_5127"               , "BIO2195127-EtOH_vs_DMSO" , 
                         "BIO2195127-EtOH_vs_U1_WT100_EtOH"       , "BIO2195127-EtOH_vs_U1_WT500_EtOH"       , "BIO2195127-CHX_vs_U1_WT100_CHX"       , "BIO2195127-CHX_vs_U1_WT500_CHX" ,  
                         "BIO2195127-EtOH_vs_U1_U1_mut1_100_EtOH" , "BIO2195127-EtOH_vs_U1_U1_mut1_500_EtOH" , "BIO2195127-CHX_vs_U1_U1_mut1_100_CHX" , "BIO2195127-CHX_vs_U1_U1_mut1_500_CHX" , 
                         "BIO2195127-EtOH_vs_U1_U1_mut2_100_EtOH" , "BIO2195127-EtOH_vs_U1_U1_mut2_500_EtOH" , "BIO2195127-CHX_vs_U1_U1_mut2_100_CHX" , "BIO2195127-CHX_vs_U1_U1_mut2_500_CHX" )

name_vector_Group5  = c( "BIO2197294-CHX_vs_DMSO"                 , "BIO2197294-CHX_vs_DMSO-CHX"             , "BIO2197294-CHX_vs_7294"               , "BIO2197294-EtOH_vs_DMSO" ,
                         "BIO2197294-EtOH_vs_U1_WT100_EtOH"       , "BIO2197294-EtOH_vs_U1_WT500_EtOH"       , "BIO2197294-CHX_vs_U1_WT100_CHX"       ,  "BIO2197294-CHX_vs_U1_WT500_CHX" ,
                         "BIO2197294-EtOH_vs_U1_U1_mut1_100_EtOH" , "BIO2197294-EtOH_vs_U1_U1_mut1_500_EtOH" , "BIO2197294-CHX_vs_U1_U1_mut1_100_CHX" ,  "BIO2197294-CHX_vs_U1_U1_mut1_500_CHX" ,  
                         "BIO2197294-EtOH_vs_U1_U1_mut2_100_EtOH" , "BIO2197294-EtOH_vs_U1_U1_mut2_500_EtOH" , "BIO2197294-CHX_vs_U1_U1_mut2_100_CHX" ,  "BIO2197294-CHX_vs_U1_U1_mut2_500_CHX" )

name_vector_Group6  = c( "Branaplam-CHX_vs_DMSO"                  , "Branaplam-CHX_vs_DMSO-CHX"              , "Branaplam-CHX_vs_Branaplam"           ,  "Branaplam-EtOH_vs_DMSO" ,
                         "Branaplam-EtOH_vs_U1_WT100_EtOH"        , "Branaplam-EtOH_vs_U1_WT500_EtOH"        , "Branaplam-CHX_vs_U1_WT100_CHX"        ,  "Branaplam-CHX_vs_U1_WT500_CHX" , 
                         "Branaplam-EtOH_vs_U1_U1_mut1_100_EtOH"  , "Branaplam-EtOH_vs_U1_U1_mut1_500_EtOH"  , "Branaplam-CHX_vs_U1_U1_mut1_100_CHX"  ,  "Branaplam-CHX_vs_U1_U1_mut1_500_CHX" , 
                         "Branaplam-EtOH_vs_U1_U1_mut2_100_EtOH"  , "Branaplam-EtOH_vs_U1_U1_mut2_500_EtOH"  , "Branaplam-CHX_vs_U1_U1_mut2_100_CHX"  ,  "Branaplam-CHX_vs_U1_U1_mut2_500_CHX" )

name_vector_Group7  = c( "U1-WT-500-BIO2197294_vs_DMSO"           ,  "U1-WT-500-BIO2197294_vs_WT-500"        , "U1-WT-500-Branaplam_vs_DMSO"          , "U1-WT-500-Branaplam_vs_WT-500" )
  

name_vector_list           = list( name_vector_Group1 , name_vector_Group2 , name_vector_Group3 , name_vector_Group4 , name_vector_Group5 , name_vector_Group6 , name_vector_Group7 ) 
names(name_vector_list)    = c( "Grp1_DSMO_and_WT" , "Grp2_WT_mut1" , "Grp3_WT_mut2" , "Grp4_Bio5127_and_mut1_mut2" , "Grp5_PTC518_and_mut1_mut2" , "Grp6_Branaplam_and_mut1_mut2" , "Grp7_PTC518_and_Branaplam_and_WT" ) # NO / should be file name, will cause path error

str(name_vector_list    )

## for plotting in loop
rotate_x <- function(data, column_to_plot, labels_vec, rot_angle, cex=0.7) {
    par(mar = c(15,15,1,1))
    plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n", ylab = "num.DSG",)
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=cex)  # CEX font size 1 a little too big for more than 6 groups , 3 is HUGE
}

### 
for(j in 1:length(name_vector_list )){
    #  = 1  
    #j = 2
    name_vector        =       name_vector_list[[j]] 
    name_vector_string = names(name_vector_list)[j] 
    
    results_long_list  = list() 
    results_short_list = list() 
    DEG_with_gn = list()
    
    for(i in 1:length(name_vector )){
        #         i=1
        name_vector[[i]]
        
        results_long_list[[i]] = results_long[ results_long$test == name_vector[i], ]   
        #results_long_list[[i]] # results_long_list[[1]] 
        results_short_list[[i]] = results_long_list[[i]][, c("UniqueID" , "logFC" , "P.Value" , "Adj.P.Value")]
        head(results_short_list[[i]])
        dim(results_short_list[[i]])
        colnames(results_short_list[[i]]) = c( "UniqueID" , paste0( name_vector[[i]] , "_", "logFC") ,  "P.Value" , "Adj.P.Value"      ) 
        
        results_short_list[[i]] = results_short_list[[i]][ which(  results_short_list[[i]]$Adj.P.Value  < 1 & ! is.na(results_short_list[[i]]$Adj.P.Value  ) ),     ] 
        DEG_with_gn[[i]] = results_short_list[[i]][ which( abs( results_short_list[[i]][ , 2 ] )  > 0.585 & results_short_list[[i]]$Adj.P.Value < 0.05  ),     ] 
        
    }
    # str(results_short_list)
    # head(results_short_list[[1]])
    # head(  DEG_with_gn )
    # head(  DEG_with_gn[[1]] )
    # tail(  DEG_with_gn[[1]] )
    #  dim(  DEG_with_gn[[1]] )
    
    ### ===== out put the data Unless p.adj == 0 or NA
    raw_working_df = Reduce(function(x, y) merge(x, y, all=TRUE, by = "UniqueID"), results_short_list)
    head(raw_working_df )

    DEG_with_name_tmp = merge(  ProteinGeneName[ ,c( "UniqueID", "Gene.Name") ] ,raw_working_df , sort = F, all.y = T)
    DEG_with_name_tmp = DEG_with_name_tmp[  order( - abs(DEG_with_name_tmp[ , 3 ] )  ) ,  ] #     DEG_with_name_tmp[  order( abs(DEG_with_name_tmp[ , 3 ] ) ) ,  ]
    head(DEG_with_name_tmp)
    fout                = paste0(name_vector_string,"_DEG_genes.csv")
    #write.csv(DEG_with_name_tmp ,          file=fout,                quote=F, row.names=F)

    ### ===== out put the data HTT, MSH3 and FOXM1 _ DIferferent FOLDER
    DEG_with_selected_genes = DEG_with_name_tmp[ ( DEG_with_name_tmp$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,                    ]
    DEG_with_selected_genes
    HTT_MSH3_FoxM1_fout = paste0(DEG_gene_out_dir, name_vector_string,"_for_HTT_MSH3_FoxM1.csv")
    #write.csv(DEG_with_selected_genes, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=F)

    ##############################    ##############################    ##############################     ####  ======== UPSet Plotting. for upset plot _logFC.  
    count_working_df = Reduce(function(x, y) merge(x, y, all=TRUE, by = "UniqueID"), DEG_with_gn ) 
    head(count_working_df )
    DEG_with_name_count = merge(  ProteinGeneName[ ,c( "UniqueID", "Gene.Name") ] , count_working_df , sort = F, all.y = T) 
    head(DEG_with_name_count)
    
    working_DEG_count_df  = DEG_with_name_count[, grep( "^UniqueID|^Gene.Name|_logFC", colnames(DEG_with_name_count) ) ] 
    
    dim(working_DEG_count_df )
    head(working_DEG_count_df )
    working_DEG_count_df[1:5, 1:5]
    
    colnames(working_DEG_count_df ) 
    
    row.names(working_DEG_count_df )  = make.names(    working_DEG_count_df$Gene.Name , unique = T)
    working_DEG_count_df$Gene.Name  = NULL
    working_DEG_count_df$UniqueID  = NULL
    
    DEG_bi = working_DEG_count_df2 = working_DEG_count_df 
    
    row.names(DEG_bi ) # row.names(DEG_bi ) =  make.names(DEG_bi$genesymbol, unique=TRUE) # DEG_bi $genesymbol = NULL
    
    DEG_bi[ is.na(DEG_bi)       ] = 0
    DEG_bi[   abs(DEG_bi) >0    ] = 1
    DEG_bi
    
    colnames(DEG_bi) 
    colnames(DEG_bi) = gsub("_logFC" , "", colnames(DEG_bi) ) # colnames(DEG_bi) = gsub("dPSI_","", gsub("_max.abs.dPSI..","",colnames(DEG_bi)))
    #colnames(DEG_bi) = c("DMSO-CHX_vs_DMSO" , "U1-WT100-CHX_vs_DMSO" , "U1-WT100-CHX_vs_DMSO-CHX" , "U1-WT100-CHX_vs_WT100" , "U1-WT100-EtOH_vs_DMSO" , "U1-WT500-CHX_vs_DMSO" , "U1-WT500-CHX_vs_DMSO-CHX" , "U1-WT500-CHX_vs_WT500" , "U1-WT500-EtOH_vs_DMSO" )
    
    head(DEG_bi)
    dim(DEG_bi)
    
    fout = paste0(DEG_plot_out_dir, "binary_",   name_vector_string,".csv")
    write.csv(DEG_bi,file = fout,quote = F)

    ####### total number of DSG ###
    y = colSums(DEG_bi)
    y = as.data.frame(y)
    colnames(y) = "count"
    y
    y  = y[which(y$count < 500), ,drop=FALSE]   
    y
    row.names(y)
    
    ####### UPset ###
    colSums(DEG_bi) 
    colnames(DEG_bi) 


    DEG_bi_working =     DEG_bi[ , colnames(DEG_bi)  %in%  row.names(y)]
    DEG_bi_working
    
    DEG_bi_working =     DEG_bi_working[!(rowSums(    DEG_bi_working) == 0), ]
    DEG_bi_working
    dim(DEG_bi_working)
    head(DEG_bi_working)
    
    
    print(name_vector_string)  #     fout = paste0(DEG_plot_out_dir, "UPset_num_DEG_",name_vector_string,".png")
    fout2 = paste0( DEG_plot_out_dir , "UPset_num_DEG_",name_vector_string,"_2.png")  #     fout = paste0(DEG_plot_out_dir, "UPset_num_DEG_",name_vector_string,".png")
    #colnames_set = colnames(DEG_bi_working)
    colnames_set
    png(fout2, height = 1000, width = 680+ncol(DEG_bi_working)*50)
    #print(  upset(DEG_bi_working ,  sets = colnames_set , mainbar.y.label = "DEG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   )
    print(  upset(DEG_bi_working ,  mainbar.y.label = "DEG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   )
    dev.off()
    
    # colnames_set = colnames(DEG_bi_working)
    # colnames_set
    # head(DEG_bi_working)
    # 
    # y  = y[which(y$count < 500), drop=FALSE]   
    # y
    # row.names(y)
    # 
    # DEG_bi_working =     DEG_bi_working[,     row.names(y)]

    # 
    # png(fout2, height = 1000, width = 680+ncol(DEG_bi_working)*50)
    # print( upset(DEG_bi_working , sets = colnames_set,  mainbar.y.label = "Gene Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))  )
    # dev.off()
    
}
# run.02.DSG_counts_minmaxtable_plot14-2023-0904-UpSetR.R
                              #setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_downstream/")
R.Version()
rm(list = ls())
#install.packages("ggVennDiagram") # 15min
library("ggVennDiagram")
library(ggplot2)
require("UpSetR") # movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header = T, sep = ";")

out_dir = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master_analysis.02.splicing_offtargets_3vs3_direct_comparison_Joon/" #From_DNANenxs_Master_analysis.02.splicing_offtargets_3vs3/
dir.create(out_dir, recursive = T)
setwd(out_dir)

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors #brightjet = colorRampPalette(c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00","#FF3535")) # make gradient of colors

######################
rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
  par(mar = c(15,15,1,1))
  plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n", ylab = "num.DSG",)
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=1) 
}
######################
minpdj_maxdPSI = function(tbtemp, colnn){
  
  ## input : BSSI formatted splice master table, optionally pre-filtered.
  ## output : for each genesymbol, write  minpadj and dPSI for                   == ==  >>  rMATS                 (I)
  ## output : for each genesymbol, write  minpadj and dPSI for                   == ==  >>  Leafcutter           (II)
  ## output : for each genesymbol, write maxPdPSI and dPSI for                   == ==  >>  majiq               (III)
  
  gn = unique(tbtemp$geneSymbol)
  
  if( sum(gn  == ".") > 0 ){ gn = gn[ !(gn  == ".") ] }  # remove un-annotated genes
  
  minmaxcum = as.data.frame( matrix(0, nrow = length(gn) , ncol = 3 ) )
  colnames(minmaxcum) = colnn 
  minmaxcum$genesymbol = gn
  
  for(i in 1:length(gn)){
    #i
    gntemp = gn[i]
    
    tbtemp2 = tbtemp[ which(tbtemp$geneSymbol ==  gntemp) , ]
    
    if(nrow(tbtemp2) ==  1){          # only one row of goi, no "min" needed
      
      dPSItemp = tbtemp2$dPSI
      if(tbtemp2$Algorithm ==  "majiq"){
        padj_PdPSI = tbtemp2$PdPSI    # for majiq, keep PdPSI.
      }else{
        padj_PdPSI = tbtemp2$FDR      # for rmats and leafcutter keep padj
      }
      
    }else{                            # if there are more than two rows of goi (gene of interest) , find either "min" padj or max PdPSI 
      
      rn_minpadj = which( tbtemp2$FDR == min( tbtemp2$FDR, na.rm = T) ) 
      tbtemp3 = tbtemp2[ rn_minpadj, ]
      
      if(nrow(tbtemp3) ==  0){        # if there's no min, it's all majiq, with multiple lines
        rn_maxpdpsi = which(tbtemp2$PdPSI  == max(tbtemp2$PdPSI, na.rm = T)  ) # ..get the max(PdPSI)
        tbtemp3 = tbtemp2[rn_maxpdpsi,]
        
        if(nrow(tbtemp3) ==  1){       # if there's only one max PdPSI
          dPSItemp = tbtemp3$dPSI
          padj_PdPSI = tbtemp3$PdPSI
        }else{                         #if there's more than on max dPSI
          dPSItemp = max(tbtemp3$dPSI) # pick the max(dPSI)
          padj_PdPSI = tbtemp3$PdPSI[1]
        }
        
      }else if(nrow(tbtemp3) ==  1){  # if there's only one min (either rmats or leafcutter, no majiq)
        dPSItemp   = tbtemp3$dPSI
        padj_PdPSI = tbtemp3$FDR
      }else{                          # more than 1 min, pick max(dPSI). both are either rmats or leafcutter. 
        dPSItemp = max(tbtemp3$dPSI)
        padj_PdPSI = tbtemp3$FDR[1]
      }
    }
    minmaxcum[i, 2] = padj_PdPSI
    minmaxcum[i, 3] = dPSItemp
    
  }
  
  return(minmaxcum)
}


din = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master_direct_Joon/"  # 09-02 From_DNANenxs_Master/"
coi_ar = dir(din, pattern = "dPSI0.1")    # all compounds, dPSI>0.1

thr_dPSI_rm = 0.3; 
thr_dPSI_lc = 0.25; 
thr_dPSI_mj = 0.3;

param = paste0("dPSIrm",thr_dPSI_rm,"_dPSIlc",thr_dPSI_lc,"_padj0.05_EdPSI",thr_dPSI_mj,"_PdPSI0.9_M50") # for fout

### minpadj per gene ####################
minmaxall = c()

for(i in 1:length(coi_ar)){
  #  i = 1
  coi = coi_ar[i] # compound of interest
  fin = dir(din,pattern = coi)
  tbtemp = read.table(paste0(din,fin) ,stringsAsFactors = F,sep = ",",header = 1)
  head(tbtemp)
  dim(tbtemp)
  
  ## apply thresholding per algorithm ##
  rm_t = tbtemp[ which( tbtemp$Algorithm ==  "rmats"      & abs(tbtemp$dPSI)       >=  thr_dPSI_rm) , ]
  lc_t = tbtemp[ which( tbtemp$Algorithm ==  "leafcutter" & abs(tbtemp$dPSI)       >=  thr_dPSI_lc) , ]
  mj_t = tbtemp[ which( tbtemp$Algorithm ==  "majiq"      &     tbtemp$Cutoff_dPSI >=  thr_dPSI_mj) , ]
  
  head(rm_t)
  head(lc_t)
  head(mj_t)
  
  tbtemp_t = rbind(rm_t, lc_t, mj_t)
  tbtemp_t = tbtemp_t[!(tbtemp_t$geneSymbol ==  ".") ,] #clean up some names
  head(tbtemp_t)
   dim(tbtemp_t)
   ## 
  # goi = "ATG16L1"; tbtemp_t[which(tbtemp_t$geneSymbol ==  goi) ,]
  
  ## extract min-pval-max-dPSI ##
  minmaxtemp = minpdj_maxdPSI( tbtemp_t, colnn = c("genesymbol",paste0("minpadj_maxPdPSI_",coi) , paste0("dPSI_",coi) ) )  ##  == == == == == == == == == == == == == == == == == == Call the function
  ## minpdj_maxdPSI = function(tbtemp, colnn) # ZG
  #minmaxtemp = minpdj_maxdPSI( tbtemp = tbtemp, colnn = c( "genesymbol" , paste0("minpadj_maxPdPSI_", coi) , paste0("dPSI_", coi) ) ) ##  == == == == == == == == == == == == == == == == == == Call the function
  #head(minmaxtemp)
  
  ## optional: clean up colnames ##
  colnames(minmaxtemp) = gsub("filteredevents_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) ) #   gsub("filteredevents_TST11742_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) )
  
  ## merge ##
  if(i ==  1){ minmaxall = minmaxtemp
  } else{  minmaxall = merge(minmaxall, minmaxtemp, by.x = "genesymbol", by.y = "genesymbol" , all = T) }
}

#minmaxall_back  =  minmaxall minmaxall = minmaxall_back  #=  minmaxall 
head(minmaxall) 
 dim(minmaxall) #[1] 1373  105     2023-09-05    [1] 1673   83

colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_maxP_", colnames(minmaxall)) # colnames(minmaxall) = gsub("filteredevents_","", colnames(minmaxall)) #              gsub("filteredevents_TST11955_","", colnames(minmaxall))
#colnames(minmaxall) = gsub("WT100","WT_100", colnames(minmaxall)) 
#colnames(minmaxall) = gsub("WT500","WT_500", colnames(minmaxall)) 

minmaxall_raw = minmaxall
rownames(minmaxall) = (minmaxall$genesymbol)            # should take care     
#rownames(minmaxall) = make.names(minmaxall$genesymbol , unique = T)            # should take care     
minmaxall$genesymbol = NULL 

colnames(minmaxall) 
row.names(minmaxall) 
minmaxall[1:10, 1:10] 
minmaxall_raw[1:10, 1:10] 

#minmaxall = minmaxall[,order(colnames(minmaxall))]
getwd()

fout = paste0("summarytable_minpadj_maxPdPSI_dPSI_",param,".csv")
write.csv(minmaxall_raw,file = fout,quote = F, row.names = F)

fout = paste0("summarytable_with_gene",param,".csv")
write.csv(minmaxall,file = fout,quote = F, row.names = T)


########################################### 08-01-2023
dim(minmaxall)
minmaxall[1:10, 1:10]

any(row.names(minmaxall) %in% "NA")
any(row.names(minmaxall) %in% "NA.")
any(row.names(minmaxall) %in% "NA..")
minmaxall = minmaxall[ !(row.names(minmaxall) %in% c("NA", "NA.", "NA..")), ]

for(i in 1:ncol(minmaxall) ){
  my_vec = minmaxall[, i]
  print(i)
  print(any( is.infinite(my_vec)) )  #   any( is.infinite(my_vec))          68
}

dim(minmaxall) # 1373 104
head(minmaxall) # 
########################################### 08-01-2023


minmaxall_dPSI  = minmaxall[, grep( "genesymbol|^dPSI", colnames(minmaxall) ) ] #drop padj, only keep dPSI, eq. delete 50% columns. <<============ colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_mP_", colnames(minmaxall) )  #####

head(minmaxall_dPSI )
dim(minmaxall_dPSI )
minmaxall_dPSI[1:5, 1:5]

 dim(minmaxall_dPSI )
colnames(minmaxall_dPSI ) ###############=====================<<<<<<<<<<<<<<< 07-16-2022

dPSI_bi = minmaxall_dPSI2 = minmaxall_dPSI 

row.names(dPSI_bi ) # row.names(dPSI_bi ) =  make.names(dPSI_bi$genesymbol, unique=TRUE) # dPSI_bi $genesymbol = NULL

dPSI_bi[ is.na(dPSI_bi)       ] = 0
dPSI_bi[   abs(dPSI_bi) >0    ] = 1
dPSI_bi

colnames(dPSI_bi) 
colnames(dPSI_bi) = gsub("dPSI_" , "", colnames(dPSI_bi) ) # colnames(dPSI_bi) = gsub("dPSI_","", gsub("_max.abs.dPSI..","",colnames(dPSI_bi)))

head(dPSI_bi)
 dim(dPSI_bi)

fout = paste0("binary_",param,".csv")
write.csv(dPSI_bi,file = fout,quote = F)


##############################
#####  plots : separate file as well ################# 
##############################

## total number of DSG ###

## vetical plot! weew

y = colSums(dPSI_bi)

fout = paste0("binary_numDSG_",param,".png")
png(fout,height = 800, width = length(y)*50) #png(fout,height = 500, width = length(y)*50)
source(pformat)
par(mar = c(30,10,1,1)) # par(mar = c(20,10,1,1))
barplot(y,border = NA,ylim = c(0,max(y)*1.2) ,ylab = "num.DSG",names.arg = names(y) ,las = 2,cex.axis = 1.5, cex.names = 1.5  ) # barplot(y,border = NA,ylim = c(0,max(y)*1.2) ,ylab = "num.DSG",names.arg = names(y) ,las = 2,cex.axis = 1.5, cex.names = 1.5  )
dev.off()

# fout = paste0("UpSet_0_top_5", ".png")
# png(fout,height = 1000, width = length(y)*50) #png(fout,height = 500, width = length(y)*50)
# upset(dPSI_bi ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
# dev.off()

###############################################################################################################
# Exp                Control
####################################
minmaxall_dPSI[1:5, 1:5] 
dPSI_bi[1:5, 1:5]


name_vector_Group1  = c( "DMSO-CHX__vs_DMSO" ,  
                         "U1-WT100-CHX_vs_DMSO" , "U1-WT100-CHX_vs_DMSO-CHX" , "U1-WT100-CHX_vs_WT100" , "U1-WT100-EtOH_vs_DMSO" , 
                         "U1-WT500-CHX_vs_DMSO" , "U1-WT500-CHX_vs_DMSO-CHX" , "U1-WT500-CHX_vs_WT500" , "U1-WT500-EtOH_vs_DMSO" )  




# name_vector_list           = list( name_vector_Group1 , name_vector_Group2 , name_vector_Group3 , name_vector_Group4 , name_vector_Group5 , name_vector_Group6 , name_vector_Group7 ) 
# names(name_vector_list)    = c( "Grp1_DSMO_and_WT" , "Grp2_WT_mut1" , "Grp3_WT_mut2" , "Grp4_Bio5127_and_mut1_mut2" , "Grp5_PTC518_and_mut1_mut2" , "Grp6_Branaplam_and_mut1_mut2" , "Grp7_PTC518_and_Branaplam_and_WT" ) # NO / should be file name, will cause path error
# 
# str(name_vector_list    )

########### ==== &&&& ( 1 ) &&&& ==== #############
########### ==== &&&& ( 2 ) &&&& ==== ############
########### ==== &&&& ( 3 ) &&&& ==== #############
# WT_EtOH (WT_100 and WT_500) --vs-- DMSO_EtOH
###################################################
# DMSO_CHX          DMSO_EtOH                                                 # [13] "DMSO_CHX_vs_DMSO_EtOH" 
# U1_WT_100_EtOH    DMSO_EtOH                                                 # [31] U1_WT_100_EtOH_vs_DMSO_EtOH                   
# U1_WT_500_EtOH    DMSO_EtOH                                                 # [39] U1_WT_500_CHX_vs_DMSO_EtOH   
# WT_CHX (WT_100 and WT_500) --vs-- DMSO_CHX
####==========================================####
###################################################
# U1_WT_100_CHX     DMSO_CHX                                                  # [32] "U1_WT_100_CHX_vs_U1_WT_100_EtOH"             
# U1_WT_500_CHX     DMSO_CHX                                                  # [38] "U1_WT_500_CHX_vs_DMSO_CHX"   
####==========================================####

########### ==== &&&& ( 4 ) &&&& ==== ############# ########### ==== &&&& ( 5 ) &&&& ==== #############
#  mut1_100 
#  mut1_500 
###################################################
U1_mut1_100_EtOH  		U1_WT_100_EtOH     
U1_mut1_100_CHX   		U1_WT_100_EtOH     
U1_mut1_100_CHX   		U1_WT_100_CHX      

U1_mut1_500_EtOH  		U1_WT_500_EtOH     
U1_mut1_500_CHX   		U1_WT_500_EtOH     
U1_mut1_500_CHX   		U1_WT_500_CHX      
####==========================================####
###################################################
head(dPSI_bi)
colnames(dPSI_bi)
dPSI_bi_working  =dPSI_bi[ , c( "U1_mut1_100_EtOH_vs_U1_WT_100_EtOH" , "U1_mut1_100_CHX_vs_U1_WT_100_EtOH" , "U1_mut1_100_CHX_vs_U1_WT_100_CHX" ,  "U1_mut1_500_EtOH_vs_U1_WT_500_EtOH" , "U1_mut1_500_CHX_vs_U1_WT_500_EtOH" , "U1_mut1_500_CHX_vs_U1_WT_500_CHX"    ),  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
dPSI_bi_working

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet2_mut1_vs_WT", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working , sets = colnames_set ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()

########### ==== &&&& ( 6 ) &&&& ==== ############# ########### ==== &&&& ( 7 ) &&&& ==== #############
#  mut2_100 #  mut2_500 
U1_mut2_100_EtOH  		U1_WT_100_EtOH     
U1_mut2_100_CHX   		U1_WT_100_EtOH     
U1_mut2_100_CHX   		U1_WT_100_CHX      
U1_mut2_500_EtOH  		U1_WT_500_EtOH     
U1_mut2_500_CHX   		U1_WT_500_EtOH     
U1_mut2_500_CHX   		U1_WT_500_CHX      
####==========================================####
head(dPSI_bi)
colnames(dPSI_bi)
dPSI_bi_working  =dPSI_bi[ , c( "U1_mut2_100_EtOH_vs_U1_WT_100_EtOH" , "U1_mut2_100_CHX_vs_U1_WT_100_EtOH" , "U1_mut2_100_CHX_vs_U1_WT_100_CHX" ,  "U1_mut2_500_EtOH_vs_U1_WT_500_EtOH" , "U1_mut2_500_CHX_vs_U1_WT_500_EtOH" , "U1_mut2_500_CHX_vs_U1_WT_500_CHX"    ),  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
dPSI_bi_working

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet3_mut2_vs_WT", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()
####==========================================####
########### ==== &&&& ( 8 ) &&&& ==== #############
#SM1-3_EtOH --vs-- WT_EtOH (WT_100 and WT_500) 
###################################################
BIO2195127_EtOH   		U1_WT_100_EtOH     
BIO2195127_EtOH   		U1_WT_500_EtOH     
BIO2195127_CHX    		U1_WT_100_CHX      
BIO2195127_CHX    		U1_WT_500_CHX      
BIO2195127_EtOH   		U1_mut1_100_EtOH   
BIO2195127_EtOH   		U1_mut1_500_EtOH   
BIO2195127_EtOH   		U1_mut2_100_EtOH   
BIO2195127_EtOH   		U1_mut2_500_EtOH   
####==========================================####
head(dPSI_bi)
colnames(dPSI_bi)

dPSI_bi_working  =dPSI_bi[ , c(  "BIO2195127_EtOH_vs_U1_WT_100_EtOH"  ,     "BIO2195127_EtOH_vs_U1_WT_500_EtOH"      , "BIO2195127_CHX_vs_U1_WT_100_CHX"       , "BIO2195127_CHX_vs_U1_WT_500_CHX" ,
                                 "BIO2195127_EtOH_vs_U1_mut1_100_EtOH" , "BIO2195127_EtOH_vs_U1_mut1_500_EtOH" , "BIO2195127_EtOH_vs_U1_mut2_100_EtOH" , "BIO2195127_EtOH_vs_U1_mut2_500_EtOH" )  ,  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
head(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet4_BIO5127", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working , sets = colnames_set ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()


########### ==== &&&& ( 9 ) &&&& ==== #############
#SM1-3_CHX --vs-- WT_CHX (WT_100 and WT_500) 
###################################################
BIO2197294_EtOH   		U1_WT_100_EtOH     
BIO2197294_EtOH   		U1_WT_500_EtOH     
BIO2197294_CHX    		U1_WT_100_CHX      
BIO2197294_CHX    		U1_WT_500_CHX     
BIO2197294_EtOH   		U1_mut1_100_EtOH   
BIO2197294_EtOH   		U1_mut1_500_EtOH   
BIO2197294_EtOH   		U1_mut2_100_EtOH   
BIO2197294_EtOH   		U1_mut2_500_EtOH   
####==========================================####

head(dPSI_bi)
colnames(dPSI_bi)

dPSI_bi_working  =dPSI_bi[ , c(  "BIO2197294_EtOH_vs_U1_WT_100_EtOH"  ,     "BIO2197294_EtOH_vs_U1_WT_500_EtOH"      , "BIO2197294_CHX_vs_U1_WT_100_CHX"       , "BIO2197294_CHX_vs_U1_WT_500_CHX" ,
                                 "BIO2197294_EtOH_vs_U1_mut1_100_EtOH" , "BIO2197294_EtOH_vs_U1_mut1_500_EtOH" , "BIO2197294_EtOH_vs_U1_mut2_100_EtOH" , "BIO2197294_EtOH_vs_U1_mut2_500_EtOH" )  ,  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
head(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet5_PTC518", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working , sets = colnames_set ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()




########### ==== &&&& ( 10 ) &&&& ==== ############
#SM1-3_CHX --vs-- mut1_EtOH (mut1_100 and mut1_500) 
###################################################
Branaplam_EtOH    		U1_WT_100_EtOH     
Branaplam_EtOH    		U1_WT_500_EtOH     
Branaplam_CHX     		U1_WT_100_CHX      
Branaplam_CHX     		U1_WT_500_CHX      
Branaplam_EtOH    		U1_mut1_100_EtOH   
Branaplam_EtOH    		U1_mut1_500_EtOH   
Branaplam_EtOH    		U1_mut2_100_EtOH   
Branaplam_EtOH    		U1_mut2_500_EtOH   
####==========================================####

head(dPSI_bi)
colnames(dPSI_bi)

dPSI_bi_working  =dPSI_bi[ , c(  "Branaplam_EtOH_vs_U1_WT_100_EtOH"  ,     "Branaplam_EtOH_vs_U1_WT_500_EtOH"      , "Branaplam_CHX_vs_U1_WT_100_CHX"       , "Branaplam_CHX_vs_U1_WT_500_CHX" ,
                                 "Branaplam_EtOH_vs_U1_mut1_100_EtOH" , "Branaplam_EtOH_vs_U1_mut1_500_EtOH" , "Branaplam_EtOH_vs_U1_mut2_100_EtOH" , "Branaplam_EtOH_vs_U1_mut2_500_EtOH" )  ,  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
head(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet6_Branaplam", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working , sets = colnames_set ,  mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()

########### ==== &&&& ( 12 ) &&&& ==== #############
# SM1-2_with_WT
####################################################
U1_WT_500_Branaplam   U1_WT_500_EtOH
U1_WT_500_BIO2197294  U1_WT_500_EtOH 

U1_WT_500_Branaplam        U1_WT_500_EtOH    
U1_WT_500_BIO2197294       U1_WT_500_EtOH    
# "U1_WT_500_BIO2197294_vs_BIO2197294_EtOH" 
# "U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH" 
# "U1_WT_500_Branaplam_vs_Branaplam_EtOH"
# "U1_WT_500_Branaplam_vs_U1_WT_500_EtOH"  
####==========================================####
# name_vector_Group7  = c( "U1-WT-500-BIO2197294_vs_DMSO"           ,  "U1-WT-500-BIO2197294_vs_WT-500"        , "U1-WT-500-Branaplam_vs_DMSO"          , "U1-WT-500-Branaplam_vs_WT-500" )
head(dPSI_bi)
colnames(dPSI_bi)
dPSI_bi_working  =dPSI_bi[ , c( "U1_WT_500_BIO2197294_vs_BIO2197294_EtOH",  "U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH" , "U1_WT_500_Branaplam_vs_Branaplam_EtOH" ,  "U1_WT_500_Branaplam_vs_U1_WT_500_EtOH"    ),  drop = F]  
dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
head(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) 

fout = paste0("UpSet_7_WT_with_SM", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working ,  sets = colnames_set , mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()



########### ==== &&&& ( 15 ) &&&& ==== #############
# All combinations comparing CHX to EtOH
####################################################    ## ============================= (0) 
DMSO_CHX                 
U1_WT_100_CHX            U1_WT_100_EtOH
U1_WT_500_CHX            U1_WT_500_EtOH
U1_mut1_100_CHX          U1_mut1_100_EtOH 
U1_mut1_500_CHX          U1_mut1_500_EtOH
U1_mut2_100_CHX          U1_mut2_100_EtOH
U1_mut2_500_CHX          U1_mut2_500_EtOH
BIO2195127_CHX           BIO2195127_EtOH
BIO2197294_CHX           BIO2197294_EtOH
Branaplam_CHX            Branaplam_EtOH

head(dPSI_bi)
colnames(dPSI_bi)
dPSI_bi_working = dPSI_bi[ , c( "DMSO_CHX_vs_DMSO_EtOH" , "U1_WT_100_CHX_vs_DMSO_EtOH"  ,   "U1_WT_500_CHX_vs_DMSO_EtOH" ,  "U1_mut1_100_CHX_vs_U1_mut1_100_EtOH" , "U1_mut1_500_CHX_vs_U1_mut1_500_EtOH" ,  "U1_mut2_100_CHX_vs_U1_mut2_100_EtOH" , "U1_mut2_500_CHX_vs_U1_mut2_500_EtOH" ,  "BIO2195127_CHX_vs_BIO2195127_EtOH", "BIO2197294_CHX_vs_BIO2197294_EtOH" , "Branaplam_CHX_vs_Branaplam_EtOH"   )]  
# 
# dPSI_bi[ , c( "DMSO_CHX_vs_DMSO_EtOH" , "U1_WT_100_CHX_vs_DMSO_EtOH"  ,   "U1_WT_500_CHX_vs_DMSO_EtOH"  )]
# dPSI_bi[ , c( "DMSO_CHX_vs_DMSO_EtOH" , "U1_WT_100_CHX_vs_DMSO_EtOH"  ,   "U1_WT_500_CHX_vs_DMSO_EtOH" ,  "U1_mut1_100_CHX_vs_U1_mut1_100_EtOH" , "U1_mut1_500_CHX_vs_U1_mut1_500_EtOH" )]
# dPSI_bi[ , c( "DMSO_CHX_vs_DMSO_EtOH" , "U1_WT_100_CHX_vs_DMSO_EtOH"  ,   "U1_WT_500_CHX_vs_DMSO_EtOH" ,  "U1_mut1_100_CHX_vs_U1_mut1_100_EtOH" , "U1_mut1_500_CHX_vs_U1_mut1_500_EtOH" , "U1_mut2_100_CHX_vs_U1_mut2_100_EtOH" , "U1_mut2_500_CHX_vs_U1_mut2_500_EtOH" )]

dPSI_bi_working = dPSI_bi_working[!(rowSums(dPSI_bi_working) == 0), ] #
head(dPSI_bi_working)

colnames_set = colnames(dPSI_bi_working) # 15 types #upset(movies , sets = colnames_set )

fout = paste0("UpSet_15", ".png")
png(fout, height = 1000, width = 680+ncol(dPSI_bi_working)*50) 
upset(dPSI_bi_working , sets = colnames_set , mainbar.y.label = "DSG Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))   
dev.off()

#########################################################################################################################################################################################################################
# [1] "BIO2195127_CHX_vs_U1_mut1_100_CHX"      
# [2] "BIO2195127_CHX_vs_U1_mut1_500_CHX"      
# [3] "BIO2195127_CHX_vs_U1_mut2_100_CHX"      
# [4] "BIO2195127_CHX_vs_U1_mut2_500_CHX"      
# [5] "BIO2195127_CHX_vs_U1_WT_100_CHX"        
# [6] "BIO2195127_CHX_vs_U1_WT_500_CHX"        
# [7] "BIO2195127_EtOH_vs_U1_mut1_100_EtOH"    
# [8] "BIO2195127_EtOH_vs_U1_mut1_500_EtOH"    
# [9] "BIO2195127_EtOH_vs_U1_mut2_100_EtOH"    
# [10] "BIO2195127_EtOH_vs_U1_mut2_500_EtOH"    
# [11] "BIO2195127_EtOH_vs_U1_WT_100_EtOH"      
# [12] "BIO2195127_EtOH_vs_U1_WT_500_EtOH"      
# [13] "BIO2197294_CHX_vs_U1_mut1_100_CHX"      
# [14] "BIO2197294_CHX_vs_U1_mut1_500_CHX"      
# [15] "BIO2197294_CHX_vs_U1_mut2_100_CHX"      
# [16] "BIO2197294_CHX_vs_U1_mut2_500_CHX"      
# [17] "BIO2197294_CHX_vs_U1_WT_100_CHX"        
# [18] "BIO2197294_CHX_vs_U1_WT_500_CHX"        
# [19] "BIO2197294_EtOH_vs_U1_mut1_100_EtOH"    
# [20] "BIO2197294_EtOH_vs_U1_mut1_500_EtOH"    
# [21] "BIO2197294_EtOH_vs_U1_mut2_100_EtOH"    
# [22] "BIO2197294_EtOH_vs_U1_mut2_500_EtOH"    
# [23] "BIO2197294_EtOH_vs_U1_WT_100_EtOH"      
# [24] "BIO2197294_EtOH_vs_U1_WT_500_EtOH"      
# [25] "Branaplam_CHX_vs_U1_mut1_100_CHX"       
# [26] "Branaplam_CHX_vs_U1_mut1_500_CHX"       
# [27] "Branaplam_CHX_vs_U1_mut2_100_CHX"       
# [28] "Branaplam_CHX_vs_U1_mut2_500_CHX"       
# [29] "Branaplam_CHX_vs_U1_WT_100_CHX"         
# [30] "Branaplam_CHX_vs_U1_WT_500_CHX"         
# [31] "Branaplam_EtOH_vs_U1_mut1_100_EtOH"     
# [32] "Branaplam_EtOH_vs_U1_mut1_500_EtOH"     
# [33] "Branaplam_EtOH_vs_U1_mut2_100_EtOH"     
# [34] "Branaplam_EtOH_vs_U1_mut2_500_EtOH"     
# [35] "Branaplam_EtOH_vs_U1_WT_100_EtOH"       
# [36] "Branaplam_EtOH_vs_U1_WT_500_EtOH"       

# [37] "U1_mut1_100_CHX_vs_U1_WT_100_CHX"       
# [38] "U1_mut1_100_CHX_vs_U1_WT_100_EtOH"      
# [39] "U1_mut1_100_EtOH_vs_U1_WT_100_EtOH"     

# [40] "U1_mut1_500_CHX_vs_U1_WT_500_CHX"       
# [41] "U1_mut1_500_CHX_vs_U1_WT_500_EtOH"      
# [42] "U1_mut1_500_EtOH_vs_U1_WT_500_EtOH"     

# [43] "U1_mut2_100_CHX_vs_U1_WT_100_CHX"       
# [44] "U1_mut2_100_CHX_vs_U1_WT_100_EtOH"      
# [45] "U1_mut2_100_EtOH_vs_U1_WT_100_EtOH"     

# [46] "U1_mut2_500_CHX_vs_U1_WT_500_CHX"       
# [47] "U1_mut2_500_CHX_vs_U1_WT_500_EtOH"      
# [48] "U1_mut2_500_EtOH_vs_U1_WT_500_EtOH"     

# [49] "U1_WT_500_BIO2197294_vs_BIO2197294_EtOH"
# [50] "U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH" 
# [51] "U1_WT_500_Branaplam_vs_Branaplam_EtOH"  
# [52] "U1_WT_500_Branaplam_vs_U1_WT_500_EtOH" 
