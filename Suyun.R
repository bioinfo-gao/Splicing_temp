R.Version()
rm(list = ls())
#setwd("/home/dhuh/project_RNAseq/TST11955_NGN2_HTT_9compound_profiling/analysis.02.run_splice_analysis_pipeline/")
#setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/analysis.03.splice_analysis/")
#dir.create("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/analysis.02.splicing_offtargets")
#setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/analysis.02.splicing_offtargets")
#setwd("/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/")

#setwd("/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/analysis.02.splicing_offtargets")
# dir.create("/camhpc/ngs/projects/TST11872/dnanexus/Suyun_Branaplam/DSG", recursive = T)
setwd("/camhpc/ngs/projects/TST11872/dnanexus/Suyun_Branaplam/DSG")
pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"

#brightjet = colorRampPalette(c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00","#FF3535")) # make gradient of colors
jet = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors


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



### define inputs ###################
#din = "/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/analysis.02.splicing_offtargets/res.01.filter_events_dPSI_padj/"
#din = "/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/analysis.02.splicing_offtargets/res.01.filter_events_dPSI_padj/iPSC_dPSI_0.1" # 2022-04-05 ZG
# din="/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG
#din="/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/analysis.02.splicing_offtargets/res.02.DSG_counts/" # 2022-04-05 ZG

#din="/camhpc/ngs/projects/TST11955/dnanexus/Suyun_Branaplam/code/" # 2023-02-07 ZG #/" Trailing /  in NEEDED
din = "/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao//code/analysis.02.splicing_offtargets/res.01.filter_events_dPSI_padj/"

#coi_ar = dir(din) # all compounds, dPSI>0.1
coi_ar = dir(din, pattern = "dPSI_0.1")    # sometimes is "dPSI0.1", "_"     # all compounds, dPSI>0.1
#coi_ar = dir( din , pattern = "filteredevents_LMI070_01x_dPSI0.1_padj0.05_M50_EdPSI0.1_PdPSI0.9.csv") # code in leafviz7.R  orkding_dirs = grep( workding_dirs , pattern = 'backup-',invert = TRUE, value = TRUE)   # exclude "backupXXX" folder not needed 
#coi_ar = grep( coi_ar , pattern = ".csv") # code in leafviz7.R  orkding_dirs = grep( workding_dirs , pattern = 'backup-',invert = TRUE, value = TRUE)   # exclude "backupXXX" folder not needed

thr_dPSI_rm = 0.3; 
thr_dPSI_lc = 0.25; 
thr_dPSI_mj = 0.3;

param = paste0("dPSIrm",thr_dPSI_rm,"_dPSIlc",thr_dPSI_lc,"_padj0.05_EdPSI",thr_dPSI_mj,"_PdPSI0.9_M50") # for fout

### minpadj per gene ####################
minmaxall = c()

for(i in 1:length(coi_ar)){
  #i = 1
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
  #minmaxtemp = minpdj_maxdPSI( tbtemp_t, colnn = c("genesymbol",paste0("minpadj_maxPdPSI_",coi) , paste0("dPSI_",coi) ) )  ##  == == == == == == == == == == == == == == == == == == Call the function
  
  ## minpdj_maxdPSI = function(tbtemp, colnn) # ZG
  minmaxtemp = minpdj_maxdPSI( tbtemp = tbtemp, colnn = c( "genesymbol" , paste0("minpadj_maxPdPSI_", coi) , paste0("dPSI_", coi) ) ) ##  == == == == == == == == == == == == == == == == == == Call the function
  head(minmaxtemp)
  
  ## optional: clean up colnames ##
  colnames(minmaxtemp) = gsub("filteredevents_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) ) #   gsub("filteredevents_TST11742_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) )
  
  ## merge ##
  if(i ==  1){ minmaxall = minmaxtemp
  } else{  minmaxall = merge(minmaxall, minmaxtemp, by.x = "genesymbol", by.y = "genesymbol" , all = T) }
}


head(minmaxall) 
 dim(minmaxall) # 539 19
 minmaxall 

# colnames(minmaxall) = gsub("filteredevents_","", colnames(minmaxall)) # gsub("filteredevents_TST11955_","", colnames(minmaxall))
# colnames(minmaxall) = gsub("_BIO",           "", colnames(minmaxall)) # gsub("filteredevents_TST11955_","", colnames(minmaxall))

rownames(minmaxall) = minmaxall$genesymbol;   minmaxall$genesymbol = NULL 
colnames(minmaxall) 
dim(minmaxall) 
minmaxall[1:10, 1:8] 

#### ========= 08-23
 ## format conc #
# conctemp = sapply( strsplit( colnames( minmaxall ) , "-") , "[[", 4)
# conctemp = gsub("uM","000nM",conctemp) #uM to nM
# 
# library(withr); 
# conctemp = sprintf("%05d",as.numeric(gsub("nM","",conctemp))) #add padding zeros
# conctemp = paste0(
#         sapply( strsplit(colnames(minmaxall) ,"-") ,"[[",1) ,"-",
#         sapply( strsplit(colnames(minmaxall) ,"-") ,"[[",2) ,"-",
#         sapply( strsplit(colnames(minmaxall) ,"-") ,"[[",3) ,"-",
#         conctemp,"nM")
# colnames(minmaxall) = conctemp

  ## order ##
#minmaxall = minmaxall[,order(colnames(minmaxall))]

## save filtered table ##
fout = paste0("summarytable_minpadj_maxPdPSI_dPSI_",param,".csv")
write.csv(minmaxall,file = fout,quote = F, row.names = F)
fout = paste0("summarytable_with_gene",param,".csv")
write.csv(minmaxall,file = fout,quote = F, row.names = T) # with gene Name !! here


dim(minmaxall) 
head(minmaxall) 

minmaxall_selected = minmaxall[  ! is.na(minmaxall$minpadj_maxPdPSI_LMI070_01x) ,  ]  # minmaxall[ ( ! is.na(minmaxall$minpadj_maxPdPSI_LMI070_01x) & ! is.na(minmaxall$minpadj_maxPdPSI_LMI070_03x)  ),  ] 

dim(minmaxall_selected )  # 155
head(minmaxall_selected )  
minmaxall_selected[, c("dPSI_LMI070_01x", "dPSI_LMI070_03x", "dPSI_LMI070_10x") ]  

minmaxall_selected1 = minmaxall_selected[  ( ( !is.na(minmaxall_selected$dPSI_LMI070_01x) ) &  ( !is.na(minmaxall_selected$dPSI_LMI070_03x)  ) &  ( !is.na(minmaxall_selected$dPSI_LMI070_10x) ) ) ,  ] 
dim(minmaxall_selected1 )  # 123
head(minmaxall_selected1 )  
minmaxall_selected1[, c("dPSI_LMI070_01x", "dPSI_LMI070_03x", "dPSI_LMI070_10x") ]  

minmaxall_selected2 = minmaxall_selected1[ abs(minmaxall_selected1$dPSI_LMI070_01x) < abs(minmaxall_selected1$dPSI_LMI070_03x) & abs(minmaxall_selected1$dPSI_LMI070_01x) < abs(minmaxall_selected1$dPSI_LMI070_10x) ,  ] 
dim(minmaxall_selected2 )  
head(minmaxall_selected2 )  
minmaxall_selected2  

minmaxall_selected3  = minmaxall_selected2[ order( -abs(minmaxall_selected2$dPSI_LMI070_01x) ),   ]  
fout = paste0("summarytable_branaplam_gene_sorted.csv")
write.csv(minmaxall_selected3, file = fout,quote = F, row.names = T) # with gene Name !! here

# scp zgao1@rstudio-prod.hpc.biogen.com:/camhpc/ngs/projects/TST11955/dnanexus/Suyun_Branaplam/DSG/summarytable_minpadj_maxPdPSI_dPSI_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv /Users/zgao1/Downloads
# The authenticity of host 'rstudio-prod.hpc.biogen.com (10.240.23.201)' can't be established.
# ECDSA key fingerprint is SHA256:79joRXB4XbqQMJDbFJkTTyYaNolLq2yT+t2rHDUmVO8.
# Are you sure you want to continue connecting (yes/no/[fingerprint])? yes
# Warning: Permanently added 'rstudio-prod.hpc.biogen.com,10.240.23.201' (ECDSA) to the list of known hosts.
# zgao1@rstudio-prod.hpc.biogen.com's password: 
#   tput: No value for $TERM and no -T specified
# summarytable_minpadj_maxPdPSI_dPSI_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv                                         100%   44KB 449.2KB/s   00:00

# scp zgao1@rstudio-prod.hpc.biogen.com:/camhpc/ngs/projects/TST11955/dnanexus/Suyun_Branaplam/DSG/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv /Users/zgao1/Downloads
# zgao1@rstudio-prod.hpc.biogen.com's password: 
# tput: No value for $TERM and no -T specified
# summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv                                                      100%   47KB 445.5KB/s   00:00   

dim(minmaxall) 
minmaxall[1:10, 1:8] 

#minmaxall_abs     = abs(minmaxall)     
minmaxall[ is.na(minmaxall)    ] = 0

colnames(minmaxall)     
minmaxall


minmaxall_df     = minmaxall[ ,  c(  "dPSI_RG7916_10x" , "dPSI_RG7916_03x"  ,"dPSI_LMI070_10x" , "dPSI_LMI070_03x"  )]
dim(minmaxall_df  )
head(minmaxall_df  )
head(minmaxall_df ,50 )


minmaxall_df[] <- lapply(minmaxall_df , function(x)  abs( as.numeric( x )  )  )    #rawdata.TST11955[] <- lapply(rawdata.TST11955, function(x) as.numeric( x )) 
minmaxall_df 

sapply(minmaxall_df , class)


minmaxall_7916     = minmaxall_df[ ,  c(  "dPSI_RG7916_10x" , "dPSI_RG7916_03x"  )]
minmaxall_I070     = minmaxall_df[ ,  c(  "dPSI_LMI070_10x" , "dPSI_LMI070_03x"  )]

minmaxall_7916   
minmaxall_I070   
  

######### 
# dPSI_temp = dPSI_selected [rowSums(dPSI_selected)>0,]
# dim(dPSI_temp)
# head(dPSI_temp)
# dPSI_temp = dPSI_temp[ order( rowSums(dPSI_temp), decreasing=T) ,]
# dPSI_temp

minmaxall_7916_wo0 = minmaxall_7916 [rowSums(minmaxall_7916) > 0, ]
dim(minmaxall_7916_wo0)
head(minmaxall_7916_wo0)

minmaxall_7916_wo0_sort = minmaxall_7916_wo0[ order( rowSums(minmaxall_7916_wo0), decreasing=T) ,]
minmaxall_7916_wo0_sort
dim(minmaxall_7916_wo0_sort)


minmaxall_sort = merge(minmaxall_7916_wo0_sort, minmaxall_I070 , by = 0, sort = F) #table_mm10ID

minmaxall_sort_2_columns =  minmaxall_sort[which(minmaxall_sort$dPSI_RG7916_03x > 0 ), ] 
minmaxall_sort_2_columns


head(minmaxall_sort_2_columns  )
dim(minmaxall_sort_2_columns  )

minmaxall_ordered_50 = minmaxall_sort_2_columns[ 1:50 ,  ]

head(minmaxall_ordered_50 )

getwd() # /camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/analysis.02.splicing_offtargets
#fout = paste0("top_50_gene_ordered",param,".csv") # 
fout = paste0("top_50_gene_ordered",param,"_for_location.csv")
write.csv(minmaxall_ordered_50 , file = fout, quote = F, row.names = T)



###########################################
minmaxall_dPSI  = minmaxall[, grep( "genesymbol|^dPSI", colnames(minmaxall) ) ] #drop padj, only keep dPSI, eq. delete 50% columns. <<============
head(minmaxall_dPSI )
 dim(minmaxall_dPSI )
colnames(minmaxall_dPSI ) ###############=====================<<<<<<<<<<<<<<< 07-16-2022

minmaxall_dPSI2 = minmaxall_dPSI

dPSI_bi         = minmaxall_dPSI2                       #ini ### make a binary  table ######################

#dPSI_bi = dPSI_bi[!is.na(dPSI_bi$genesymbol) , ]
# if( sum( is.na(dPSI_bi$genesymbol) )>0){ 
#   dPSI_bi = dPSI_bi[!is.na(dPSI_bi$genesymbol) , ]
#   }
dim(dPSI_bi)
head(dPSI_bi)

dPSI_bi[ is.na(dPSI_bi)    ] = 0
dPSI_bi[   abs(dPSI_bi) >0 ] = 1


colnames(dPSI_bi) 
colnames(dPSI_bi) = gsub("dPSI_" , "", colnames(dPSI_bi) ) # colnames(dPSI_bi) = gsub("dPSI_","", gsub("_max.abs.dPSI..","",colnames(dPSI_bi)))
# colnames(dPSI_bi) = paste0(sapply(strsplit(colnames(dPSI_bi) ,"-") ,"[[",2) ,"-",
#         sapply(strsplit(colnames(dPSI_bi) ,"-") ,"[[",3) ,"-",
#         sapply(strsplit(colnames(dPSI_bi) ,"-") ,"[[",1) ,"-",
#         sapply(strsplit(colnames(dPSI_bi) ,"-") ,"[[",4))
# 
# dPSI_bi = dPSI_bi[,order(colnames(dPSI_bi))]

head(dPSI_bi)
 dim(dPSI_bi)

## write ## #fout = paste0("binary_",param,"x_rmlf_x1_mq.csv")
fout = paste0("binary_",param,".csv")
write.csv(dPSI_bi,file = fout,quote = F)



##############################
#####  plots : separate file as well #################
##############################

## total number of DSG ###
y = colSums(dPSI_bi)

fout = paste0("binary_numDSG_",param,".png")
#png(fout,height = 500, width = length(y)*50)
png(fout,height = 500, width = length(y)*200)
source(pformat)
par(mar = c(20,10,1,1))
barplot(y,border = NA,ylim = c(0,max(y)*1.2) ,ylab = "num.DSG",names.arg = names(y) ,las = 2,cex.axis = 1.5, cex.names = 1.5  )
dev.off()


##########
## HTT ##
dPSI_bi[rownames(dPSI_bi) ==  "HTT",]

# 
# #  ## is lower-dose always a subset of higher dose? #
# allcompound = 1
# 
# if(allcompound ==  1){
#   
#   compounds = c("allcomp")
#   dPSI_temp = dPSI_bi
#   
# }else{
#   
#   temp = colnames(dPSI_bi)[-grep("fibroblast",colnames(dPSI_bi))] #remove "fibroblast"
#   compounds = c("RGX71083","Risdiplam","TEC1" )
#   dPSI_temp = dPSI_bi[ ,grep(comp,colnames(dPSI_bi))]
#   
# }
# 
# 
# for(comp in compounds){
#     dPSI_temp = dPSI_temp[rowSums(dPSI_temp)>0,]
#     dPSI_temp = dPSI_temp[order(rowSums(dPSI_temp) ,decreasing = T) ,]
#   
#     colorr =  c(rgb(1,1,1) ,rgb(1,0,0))
#     fout = paste0("binary_",comp,"_",param,".png")
#     png(fout,height = max(400,nrow(dPSI_temp)*20) , width = 500)
#     source(pformat)
#     par(mar = c(10,3,1,1))
#     hmin = as.matrix(dPSI_temp);
#     colnames(hmin) = gsub(comp,"",colnames(dPSI_temp));
#     rownames(hmin) = rownames(dPSI_temp)
#     heatmap.2(hmin, col = colorr,Rowv = F,Colv = F,main = comp,
#               density.info = "none", trace = "none", #dendrogram = c("row") ,
#               symm = F,symkey = T,symbreaks = T, scale = "none",margins = c(12,12) ,
#               sepwidth = c(0.001,0.001) ,sepcolor = "grey",colsep = 1:ncol(hmin) ,rowsep = 1:nrow(hmin)) #add boundary between slide
#     dev.off()
# }

# 
# ### upset plot ####
# library(UpSetR)
# upsetin = dPSI_bi
# 
# fout = paste0("upset_",exc_status,"_dPSI0.3.png")
# 
# png(fout,height = 1700, width = 4000)
# source(pformat)
# upset(upsetin,order.by = "freq",point.size = 3,line.size = 1,text.scale = 2, nset = ncol(upsetin))
# dev.off()
# 
# 
