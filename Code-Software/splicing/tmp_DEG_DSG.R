# HTT_MSH3_FoxM1_df = rbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
# Stack(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbindlist(list(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2), fill = F) sjmisc::add_rows(), which uses dplyr::bind_rows()
#setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_downstream/")
R.Version()
rm(list = ls())
#install.packages("ggVennDiagram") # 15min
#install.packages("sjmisc")
#install.packages("stacks") library(Stack) library(stacks) #library(sjmisc)
library("ggVennDiagram")
library(ggplot2)
library(data.table)
library(gtools)
# require(devtools)
# install_version("ffbase",  dependencies = T)
# install_version("Stack", version = "2.0.1", repos = "http://cran.us.r-project.org", dependencies = T)
#library(ffbase)
#library(Stack)
out_dir = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master_analysis.02.splicing_offtargets_3vs3/"
dir.create(out_dir, recursive = T)
setwd(out_dir)

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors #brightjet = colorRampPalette(c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00","#FF3535")) # make gradient of colors

######################
# rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
#     par(mar = c(15,15,1,1))
#     plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n", ylab = "num.DSG",)
#     text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=1) 
# }
# ######################
# rbind.all.columns <- function(x, y) {
#     x.diff <- setdiff(colnames(x), colnames(y))
#     y.diff <- setdiff(colnames(y), colnames(x))
#     x[, c(as.character(y.diff))] <- NA
#     y[, c(as.character(x.diff))] <- NA
#     return(rbind(x, y))
# }
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


din = "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master/"
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
dim(minmaxall) #[1] 1673   83

colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_maxP_", colnames(minmaxall)) # colnames(minmaxall) = gsub("filteredevents_","", colnames(minmaxall)) #              gsub("filteredevents_TST11955_","", colnames(minmaxall))
colnames(minmaxall) = gsub("WT100","WT_100", colnames(minmaxall)) 
colnames(minmaxall) = gsub("WT500","WT_500", colnames(minmaxall)) 

minmaxall_raw = minmaxall
rownames(minmaxall) = make.names(minmaxall$genesymbol , unique = T)            # should take care     
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

dim(minmaxall) # 1672 82
head(minmaxall) # 1672 82
########################################### 08-01-2023
#minmaxall  = minmaxall[, grep( "genesymbol|^dPSI", colnames(minmaxall) ) ] #drop padj, only keep dPSI, eq. delete 50% columns. <<============ colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_mP_", colnames(minmaxall) )  #####
minmaxall[1:5, 1:5]

##############################
#####  plots : separate file as well #################
##############################

## total number of DSG ###
y = colSums(dPSI_bi)

fout = paste0("binary_numDSG_",param,".png")
png(fout,height = 800, width = length(y)*50) #png(fout,height = 500, width = length(y)*50)
source(pformat)
par(mar = c(30,10,1,1)) # par(mar = c(20,10,1,1))
barplot(y,border = NA,ylim = c(0,max(y)*1.2) ,ylab = "num.DSG",names.arg = names(y) ,las = 2,cex.axis = 1.5, cex.names = 1.5  ) # barplot(y,border = NA,ylim = c(0,max(y)*1.2) ,ylab = "num.DSG",names.arg = names(y) ,las = 2,cex.axis = 1.5, cex.names = 1.5  )
dev.off()

###########################################################################
load("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/EA20230815_0/TST12145.RData") 
# ls() # "comp_info"       "data_long"       "data_results"    "data_wide"       "MetaData"        "ProteinGeneName" "results_long"   
# dim(results_long)   
head(results_long)   
DEG_full = merge(   ProteinGeneName[ ,c( "UniqueID", "Gene.Name") ] , results_long , sort = F, all.x = T) # 1min
head(DEG_full)
dim(DEG_full)

DMSO_CHX_vs_DMSO_EtOH_dPSI                   = minmaxall[ ,  c("Padj_maxP_DMSO_CHX_vs_DMSO_EtOH"                   , "dPSI_DMSO_CHX_vs_DMSO_EtOH"), drop = F] 
U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_U1_WT_100_EtOH_vs_DMSO_EtOH"             , "dPSI_U1_WT_100_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_U1_WT_500_EtOH_vs_DMSO_EtOH"             , "dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_WT_100_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_U1_WT_100_CHX_vs_DMSO_CHX"               , "dPSI_U1_WT_100_CHX_vs_DMSO_CHX") , drop = F]  
U1_WT_500_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_U1_WT_500_CHX_vs_DMSO_CHX"               , "dPSI_U1_WT_500_CHX_vs_DMSO_CHX") , drop = F] 

U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut1_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut1_100_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut1_500_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut1_500_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_mut1_100_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut1_100_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut1_100_CHX_vs_DMSO_CHX") , drop = F] 
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut1_500_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut1_500_CHX_vs_DMSO_CHX") , drop = F] 
U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut2_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut2_100_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut2_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut2_100_EtOH_vs_DMSO_EtOH") , drop = F] 
U1_mut2_100_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut2_100_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut2_100_CHX_vs_DMSO_CHX") , drop = F] 
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut2_500_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut2_500_CHX_vs_DMSO_CHX") , drop = F] 

BIO2195127_EtOH_vs_DMSO_EtOH_dPSI            = minmaxall[  , c("Padj_maxP_BIO2195127_EtOH_vs_DMSO_EtOH"            , "dPSI_BIO2195127_EtOH_vs_DMSO_EtOH") , drop = F] 
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI            = minmaxall[  , c("Padj_maxP_BIO2197294_EtOH_vs_DMSO_EtOH"            , "dPSI_BIO2197294_EtOH_vs_DMSO_EtOH") , drop = F] 
Branaplam_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_Branaplam_EtOH_vs_DMSO_EtOH"             , "dPSI_Branaplam_EtOH_vs_DMSO_EtOH") , drop = F] 
BIO2195127_CHX_vs_DMSO_CHX_dPSI              = minmaxall[  , c("Padj_maxP_BIO2195127_CHX_vs_DMSO_CHX"              , "dPSI_BIO2195127_CHX_vs_DMSO_CHX") , drop = F] 
BIO2197294_CHX_vs_DMSO_CHX_dPSI              = minmaxall[  , c("Padj_maxP_BIO2197294_CHX_vs_DMSO_CHX"              , "dPSI_BIO2197294_CHX_vs_DMSO_CHX") , drop = F] 
Branaplam_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_Branaplam_CHX_vs_DMSO_CHX"               , "dPSI_Branaplam_CHX_vs_DMSO_CHX") , drop = F] 

U1_WT_500_BIO2197294_vs_DMSO_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_U1_WT_500_BIO2197294_vs_DMSO_EtOH"       , "dPSI_U1_WT_500_BIO2197294_vs_DMSO_EtOH") , drop = F] 
U1_WT_500_Branaplam_vs_DMSO_EtOH_dPSI        = minmaxall[  , c("Padj_maxP_U1_WT_500_Branaplam_vs_DMSO_EtOH"        , "dPSI_U1_WT_500_Branaplam_vs_DMSO_EtOH") , drop = F] 
U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH_dPSI  = minmaxall[  , c("Padj_maxP_U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH"  , "dPSI_U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH") , drop = F] 
U1_WT_500_Branaplam_vs_U1_WT_500_EtOH_dPSI   = minmaxall[  , c("Padj_maxP_U1_WT_500_Branaplam_vs_U1_WT_500_EtOH"   , "dPSI_U1_WT_500_Branaplam_vs_U1_WT_500_EtOH") , drop = F] 

U1_WT_100_CHX_vs_U1_WT_100_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_U1_WT_100_CHX_vs_U1_WT_100_EtOH"         , "dPSI_U1_WT_100_CHX_vs_U1_WT_100_EtOH") , drop = F] 
U1_WT_500_CHX_vs_U1_WT_500_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_U1_WT_500_CHX_vs_U1_WT_500_EtOH"         , "dPSI_U1_WT_500_CHX_vs_U1_WT_500_EtOH") , drop = F] 
U1_mut1_100_CHX_vs_U1_mut1_100_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut1_100_CHX_vs_U1_mut1_100_EtOH"     , "dPSI_U1_mut1_100_CHX_vs_U1_mut1_100_EtOH") , drop = F] 
U1_mut1_500_CHX_vs_U1_mut1_500_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut1_500_CHX_vs_U1_mut1_500_EtOH"     , "dPSI_U1_mut1_500_CHX_vs_U1_mut1_500_EtOH") , drop = F] 
U1_mut2_100_CHX_vs_U1_mut2_100_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut2_100_CHX_vs_U1_mut2_100_EtOH"     , "dPSI_U1_mut2_100_CHX_vs_U1_mut2_100_EtOH") , drop = F] 
U1_mut2_500_CHX_vs_U1_mut2_500_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut2_500_CHX_vs_U1_mut2_500_EtOH"     , "dPSI_U1_mut2_500_CHX_vs_U1_mut2_500_EtOH") , drop = F] 

BIO2195127_CHX_vs_BIO2195127_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_BIO2195127_CHX_vs_BIO2195127_EtOH"       , "dPSI_BIO2195127_CHX_vs_BIO2195127_EtOH") , drop = F] 
BIO2197294_CHX_vs_BIO2197294_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_BIO2197294_CHX_vs_BIO2197294_EtOH"       , "dPSI_BIO2197294_CHX_vs_BIO2197294_EtOH") , drop = F] 
Branaplam_CHX_vs_Branaplam_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_Branaplam_CHX_vs_Branaplam_EtOH"         , "dPSI_Branaplam_CHX_vs_Branaplam_EtOH") , drop = F] 
head(Branaplam_CHX_vs_Branaplam_EtOH_dPSI    )
########### Exp =========== Control ###########################################################################
########### ==== &&&& ( 1 ) &&&& ==== #############         some Changes in HTT
# DMSO_CHX --vs-- DMSO_EtOH
# DMSO_CHX          DMSO_EtOH 
####==========================================####
# DMSO_CHX_vs_DMSO_EtOH_dPSI[which(! is.na( DMSO_CHX_vs_DMSO_EtOH_dPSI$dPSI_DMSO_CHX_vs_DMSO_EtOH) ) , ] 
# U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]

DMSO_CHX_vs_DMSO_EtOH_dPSI   
dim(DMSO_CHX_vs_DMSO_EtOH_dPSI   )
head(DMSO_CHX_vs_DMSO_EtOH_dPSI   )

DMSO_CHX_vs_DMSO_EtOH_dPSI[which( row.names(DMSO_CHX_vs_DMSO_EtOH_dPSI) == "HTT" ) , ] 
DMSO_CHX_vs_DMSO_EtOH_dPSI[which( row.names(DMSO_CHX_vs_DMSO_EtOH_dPSI) == "MSH3" ) , ] 
DMSO_CHX_vs_DMSO_EtOH_dPSI[which( row.names(DMSO_CHX_vs_DMSO_EtOH_dPSI) == "HTT" ) , ] 
DMSO_CHX_vs_DMSO_EtOH_dPSI[which(! is.na( DMSO_CHX_vs_DMSO_EtOH_dPSI$dPSI_DMSO_CHX_vs_DMSO_EtOH) ) , ] 

DEG = DEG_full[DEG_full$test == "DMSO-CHX__vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG

DMSO_CHX_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, DMSO_CHX_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
DMSO_CHX_vs_DMSO_EtOH_dPSI_DEG

HTT_MSH3_FoxM1_df = DMSO_CHX_vs_DMSO_EtOH_dPSI_DEG[ ( DMSO_CHX_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("1_DMSO-CHX_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)

getwd()
########### ==== &&&& ( 2 ) &&&& ==== ############    do.call("smartbind", mget(ls(pattern = "^data.")))
# WT_EtOH (WT_100 and WT_500) --vs-- DMSO_EtOH
# U1_WT_100_EtOH    DMSO_EtOH  
# U1_WT_500_EtOH    DMSO_EtOH 
####==========================================####
U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI # 

U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_100_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]

DEG = DEG_full[DEG_full$test == "U1-WT100-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI_DEG

HTT_MSH3_FoxM1_df1 = U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-WT500-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI_DEG

HTT_MSH3_FoxM1_df2 = U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2

HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_df

HTT_MSH3_FoxM1_fout = paste0("2_U1-WT_EtOHvs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 3 ) &&&& ==== #############           
# WT_CHX (WT_100 and WT_500) --vs-- DMSO_CHX
###################################################
U1_WT_100_CHX     DMSO_CHX                                                  # [32] "U1_WT_100_CHX_vs_U1_WT_100_EtOH"             
U1_WT_500_CHX     DMSO_CHX                                                  # [38] "U1_WT_500_CHX_vs_DMSO_CHX"   
####==========================================####
U1_WT_100_CHX_vs_DMSO_CHX_dPSI  
U1_WT_500_CHX_vs_DMSO_CHX_dPSI  

U1_WT_100_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_WT_100_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_WT_100_CHX_vs_DMSO_CHX) ) , ]
U1_WT_500_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_WT_500_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_WT_500_CHX_vs_DMSO_CHX) ) , ]

DEG = DEG_full[DEG_full$test == "U1-WT100-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_WT_100_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_WT_100_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_WT_100_CHX_vs_DMSO_CHX_dPSI_DEG

HTT_MSH3_FoxM1_df1 = U1_WT_100_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_WT_100_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-WT500-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_WT_500_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_WT_500_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_WT_500_CHX_vs_DMSO_CHX_dPSI_DEG

HTT_MSH3_FoxM1_df2 = U1_WT_500_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_WT_500_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2

HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("3_U1-WT_CHXvs_DMSO-CHX","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)
## Done+++

########### ==== &&&& ( 4 ) &&&& ==== #############     No Changes in HTT
#  mut1_100 
###################################################
U1_mut1_100_EtOH  DMSO_EtOH		U1_WT_100_EtOH     DMSO_EtOH
U1_mut1_100_CHX   DMSO_CHX		U1_WT_100_EtOH     DMSO_EtOH
U1_mut1_100_CHX   DMSO_CHX		U1_mut1_100_EtOH   DMSO_EtOH	
U1_mut1_100_CHX   DMSO_CHX		U1_WT_100_CHX      DMSO_CHX
####==========================================####  
sort(  as.character( unique(DEG_full$test) ) )
# [51] "U1-mut1-100-CHX_vs_DMSO-CHX"            "U1-mut1-100-CHX_vs_mut1-100"           
# [53] "U1-mut1-100-CHX_vs_WT_CHX"              "U1-mut1-100-CHX_vs_WT_EtOH"            
# [55] "U1-mut1-100-EtOH_vs_DMSO"               "U1-mut1-100-EtOH_vs_WT_EtOH"       
head(U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI )
head(U1_mut1_100_CHX_vs_DMSO_CHX_dPSI  )
U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut1_100_EtOH_vs_DMSO_EtOH) ) , ]
U1_mut1_100_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut1_100_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut1_100_CHX_vs_DMSO_CHX) ) , ]


DEG = DEG_full[DEG_full$test == "U1-mut1-100-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-mut1-100-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut1_100_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_mut1_100_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut1_100_CHX_vs_DMSO_CHX_dPSI_DEG
HTT_MSH3_FoxM1_df2 = U1_mut1_100_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_mut1_100_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_dfå
HTT_MSH3_FoxM1_fout = paste0("4_U1-mut1_100_EtOHvs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 5 ) &&&& ==== #############              some Changes in HTT
#  mut1_500 
###################################################
U1_mut1_500_EtOH  DMSO_EtOH		U1_WT_500_EtOH     DMSO_EtOH
U1_mut1_500_CHX   DMSO_CHX		U1_WT_500_EtOH     DMSO_EtOH
U1_mut1_500_CHX   DMSO_CHX		U1_mut1_500_EtOH   DMSO_EtOH	
U1_mut1_500_CHX   DMSO_CHX		U1_WT_500_CHX      DMSO_CHX
####==========================================####   sort(  as.character( unique(DEG_full$test) ) )
U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI 
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI  
sort(  as.character( unique(DEG_full$test) ) )
# [51] "U1-mut1-500-CHX_vs_DMSO-CHX"            "U1-mut1-500-CHX_vs_mut1-500"           
# [53] "U1-mut1-500-CHX_vs_WT_CHX"              "U1-mut1-500-CHX_vs_WT_EtOH"            
# [55] "U1-mut1-500-EtOH_vs_DMSO"               "U1-mut1-500-EtOH_vs_WT_EtOH"       å
head(U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI )
head(U1_mut1_500_CHX_vs_DMSO_CHX_dPSI  )
U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut1_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut1_500_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut1_500_CHX_vs_DMSO_CHX) ) , ]
dim(U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut1_500_EtOH_vs_DMSO_EtOH) ) , ]) # 140
dim(U1_mut1_500_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut1_500_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut1_500_CHX_vs_DMSO_CHX) ) , ])       # 176


DEG = DEG_full[DEG_full$test == "U1-mut1-500-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-mut1-500-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_mut1_500_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI_DEG
HTT_MSH3_FoxM1_df2 = U1_mut1_500_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_mut1_500_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_df

HTT_MSH3_FoxM1_fout = paste0("5_U1-mut1_500_EtOHvs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


# unique(DEG_full$test)
########### ==== &&&& ( 6 ) &&&& ==== ############# No Changes in HTT
#  mut2_100 
###################################################
U1_mut2_100_EtOH  DMSO_EtOH		U1_WT_100_EtOH     DMSO_EtOH
U1_mut2_100_CHX   DMSO_CHX		U1_WT_100_EtOH     DMSO_EtOH
U1_mut2_100_CHX   DMSO_CHX		U1_mut2_100_EtOH   DMSO_EtOH	
U1_mut2_100_CHX   DMSO_CHX		U1_WT_100_CHX      DMSO_CHX
####==========================================####     sort(  as.character( unique(DEG_full$test) ) )
head(U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI )
head(U1_mut2_100_CHX_vs_DMSO_CHX_dPSI  )
U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut2_100_EtOH_vs_DMSO_EtOH) ) , ] # 2 
U1_mut2_100_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut2_100_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut2_100_CHX_vs_DMSO_CHX) ) , ]       # 4


DEG = DEG_full[DEG_full$test == "U1-mut2-100-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-mut2-100-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_100_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_mut2_100_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_100_CHX_vs_DMSO_CHX_dPSI_DEG
HTT_MSH3_FoxM1_df2 = U1_mut2_100_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_mut2_100_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2

HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_df

HTT_MSH3_FoxM1_fout = paste0("6_U1-mut2_100_EtOHvs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 7 ) &&&& ==== #############        No Changes in HTT
#  mut2_500 
###################################################
U1_mut2_500_EtOH  DMSO_EtOH		U1_WT_500_EtOH     DMSO_EtOH
U1_mut2_500_CHX   DMSO_CHX		U1_WT_500_EtOH     DMSO_EtOH
U1_mut2_500_CHX   DMSO_CHX		U1_mut2_500_EtOH   DMSO_EtOH	
U1_mut2_500_CHX   DMSO_CHX		U1_WT_500_CHX      DMSO_CHX
####==========================================####
sort(  as.character( unique(DEG_full$test) ) )
# [51] "U1-mut2-500-CHX_vs_DMSO-CHX"            "U1-mut2-500-CHX_vs_mut2-500"           
# [53] "U1-mut2-500-CHX_vs_WT_CHX"              "U1-mut2-500-CHX_vs_WT_EtOH"            
# [55] "U1-mut2-500-EtOH_vs_DMSO"               "U1-mut2-500-EtOH_vs_WT_EtOH"       
head(U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI )
head(U1_mut2_500_CHX_vs_DMSO_CHX_dPSI  )
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut2_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut2_500_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut2_500_CHX_vs_DMSO_CHX) ) , ]

dim(U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_mut2_500_EtOH_vs_DMSO_EtOH) ) , ]) # 0
dim(U1_mut2_500_CHX_vs_DMSO_CHX_dPSI[which(! is.na( U1_mut2_500_CHX_vs_DMSO_CHX_dPSI$dPSI_U1_mut2_500_CHX_vs_DMSO_CHX) ) , ])       # 104


DEG = DEG_full[DEG_full$test == "U1-mut2-500-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-mut2-500-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_mut2_500_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG
HTT_MSH3_FoxM1_df2 = U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 
HTT_MSH3_FoxM1_df

HTT_MSH3_FoxM1_fout = paste0("7_U1-mut2_500_EtOHvs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 8 ) &&&& ==== #############
#SM1-3_EtOH --vs-- WT_EtOH (WT_100 and WT_500) 
###################################################
BIO2195127_EtOH   DMSO_EtOH		U1_WT_100_EtOH     DMSO_EtOH
BIO2197294_EtOH   DMSO_EtOH		U1_WT_100_EtOH     DMSO_EtOH
Branaplam_EtOH    DMSO_EtOH		U1_WT_100_EtOH     DMSO_EtOH

BIO2195127_EtOH   DMSO_EtOH		U1_WT_500_EtOH     DMSO_EtOH
BIO2197294_EtOH   DMSO_EtOH		U1_WT_500_EtOH     DMSO_EtOH
Branaplam_EtOH    DMSO_EtOH		U1_WT_500_EtOH     DMSO_EtOH
####==========================================####  sort(  as.character( unique(DEG_full$test) ) )
sort(  as.character( unique(DEG_full$test) ) )
# [1] "BIO2195127-CHX_vs_5127"                 "BIO2195127-CHX_vs_DMSO"                
# [3] "BIO2195127-CHX_vs_DMSO-CHX"             "BIO2195127-CHX_vs_U1_U1_mut1_100_CHX"  
# [5] "BIO2195127-CHX_vs_U1_U1_mut1_500_CHX"   "BIO2195127-CHX_vs_U1_U1_mut2_100_CHX"  
# [7] "BIO2195127-CHX_vs_U1_U1_mut2_500_CHX"   "BIO2195127-CHX_vs_U1_WT100_CHX"        
# [9] "BIO2195127-CHX_vs_U1_WT500_CHX"         "BIO2195127-EtOH_vs_DMSO"               
# [11] "BIO2195127-EtOH_vs_U1_U1_mut1_100_EtOH" "BIO2195127-EtOH_vs_U1_U1_mut1_500_EtOH"
# [13] "BIO2195127-EtOH_vs_U1_U1_mut2_100_EtOH" "BIO2195127-EtOH_vs_U1_U1_mut2_500_EtOH"
# [15] "BIO2195127-EtOH_vs_U1_WT100_EtOH"       "BIO2195127-EtOH_vs_U1_WT500_EtOH"      
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI  
Branaplam_EtOH_vs_DMSO_EtOH_dPSI 

BIO2195127_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( BIO2195127_EtOH_vs_DMSO_EtOH_dPSI$dPSI_BIO2195127_EtOH_vs_DMSO_EtOH) ) , ]
dim(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( BIO2195127_EtOH_vs_DMSO_EtOH_dPSI$dPSI_BIO2195127_EtOH_vs_DMSO_EtOH) ) , ]) # 23
dim(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI   )
head(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI   )

BIO2195127_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI) == "HTT" ) , ] 
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI) == "MSH3" ) , ] 
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2195127_EtOH_vs_DMSO_EtOH_dPSI) == "FOXM1" ) , ] 

BIO2197294_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( BIO2197294_EtOH_vs_DMSO_EtOH_dPSI$dPSI_BIO2197294_EtOH_vs_DMSO_EtOH) ) , ]
dim(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( BIO2197294_EtOH_vs_DMSO_EtOH_dPSI$dPSI_BIO2197294_EtOH_vs_DMSO_EtOH) ) , ]) # 17
dim(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI   )
head(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI   )
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI) == "HTT" ) , ] 
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI) == "MSH3" ) , ] 
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(BIO2197294_EtOH_vs_DMSO_EtOH_dPSI) == "FOXM1" ) , ] 

Branaplam_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( Branaplam_EtOH_vs_DMSO_EtOH_dPSI$dPSI_Branaplam_EtOH_vs_DMSO_EtOH) ) , ]
dim(Branaplam_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( Branaplam_EtOH_vs_DMSO_EtOH_dPSI$dPSI_Branaplam_EtOH_vs_DMSO_EtOH) ) , ]) # 23
dim(Branaplam_EtOH_vs_DMSO_EtOH_dPSI   )
head(Branaplam_EtOH_vs_DMSO_EtOH_dPSI   )
Branaplam_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(Branaplam_EtOH_vs_DMSO_EtOH_dPSI) == "HTT" ) , ] 
Branaplam_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(Branaplam_EtOH_vs_DMSO_EtOH_dPSI) == "MSH3" ) , ] 
Branaplam_EtOH_vs_DMSO_EtOH_dPSI[which( row.names(Branaplam_EtOH_vs_DMSO_EtOH_dPSI) == "FOXM1" ) , ] 


DEG = DEG_full[DEG_full$test == "BIO2195127-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, BIO2195127_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = BIO2195127_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( BIO2195127_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "BIO2197294-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, BIO2197294_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = BIO2197294_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( BIO2197294_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


DEG = DEG_full[DEG_full$test == "Branaplam-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
Branaplam_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, Branaplam_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
Branaplam_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = Branaplam_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( Branaplam_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 


HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("8_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)

########### ==== &&&& ( 9 ) &&&& ==== #############
#SM1-3_CHX --vs-- WT_CHX (WT_100 and WT_500) 
###################################################
BIO2195127_CHX    DMSO_CHX		U1_WT_100_CHX      DMSO_CHX
BIO2197294_CHX    DMSO_CHX		U1_WT_100_CHX      DMSO_CHX
Branaplam_CHX     DMSO_CHX		U1_WT_100_CHX      DMSO_CHX

BIO2195127_CHX    DMSO_CHX		U1_WT_500_CHX      DMSO_CHX
BIO2197294_CHX    DMSO_CHX		U1_WT_500_CHX      DMSO_CHX
Branaplam_CHX     DMSO_CHX		U1_WT_500_CHX      DMSO_CHX
####==========================================####
sort(  as.character( unique(DEG_full$test) ) )
BIO2195127_CHX_vs_DMSO_CHX_dPSI 
BIO2197294_CHX_vs_DMSO_CHX_dPSI  
Branaplam_CHX_vs_DMSO_CHX_dPSI  
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]




HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("9_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


unique(DEG_full$test)
########### ==== &&&& ( 10 ) &&&& ==== ############
#SM1-3_CHX --vs-- mut1_EtOH (mut1_100 and mut1_500) 
###################################################
BIO2195127_EtOH   DMSO_EtOH		U1_mut1_100_EtOH   DMSO_EtOH
BIO2197294_EtOH   DMSO_EtOH		U1_mut1_100_EtOH   DMSO_EtOH
Branaplam_EtOH    DMSO_EtOH		U1_mut1_100_EtOH   DMSO_EtOH

BIO2195127_EtOH   DMSO_EtOH		U1_mut1_500_EtOH   DMSO_EtOH
BIO2197294_EtOH   DMSO_EtOH		U1_mut1_500_EtOH   DMSO_EtOH
Branaplam_EtOH    DMSO_EtOH		U1_mut1_500_EtOH   DMSO_EtOH
####==========================================####
BIO2195127_EtOH_vs_DMSO_EtOH_dPSI
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI
Branaplam_EtOH_vs_DMSO_EtOH_dPSI

U1_mut1_100_CHX_vs_DMSO_CHX_dPSI 
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI 

U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]





HTT_MSH3_FoxM1_df = _vs_DMSO_EtOH_dPSI_DEG[ ( _vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df



HTT_MSH3_FoxM1_fout = paste0("10_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 11 ) &&&& ==== #############
#SM1-3_EtOH --vs-- mut1_CHX (mut1_100 and mut1_500) 
####################################################
BIO2195127_CHX    DMSO_CHX		U1_mut1_100_CHX      DMSO_CHX
BIO2197294_CHX    DMSO_CHX		U1_mut1_100_CHX      DMSO_CHX
Branaplam_CHX     DMSO_CHX		U1_mut1_100_CHX      DMSO_CHX

BIO2195127_CHX    DMSO_CHX		U1_mut1_500_CHX      DMSO_CHX
BIO2197294_CHX    DMSO_CHX		U1_mut1_500_CHX      DMSO_CHX
Branaplam_CHX     DMSO_CHX		U1_mut1_500_CHX      DMSO_CHX
####==========================================####
BIO2195127_CHX_vs_DMSO_CHX_dPSI
BIO2197294_CHX_vs_DMSO_CHX_dPSI
Branaplam_CHX_vs_DMSO_CHX_dPSI

U1_mut1_100_CHX_vs_DMSO_CHX_dPSI 
U1_mut1_500_CHX_vs_DMSO_CHX_dPSI 

U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]

DEG = DEG_full[DEG_full$test == "",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG




HTT_MSH3_FoxM1_df

HTT_MSH3_FoxM1_fout = paste0("11_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)

########### ==== &&&& ( 12 ) &&&& ==== #############
# SM1-2_with_WT
####################################################
U1_WT_500_Branaplam   DMSO_EtOH
U1_WT_500_BIO2197294  DMSO_EtOH 

U1_WT_500_Branaplam   U1_WT_500_EtOH
U1_WT_500_BIO2197294  U1_WT_500_EtOH 

U1_WT_500_Branaplam   DMSO_EtOH     U1_WT_500_EtOH    DMSO_EtOH
U1_WT_500_BIO2197294  DMSO_EtOH     U1_WT_500_EtOH    DMSO_EtOH
####==========================================####

U1_WT_500_BIO2197294_vs_DMSO_EtOH_dPSI  
U1_WT_500_Branaplam_vs_DMSO_EtOH_dPSI  
U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH_dPSI 
U1_WT_500_Branaplam_vs_U1_WT_500_EtOH_dPSI  

U1_WT_500_BIO2197294_vs_DMSO_EtOH_dPSI
U1_WT_500_Branaplam_vs_DMSO_EtOH_dPSI

U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI 

U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]






HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("12_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)



########### ==== &&&& ( 13 ) &&&& ==== ############
#SM1-3_EtOH --vs-- mut2_EtOH (mut2_100 and mut2_500) 
###################################################
BIO2195127_EtOH   DMSO_EtOH		U1_mut2_100_EtOH   DMSO_EtOH        # (13.1) SM1 vs mut2 in EtOH (mut2-100)
BIO2197294_EtOH   DMSO_EtOH		U1_mut2_100_EtOH   DMSO_EtOH        # (13.2) SM2 vs mut2 in EtOH
Branaplam_EtOH    DMSO_EtOH		U1_mut2_100_EtOH   DMSO_EtOH        # (13.3) SM3 vs mut2 in EtOH

BIO2195127_EtOH   DMSO_EtOH		U1_mut2_500_EtOH   DMSO_EtOH
BIO2197294_EtOH   DMSO_EtOH		U1_mut2_500_EtOH   DMSO_EtOH
Branaplam_EtOH    DMSO_EtOH		U1_mut2_500_EtOH   DMSO_EtOH
####==========================================####

BIO2195127_EtOH_vs_DMSO_EtOH_dPSI
BIO2197294_EtOH_vs_DMSO_EtOH_dPSI
BIO2197294_EtOH_vs_DMSO_EtOH_genes

U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI 


U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]



HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("13_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)


########### ==== &&&& ( 14 ) &&&& ==== #############
#SM1-3_CHX --vs-- mut2_CHX (mut2_100 and mut2_500) 
####################################################
BIO2195127_CHX    DMSO_CHX		U1_mut2_100_CHX      DMSO_CHX      # (14.1) SM1 vs mut2 in CHX (mut2-100)
BIO2197294_CHX    DMSO_CHX		U1_mut2_100_CHX      DMSO_CHX      # (14.2) SM2 vs mut2 in CHX
Branaplam_CHX     DMSO_CHX		U1_mut2_100_CHX      DMSO_CHX      # (14.3) SM3 vs mut2 in CHX

BIO2195127_CHX    DMSO_CHX		U1_mut2_500_CHX      DMSO_CHX
BIO2197294_CHX    DMSO_CHX		U1_mut2_500_CHX      DMSO_CHX
Branaplam_CHX     DMSO_CHX		U1_mut2_500_CHX      DMSO_CHX
####==========================================####
BIO2195127_CHX_vs_DMSO_CHX_dPSI
BIO2197294_CHX_vs_DMSO_CHX_dPSI
Branaplam_CHX_vs_DMSO_CHX_dPSI


HTT_MSH3_FoxM1_df



HTT_MSH3_FoxM1_fout = paste0("14_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)



########### ==== &&&& ( 15 ) &&&& ==== #############
# All combinations comparing CHX to EtOH
####################################################    ## ============================= (0) 
DMSO_CHX                 DMSO_EtOH
U1_WT_100_CHX            U1_WT_100_EtOH
U1_WT_500_CHX            U1_WT_500_EtOH
U1_mut1_100_CHX          U1_mut1_100_EtOH 
U1_mut1_500_CHX          U1_mut1_500_EtOH
U1_mut2_100_CHX          U1_mut2_100_EtOH
U1_mut2_500_CHX          U1_mut2_500_EtOH
BIO2195127_CHX           BIO2195127_EtOH
BIO2197294_CHX           BIO2197294_EtOH
Branaplam_CHX            Branaplam_EtOH
####==========================================####

U1_WT_100_CHX_vs_U1_WT_100_EtOH_dPSI  
U1_WT_500_CHX_vs_U1_WT_500_EtOH_dPSI 

U1_mut1_100_CHX_vs_U1_mut1_100_EtOH_dPSI

U1_mut2_100_CHX_vs_U1_mut2_100_EtOH_dPSI
U1_mut2_500_CHX_vs_U1_mut2_500_EtOH_dPSI 
BIO2195127_CHX_vs_BIO2195127_EtOH_dPSI
BIO2197294_CHX_vs_BIO2197294_EtOH_dPSI  
Branaplam_CHX_vs_Branaplam_EtOH_dPSI 


U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]
U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI[which(! is.na( U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI$dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH) ) , ]


DEG = DEG_full[DEG_full$test == "U1-mut2-500-EtOH_vs_DMSO",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG = merge( DEG, U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG
HTT_MSH3_FoxM1_df1 = U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG[ ( U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df1

DEG = DEG_full[DEG_full$test == "U1-mut2-500-CHX_vs_DMSO-CHX",  c("UniqueID" , "Gene.Name", "test" , "logFC" , "P.Value" , "Adj.P.Value"), ] 
DEG
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG = merge( DEG, U1_mut2_500_CHX_vs_DMSO_CHX_dPSI, by.x = "Gene.Name", by.y = 0,  all.x = T)     
U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG
HTT_MSH3_FoxM1_df2 = U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG[ ( U1_mut2_500_CHX_vs_DMSO_CHX_dPSI_DEG$Gene.Name  %in%  c( "HTT" , "MSH3" ,"FOXM1" ) ) ,    ]
HTT_MSH3_FoxM1_df2


HTT_MSH3_FoxM1_df =  smartbind(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) # rbind.all.columns(HTT_MSH3_FoxM1_df1, HTT_MSH3_FoxM1_df2) 



HTT_MSH3_FoxM1_df
HTT_MSH3_FoxM1_fout = paste0("15_vs_DMSO","_for_HTT_MSH3_FoxM1.csv")
write.csv(HTT_MSH3_FoxM1_df, file=HTT_MSH3_FoxM1_fout, quote=F, row.names=T)




####==========================================####

#########################################################################################################################################################################################################################
# DMSO_CHX_vs_DMSO_EtOH_dPSI                   = minmaxall[ ,  c("Padj_maxP_DMSO_CHX_vs_DMSO_EtOH"                   , "dPSI_DMSO_CHX_vs_DMSO_EtOH"), drop = F] 
# 
# U1_WT_100_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_U1_WT_100_EtOH_vs_DMSO_EtOH"             , "dPSI_U1_WT_100_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_WT_500_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_U1_WT_500_EtOH_vs_DMSO_EtOH"             , "dPSI_U1_WT_500_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_WT_100_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_U1_WT_100_CHX_vs_DMSO_CHX"               , "dPSI_U1_WT_100_CHX_vs_DMSO_CHX") , drop = F]  
# U1_WT_500_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_U1_WT_500_CHX_vs_DMSO_CHX"               , "dPSI_U1_WT_500_CHX_vs_DMSO_CHX") , drop = F] 
# 
# U1_mut1_100_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut1_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut1_100_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_mut1_500_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut1_500_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut1_500_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_mut1_100_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut1_100_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut1_100_CHX_vs_DMSO_CHX") , drop = F] 
# U1_mut1_500_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut1_500_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut1_500_CHX_vs_DMSO_CHX") , drop = F] 
# 
# U1_mut2_100_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut2_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut2_100_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_mut2_500_EtOH_vs_DMSO_EtOH_dPSI           = minmaxall[  , c("Padj_maxP_U1_mut2_100_EtOH_vs_DMSO_EtOH"           , "dPSI_U1_mut2_100_EtOH_vs_DMSO_EtOH") , drop = F] 
# U1_mut2_100_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut2_100_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut2_100_CHX_vs_DMSO_CHX") , drop = F] 
# U1_mut2_500_CHX_vs_DMSO_CHX_dPSI             = minmaxall[  , c("Padj_maxP_U1_mut2_500_CHX_vs_DMSO_CHX"             , "dPSI_U1_mut2_500_CHX_vs_DMSO_CHX") , drop = F] 
# 
# BIO2195127_EtOH_vs_DMSO_EtOH_dPSI            = minmaxall[  , c("Padj_maxP_BIO2195127_EtOH_vs_DMSO_EtOH"            , "dPSI_BIO2195127_EtOH_vs_DMSO_EtOH") , drop = F] 
# BIO2197294_EtOH_vs_DMSO_EtOH_dPSI            = minmaxall[  , c("Padj_maxP_BIO2197294_EtOH_vs_DMSO_EtOH"            , "dPSI_BIO2197294_EtOH_vs_DMSO_EtOH") , drop = F] 
# Branaplam_EtOH_vs_DMSO_EtOH_dPSI             = minmaxall[  , c("Padj_maxP_Branaplam_EtOH_vs_DMSO_EtOH"             , "dPSI_Branaplam_EtOH_vs_DMSO_EtOH") , drop = F] 
# 
# BIO2195127_CHX_vs_DMSO_CHX_dPSI              = minmaxall[  , c("Padj_maxP_BIO2195127_CHX_vs_DMSO_CHX"              , "dPSI_BIO2195127_CHX_vs_DMSO_CHX") , drop = F] 
# BIO2197294_CHX_vs_DMSO_CHX_dPSI              = minmaxall[  , c("Padj_maxP_BIO2197294_CHX_vs_DMSO_CHX"              , "dPSI_BIO2197294_CHX_vs_DMSO_CHX") , drop = F] 
# Branaplam_CHX_vs_DMSO_CHX_dPSI               = minmaxall[  , c("Padj_maxP_Branaplam_CHX_vs_DMSO_CHX"               , "dPSI_Branaplam_CHX_vs_DMSO_CHX") , drop = F] 
# 
# U1_WT_500_BIO2197294_vs_DMSO_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_U1_WT_500_BIO2197294_vs_DMSO_EtOH"       , "dPSI_U1_WT_500_BIO2197294_vs_DMSO_EtOH") , drop = F] 
# U1_WT_500_Branaplam_vs_DMSO_EtOH_dPSI        = minmaxall[  , c("Padj_maxP_U1_WT_500_Branaplam_vs_DMSO_EtOH"        , "dPSI_U1_WT_500_Branaplam_vs_DMSO_EtOH") , drop = F] 
# U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH_dPSI  = minmaxall[  , c("Padj_maxP_U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH"  , "dPSI_U1_WT_500_BIO2197294_vs_U1_WT_500_EtOH") , drop = F] 
# U1_WT_500_Branaplam_vs_U1_WT_500_EtOH_dPSI   = minmaxall[  , c("Padj_maxP_U1_WT_500_Branaplam_vs_U1_WT_500_EtOH"   , "dPSI_U1_WT_500_Branaplam_vs_U1_WT_500_EtOH") , drop = F] 
# 
# U1_WT_100_CHX_vs_U1_WT_100_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_U1_WT_100_CHX_vs_U1_WT_100_EtOH"         , "dPSI_U1_WT_100_CHX_vs_U1_WT_100_EtOH") , drop = F] 
# U1_WT_500_CHX_vs_U1_WT_500_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_U1_WT_500_CHX_vs_U1_WT_500_EtOH"         , "dPSI_U1_WT_500_CHX_vs_U1_WT_500_EtOH") , drop = F] 
# U1_mut1_100_CHX_vs_U1_mut1_100_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut1_100_CHX_vs_U1_mut1_100_EtOH"     , "dPSI_U1_mut1_100_CHX_vs_U1_mut1_100_EtOH") , drop = F] 
# U1_mut1_500_CHX_vs_U1_mut1_500_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut1_500_CHX_vs_U1_mut1_500_EtOH"     , "dPSI_U1_mut1_500_CHX_vs_U1_mut1_500_EtOH") , drop = F] 
# U1_mut2_100_CHX_vs_U1_mut2_100_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut2_100_CHX_vs_U1_mut2_100_EtOH"     , "dPSI_U1_mut2_100_CHX_vs_U1_mut2_100_EtOH") , drop = F] 
# U1_mut2_500_CHX_vs_U1_mut2_500_EtOH_dPSI     = minmaxall[  , c("Padj_maxP_U1_mut2_500_CHX_vs_U1_mut2_500_EtOH"     , "dPSI_U1_mut2_500_CHX_vs_U1_mut2_500_EtOH") , drop = F] 
# 
# BIO2195127_CHX_vs_BIO2195127_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_BIO2195127_CHX_vs_BIO2195127_EtOH"       , "dPSI_BIO2195127_CHX_vs_BIO2195127_EtOH") , drop = F] 
# BIO2197294_CHX_vs_BIO2197294_EtOH_dPSI       = minmaxall[  , c("Padj_maxP_BIO2197294_CHX_vs_BIO2197294_EtOH"       , "dPSI_BIO2197294_CHX_vs_BIO2197294_EtOH") , drop = F] 
# Branaplam_CHX_vs_Branaplam_EtOH_dPSI         = minmaxall[  , c("Padj_maxP_Branaplam_CHX_vs_Branaplam_EtOH"         , "dPSI_Branaplam_CHX_vs_Branaplam_EtOH") , drop = F] 
# 
