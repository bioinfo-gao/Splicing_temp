library(tidyverse)

data_1 <- structure( list( GO.term = c("BP", "MF", "CC", "KEGG"),          
                           Count = c(163L,  48L, 58L, 27L), 
                           Enrichment = c(0.008, 0.007, 0.008, 0.008),
                           P.value = c(0.37,   0.33, 0.39, 0.43)    ),
                     class = "data.frame",
                     row.names = c(NA, 4L)
                     )
data_2 <- structure( list( GO.term = c("BP", "MF", "CC", "KEGG"),
                           Count = c(167L, 50L, 50L, 23L), 
                           Enrichment = c(0.01, 0.008, 0.006, 0.01),
                           P.value = c(0.31, 0.29, 0.34, 0.37 ) ),
                     class = "data.frame",
                     row.names = c(NA, 4L)
                     )

data_3 <- structure( list( GO.term = c("BP", "MF", "CC", "KEGG"),
                           Count = c(123L, 44L, 50L, 14L),
                           Enrichment = c(0.009, 0.01, 0.007, 0.009),
                           P.value = c(0.22, 0.22, 0.24, 0.28) ) ,
                     class = "data.frame",
                     row.names = c(NA, 4L)
                     )
data_1
data_2
data_3
bind_rows('data1' = data_1,
          'data2' = data_2,
          'data3' = data_3,  .id = 'type') %>%
    ggplot( aes( x = Enrichment , 
                 y = GO.term, 
                 col = P.value, 
                 size = Count, 
                 shape = type ) ) + 
    geom_point()


#setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_downstream/")
R.Version()
#rm(list = ls()) #install.packages("ggVennDiagram") # 15min
library("ggVennDiagram")
library(ggplot2)
require("UpSetR") # movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header = T, sep = ";")
library(dplyr)
library(tidyverse)
# din     = "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/From_DNANenxs_Master/"
# out_dir = "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/splicing_offtargets_distrtibution/"

din     = "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/From_DNANenxs_Master/"
out_dir = "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/splicing_offtargets_distrtibution/"

dir.create(out_dir, recursive = T)
setwd(out_dir)

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors #brightjet = colorRampPalette(c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00","#FF3535")) # make gradient of colors

######################
minpdj_maxdPSI_detail = function(tbtemp, colnn,  Algorithm ){
    
    gn = unique(tbtemp$geneSymbol)
    if( sum(gn  == ".") > 0 ){ gn = gn[ !(gn  == ".") ] }  # remove un-annotated genes
    
    minmaxcum = as.data.frame( matrix(0, nrow = length(gn) , ncol = 4 ) )
    colnames(minmaxcum) = colnn 
    minmaxcum$genesymbol = gn
    minmaxcum$Algorithm = Algorithm
    minmaxcum
    
    for(i in 1:length(gn) ){
        #             i =1 
        i
        #i = i+1
        gntemp = gn[i]
        
        tbtemp2 = tbtemp[ which(tbtemp$geneSymbol ==  gntemp) , ]
        
        print(i)
        print( tbtemp2)
        
        tbtemp3 = tbtemp2[which( tbtemp2$FDR_equal == min( tbtemp2$FDR_equal, na.rm = T) ), ]  # for leafcutter and rmats
        tbtemp3
        
        #     dPSItemp = max(tbtemp3$abs_dPSI)  
        #     padj_PdPSI = tbtemp3$FDR_equal[1]
        #     }
        # }
        
        minmaxcum[i, 3] = tbtemp3$FDR_equal[1]
        minmaxcum[i, 4] = max(tbtemp3$abs_dPSI)  
        
    }
    
    return(minmaxcum)
}

coi_ar = dir(din, pattern = "dPSI0.1")    # all compounds, dPSI>0.1 #coi_ar = dir(din)

thr_dPSI_rm = 0.3; 
thr_dPSI_lc = 0.25; 
thr_dPSI_mj = 0.3;

param = paste0("dPSIrm",thr_dPSI_rm,"_dPSIlc",thr_dPSI_lc,"_padj0.05_EdPSI",thr_dPSI_mj,"_PdPSI0.9_M50") # for fout

### minpadj per gene ####################
minmaxall = c()

for(i in 1:length(coi_ar)){
    #    i = 1
    coi = coi_ar[i] # compound of interest
    fin = dir(din , pattern = coi)
    tbtemp = read.table(paste0(din,fin) ,stringsAsFactors = F,sep = ",",header = 1)
    head(tbtemp)
    dim(tbtemp)
    
    tbtemp$abs_dPSI  =  abs(tbtemp$dPSI)
    tbtemp$FDR_equal =      tbtemp$FDR    # as default
    tbtemp[which(tbtemp$Algorithm == "majiq_v2"), ]$FDR_equal =  1- tbtemp[which(tbtemp$Algorithm == "majiq_v2"), ]$PdPSI   #   tbtemp2$FDR_equal = (1-tbtemp2$PdPSI) 
    
    ## apply thresholding per algorithm ##
    rm_t = tbtemp[ which( tbtemp$Algorithm ==  "rmats_turbo" & abs(tbtemp$dPSI)       >=  thr_dPSI_rm) , ] #   rm_t = tbtemp[ which( tbtemp$Algorithm ==  "rmats"      & abs(tbtemp$dPSI)       >=  thr_dPSI_rm) , ] 
    lc_t = tbtemp[ which( tbtemp$Algorithm ==  "leafcutter"  & abs(tbtemp$dPSI)       >=  thr_dPSI_lc) , ]
    mj_t = tbtemp[ which( tbtemp$Algorithm ==  "majiq_v2"    &     tbtemp$Cutoff_dPSI >=  thr_dPSI_mj) , ] #   mj_t = tbtemp[ which( tbtemp$Algorithm ==  "majiq"      &     tbtemp$Cutoff_dPSI >=  thr_dPSI_mj) , ]
    
    head(rm_t)
    head(lc_t)
    head(mj_t) #   tail(mj_t)
    
    rm_r = rm_t[, c( "geneSymbol",  "GeneID", "Algorithm", "FDR", "PdPSI", "FDR_equal", "dPSI", "abs_dPSI")] 
    lc_r = lc_t[, c( "geneSymbol",  "GeneID", "Algorithm", "FDR", "PdPSI", "FDR_equal", "dPSI", "abs_dPSI")] 
    mj_r = mj_t[, c( "geneSymbol",  "GeneID", "Algorithm", "FDR", "PdPSI", "FDR_equal", "dPSI", "abs_dPSI")] #   tbtemp_t[, c( "chr", "strand", "geneSymbol",  "GeneID", "Event", "Algorithm", "FDR", "PValue", "dPSI", "M", "PSI_1", "PSI_2")] # 199 18
    
    head(rm_r)
    head(lc_r)
    head(mj_r) #   tail(mj_t)
    
    coi_name = gsub( "filteredevents_" ,  "" , coi )
    coi_name = gsub( "_dPSI0.1_padj0.05_M50_EdPSI0.1_PdPSI0.9.csv" ,  "" , coi_name )
    coi_name = gsub( "BIO_" ,  "" , coi_name )
    coi_name
    
    rm_r_temp = minpdj_maxdPSI_detail( tbtemp = rm_r , colnn = c("genesymbol" , "Algorithm", paste0("minpadj_maxPdPSI_",coi_name) , paste0("dPSI_",coi_name) ) , Algorithm = "rmats_turbo" )  ##  == == == == == == == == == == == == == == == == == == Call the function
    lc_r_temp = minpdj_maxdPSI_detail( tbtemp = lc_r , colnn = c("genesymbol" , "Algorithm", paste0("minpadj_maxPdPSI_",coi_name) , paste0("dPSI_",coi_name) ) , Algorithm = "leafcutter"  )  ##  == == == == == == == == == == == == == == == == == == Call the function
    mj_r_temp = minpdj_maxdPSI_detail( tbtemp = mj_r , colnn = c("genesymbol" , "Algorithm", paste0("minpadj_maxPdPSI_",coi_name) , paste0("dPSI_",coi_name) ) , Algorithm = "majiq_v2"    )  ##  == == == == == == == == == == == == == == == == == == Call the function
    
    rm_r_temp
    lc_r_temp
    mj_r_temp
    dim( rm_r_temp) 
    dim(lc_r_temp)
    dim(mj_r_temp)
    
    
    DSG_Called = bind_rows('rmats_turbo' = rm_r_temp,
                           'leafcutter'  = lc_r_temp,
                           'majiq_v2'    = mj_r_temp,  .id = 'type') 
    
    DSG_Called %>%   ggplot( aes( x = dPSI_SH_1949634_10x_vs_SH_DMSO, 
                                  y = log( minpadj_maxPdPSI_SH_1949634_10x_vs_SH_DMSO ), 
                                  color = type ,
                                  shape = type ) ) +     
        geom_point() + ggtitle( paste0( "DSG genes in     " , coi_name ) ) + theme(plot.title = element_text(hjust = 0.5))
    
    bind_rows('rmats_turbo' = rm_r_temp,
              'leafcutter'  = lc_r_temp,
              'majiq_v2'    = mj_r_temp,  .id = 'type') %>%
        ggplot( aes( x = dPSI_SH_1949634_10x_vs_SH_DMSO, 
                     y = log( minpadj_maxPdPSI_SH_1949634_10x_vs_SH_DMSO ), 
                     color = type ,
                     shape = type ) ) +     
        geom_point() + ggtitle( paste0( "DSG genes in     " , coi_name ) ) + theme(plot.title = element_text(hjust = 0.5))      #+       #ylim( 0, 0.001)
    
    
    ggplot(DSG_Called, aes(fill = type, y= dPSI_SH_1949634_10x_vs_SH_DMSO, x = type)) +
        geom_bar(position='dodge', stat='identity')
    
    DSG_Called
    #p<-ggplot(df, aes(x=Category, y=Mean, fill=Quality)) + 
    ggplot(DSG_Called, aes(fill = type, y= mean(dPSI_SH_1949634_10x_vs_SH_DMSO), x = type)) +
        geom_point() #+
    geom_errorbar(aes(ymin= mean() -sd, ymax=Mean+sd), width=.2,
                  position=position_dodge(0.05))
    
    DSG_Called_means <- DSG_Called %>% 
        group_by(type) %>% 
        summarize(dPSI_mean = mean(dPSI_SH_1949634_10x_vs_SH_DMSO) ,
                  dPSI_sd   = sd(dPSI_SH_1949634_10x_vs_SH_DMSO)                    ) 
    
    DSG_Called_means
    
    # Creating barplots of means
    ggplot(DSG_Called_means, aes(x = type, y = dPSI_mean)) +
        geom_point() +      
        geom_errorbar(aes( ymin= (dPSI_mean - dPSI_sd), ymax= (dPSI_mean + dPSI_sd) , width=.2 ) ) #          geom_errorbar(aes( ymin= (dPSI_mean - dPSI_sd), ymax= (dPSI_mean + dPSI_sd) , width=.2,  position = position_dodge( 0.05 ) ) )
    
    ggplot(DSG_Called, aes(x = type, y = dPSI_SH_1949634_10x_vs_SH_DMSO)) +
        geom_boxplot() +
        geom_jitter()
    #geom_errorbar(aes( ymin= (dPSI_mean - dPSI_sd), ymax= (dPSI_mean + dPSI_sd) , width=.2 ) ) #          geom_errorbar(aes( ymin= (dPSI_mean - dPSI_sd), ymax= (dPSI_mean + dPSI_sd) , width=.2,  position = position_dodge( 0.05 ) ) )
    
    #geom_bar(stat="identity")      
    
    tbtemp_t = rbind(rm_t, lc_t, mj_t)
    tbtemp_t = tbtemp_t[!(tbtemp_t$geneSymbol ==  ".") ,] #clean up some names
    head(tbtemp_t)
    dim(tbtemp_t) # 199 18
    
    # goi = "SMN2"; tbtemp_t[which(tbtemp_t$geneSymbol ==  goi) ,]   # goi = "ATG16L1"; tbtemp_t[which(tbtemp_t$geneSymbol ==  goi) ,]
    
    ## extract min-pval-max-dPSI ##
    #minmaxtemp = minpdj_maxdPSI( tbtemp_t, colnn = c("genesymbol",paste0("minpadj_maxPdPSI_",coi) , paste0("dPSI_",coi) ) )  ##  == == == == == == == == == == == == == == == == == == Call the function
    
    head(minmaxtemp)
    
    # ## optional: clean up colnames ##
    # colnames(minmaxtemp) = gsub("filteredevents_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) ) #   gsub("filteredevents_TST11742_","", sapply( strsplit( colnames(minmaxtemp) , paste0("_dPSI") ) , "[[", 1) )
    # 
    # ## merge ##
    # if(i ==  1){ minmaxall = minmaxtemp
    # } else{  minmaxall = merge(minmaxall, minmaxtemp, by.x = "genesymbol", by.y = "genesymbol" , all = T) }
    
    minmaxtemp = tbtemp_t
    head(minmaxtemp)
    
}

head(minmaxall) 
dim(minmaxall) # [1] [1] 2331   85 for GM GM09677c SMA fibroblasts  cells  *** 2195   61 for SySy5Y cells                                  [2] 162  60 GM cells [1] 249   84 GM cells 

colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_maxP_", colnames(minmaxall)) # colnames(minmaxall) = gsub("filteredevents_","", colnames(minmaxall)) #              gsub("filteredevents_TST11955_","", colnames(minmaxall))
#colnames(minmaxall) = gsub("WT100","WT_100", colnames(minmaxall))  #colnames(minmaxall) = gsub("WT500","WT_500", colnames(minmaxall)) 

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

dim(minmaxall) #  2331 82 GM cells
head(minmaxall) 
########################################### 08-01-2023

minmaxall_dPSI  = minmaxall[, grep( "genesymbol|^dPSI", colnames(minmaxall) ) ] #drop padj, only keep dPSI, eq. delete 50% columns. <<============ colnames(minmaxall) = gsub("minpadj_maxPdPSI_","Padj_mP_", colnames(minmaxall) )  #####

head(minmaxall_dPSI )
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

# 1) all samples 
y = colSums(dPSI_bi)
names(y) = gsub( "_3x" , "_03x", names(y) )  # y %>% select(sub('_ln$', '', filter_vector))
y = y[order(names(y), decreasing=T)]

fout=paste0("binary_numDSG_",param,"vertica.png")
colorr=c( rep(rgb(0.7,0.7,0.7),4), rep(rgb(0.3,0.3,0.3),4) )    #"violet" , "grey"  ) # blue, navy, grey "indianred" , "violet"   "cyan" , "violet"

y

# 2) samples to DMSO
y = colSums(dPSI_bi)
names(y) = gsub( "_3x" , "_03x", names(y) )  # y %>% select(sub('_ln$', '', filter_vector))
y = y[order(names(y), decreasing=T)]

y = y[ grepl( "vs_.._DMSO" , names(y) ) ] # y %>% select(sub('_ln$', '', filter_vector)) # y = y[ grepl( "vs_SH_DMSO" , names(y) ) ] # y %>% select(sub('_ln$', '', filter_vector))

fout=paste0("binary_numDSG_",param,"vertical_DMSO.png")
colorr=c( "cyan" , "grey"  ) # blue, navy, grey "indianred" , "violet"   "cyan" , "violet"  colorr=c( rep(rgb(0.7,0.7,0.7),2),rep(rgb(0.3,0.3,0.3),2)   )

# 2) samples to Risdiplam  1949634
y = colSums(dPSI_bi)
names(y) = gsub( "_3x" , "_03x", names(y) )  # y %>% select(sub('_ln$', '', filter_vector))
y = y[order(names(y), decreasing=T)]
y = y[ grepl( "vs_.._BIO_1949634" , names(y) ) ] # y %>% select(sub('_ln$', '', filter_vector))
fout=paste0("binary_numDSG_",param,"vertical_Risdiplam.png")
colorr=c( "violet" , "grey"  ) # blue, navy, grey "indianred" , "violet"   "cyan" , "violet"

#### plotting
png(fout,width=800, height=length(y)*35) # png(fout,width=800, height=length(y)*25)
source(pformat)
par(mar=c(10,28,1,1))
bp=barplot(y,  xlim = c(0, max(y)*1.2 ), border=NA,  names.arg=gsub("dPSI","", names(y) ),
           las=2,  cex.axis=0.8,  cex.names=1.0, #cex.axis	 #           expansion factor for numeric axis labels.            las=2,  cex.axis=1.2,  cex.names=1.5 font size
           horiz=T, col=colorr,
           sub = "number(DSG)", font.sub = 4,  cex.sub = 2.0)  #http://howtoinr.weebly.com/customize-labels1.html col.sub = "Red",
text(y, bp, labels = y, pos=4, cex=1.2, offset = 0.3) # lable the differenrtial gene count on the top of BAR
dev.off()


fout = paste0("UpSet_0_top_5", ".png")
png(fout,height = 1000, width = length(y)*75) #png(fout,height = 500, width = length(y)*50)
upset(dPSI_bi ,  mainbar.y.label = "Differential Splicing Differential Splicing Gene Overlaps", point.size = 3.5, text.scale = c(2.5, 2, 1.5, 1.5, 2, 2))    # https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
dev.off()

ds <- read.csv("Crop_recommendation.csv", header = TRUE)

ggplot(ds, aes(x=label, y=temperature)) + geom_boxplot() 



# create a new dataframe crop_means
crop_means <- ds %>% 
    group_by(label) %>% 
    summarize(mean_temperature=mean(temperature)) 
crop_means

# Creating barplots of means
ggplot(crop_means, aes(x=label, y=mean_temperature)) +
    geom_bar(stat="identity") 




# https://epirhandbook.com/en/ggplot-basics.html
pacman::p_load(
  tidyverse,      # includes ggplot2 and other data management tools
  rio,            # import/export
  here,           # file locator
  stringr         # working with characters   
)

setwd("/edgehpc/dept/compbio/users/zgao1/Coding_system_tool_utilities/tool_utilities/plotting")
# wget https://github.com/appliedepi/epirhandbook_eng/raw/master/data/case_linelists/linelist_cleaned.rds

linelist <- rio::import("linelist_cleaned.rds")
linelist


linelist <- linelist %>%
  mutate(
    gender_disp = case_when(gender == "m" ~ "Male",        # m to Male 
                            gender == "f" ~ "Female",      # f to Female,
                            is.na(gender) ~ "Unknown"),    # NA to Unknown
    
    outcome_disp = replace_na(outcome, "Unknown")          # replace NA outcome with "unknown"
  )

symptoms_data <- linelist %>% 
  select(c(case_id, fever, chills, cough, aches, vomit))

symptoms_data_long <- symptoms_data %>%    # begin with "mini" linelist called symptoms_data
  
  pivot_longer(
    cols = -case_id,                       # pivot all columns except case_id (all the symptoms columns)
    names_to = "symptom_name",             # assign name for new column that holds the symptoms
    values_to = "symptom_is_present") %>%  # assign name for new column that holds the values (yes/no)
  
  mutate(symptom_is_present = replace_na(symptom_is_present, "unknown")) # convert NA to "unknown"

# # plot data from my_data columns as red points
# ggplot(data = my_data)+                   # use the dataset "my_data"
#   geom_point(                             # add a layer of points (dots)
#     mapping = aes(x = col1, y = col2),    # "map" data column to axes
#     color = "red")+                       # other specification for the geom
#   labs()+                                 # here you add titles, axes labels, etc.
#   theme()                                 # here you adjust color, font, size etc of non-data plot elements (axes, title, etc.) 


age_by_wt <- ggplot(
  data = linelist,   # set data
  mapping = aes(     # map aesthetics to column values
    x = age,           # map x-axis to age            
    y = wt_kg,         # map y-axis to weight
    color = age))+     # map color to age
  geom_point()+           # display data as points
  labs(
    title = "Age and weight distribution",
    subtitle = "Fictional Ebola outbreak, 2014",
    x = "Age in years",
    y = "Weight in kilos",
    color = "Age",
    caption = stringr::str_glue("Data as of {max(linelist$date_hospitalisation, na.rm=T)}"))

age_by_wt

class(linelist$hospital)

library("data.table")
library(janitor)
linelist %>% 
  tabyl(hospital)

dev.off()
dim(linelist)
ggplot(data = linelist) # a blank canvas !!!!!!!!!!!!!!

# A) Outcomes in all cases
ggplot(linelist %>% drop_na(outcome)) + 
  geom_bar(aes(y = fct_rev(hospital)), width = 0.7) +
  theme_minimal()+
  labs(title = "A) Number of cases by hospital",
       y = "Hospital")


# B) Outcomes in all cases by hosptial
ggplot(linelist %>% drop_na(outcome)) + 
  geom_bar(aes(y = fct_rev(hospital), fill = outcome), width = 0.7) +
  theme_minimal()+
  theme(legend.position = "bottom") +
  labs(title = "B) Number of recovered and dead Ebola cases, by hospital",
       y = "Hospital")


outcomes2 <- linelist %>% 
  drop_na(outcome) %>% 
  count(hospital, outcome) %>%  # get counts by hospital and outcome
  group_by(hospital) %>%        # Group so proportions are out of hospital total
  mutate(proportion = n/sum(n)*100) # calculate proportions of hospital total

head(outcomes2) # Preview data


ggplot(outcomes2) +  
  geom_col(
    mapping = aes(
      x = proportion,                 # show pre-calculated proportion values
      y = fct_rev(hospital),          # reverse level order so missing/other at bottom
      fill = outcome),                # stacked by outcome
    width = 0.5)+                    # thinner bars (out of 1)
  theme_minimal() +                  # Minimal theme 
  theme(legend.position = "bottom")+
  labs(subtitle = "Number of recovered and dead Ebola cases, by hospital",
       fill = "Outcome",             # legend title
       y = "Count",                  # y axis title
       x = "Hospital of admission")+ # x axis title
  scale_fill_manual(                 # adding colors manually
    values = c("Death"= "#3B1c8C",
               "Recover" = "#21908D" )) 




R.Version()
#rm(list = ls()) #install.packages("ggVennDiagram") # 15min
#library("ggVennDiagram")
library(ggplot2)
require("UpSetR") # movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), header = T, sep = ";")
library(dplyr)

getwd()

data1 <- 
    diamonds %>% 
    group_by(clarity) %>% 
    summarize(m = mean(price))  

## creating dataset #2
data2 <- 
    diamonds %>% 
    group_by(clarity, cut) %>% 
    summarize(m = mean(price))

## graphing data points from 2 different datasets on one graph
ggplot() +
    geom_point(data = data1, aes(x = clarity, y = m), color = "blue") + # must include argument label "data"
    geom_point(data = data2, aes(x = clarity, y = m))

barplot1=c(10,2,5,4,6,5,8,10,5,9)
barplot2=c(9,5,6,4,7,1,2,6,2,6)
barplot3=c(4,2,9,4,3,5,7,10,10,3)
data <- data.frame(barplot1,barplot2,barplot3)

# plotting multiple bar plots
barplot(as.matrix(data),
        main="Multiple Bar Plots",
        
        # setting y label only
        # because x-label will be our
        # barplots name
        ylab="Count",
        
        # to plot the bars vertically
        beside=TRUE,
        #beside=F,
)

# https://ssc.wisc.edu/sscc/pubs/dvr/index.html#download-the-data
# only works in R4
library(ggplot2)
library(dplyr)
library(scales)

#acs <- readRDS("acs.rds")
acs <- readRDS(url("https://sscc.wisc.edu/sscc/pubs/dvr/acs.rds"))
acs

acs |> 
  group_by(race, edu) |> 
  summarize(n = n()) |> 
  mutate(perc = 100*n/sum(n)) |> 
  ggplot(aes(x = race, y = perc)) +
  geom_bar(stat = "identity") +
  facet_grid(~ edu) + 
  coord_flip()


# https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r 
library(tidyverse)

dataframe1 = read.table(text="sl zone   meangpp
                        1     1 4.050161
                        2     2 7.729265
                        3     3 3.408220
                        4     4 4.884040
                        5     5 4.258422
                        6     6 2.906374
                        7     7 2.241984
                        8     8 4.703197
                        9     9 3.617657
                        10   10 2.712997
                        11   12 3.589406", header=T)
dataframe2 = read.table(text="sl zone   meangpp
                        1     1 5.4153407
                        2     2 4.2429236
                        3     3 4.5719178
                        4     4 3.1215946
                        5     5 4.9222054
                        6     6 3.0384872
                        7     7 1.9293729
                        8     8 8.9709741
                        9     9 7.8904906
                        10   10 6.6410986
                        11   12 5.5011823", header=T)

df <- bind_rows("dataframe1" = dataframe1, "dataframe2" = dataframe2, .id = "groups")

df %>% 
  ggplot(aes(x=factor(zone), y=meangpp, fill = groups)) + 
  geom_col(position = position_dodge())


dataframe3 = read.table(text="sl zone   meannpp
                        1     1 5.4153407
                        2     2 4.2429236
                        3     3 4.5719178
                        4     4 3.1215946
                        5     5 4.9222054
                        6     6 3.0384872
                        7     7 1.9293729
                        8     8 8.9709741
                        9     9 7.8904906
                        10   10 6.6410986
                        11   12 5.5011823", header=T)
dataframe4 = read.table(text="sl zone   meannpp
                        1     1 4.050161
                        2     2 7.729265
                        3     3 3.408220
                        4     4 4.884040
                        5     5 4.258422
                        6     6 2.906374
                        7     7 2.241984
                        8     8 4.703197
                        9     9 3.617657
                        10   10 2.712997
                        11   12 3.589406", header=T)

df <- bind_rows("dataframe3" = dataframe3, "dataframe4" = dataframe4, .id = "groups")

df %>% 
  ggplot(aes(x=factor(zone), y=meannpp, fill = groups)) + 
  geom_col(position = position_dodge())



library(patchwork)

df <- bind_rows(dataframe1, dataframe2, .id = "id")

axis_labels <- c("first", "second", "third", "D", "E", "F", "G", "H", "I", "J", "K", "L")
axis_labels <- setNames(axis_labels, 1:12)

p1 <- ggplot(df, aes(x=factor(zone), y=meangpp, fill = id)) + 
  geom_col(position = position_dodge()) +
  scale_fill_discrete(labels = c("1" = "dataframe 1 & 3", "2" = "dataframe 2 & 4")) +
  scale_x_discrete(labels = axis_labels) +
  theme(axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(plot.margin = unit(rep(0, 4), "pt"))  

df <- bind_rows(dataframe3, dataframe4, .id = "id")

p2 <- ggplot(df, aes(x=factor(zone), y=meannpp, fill = id)) + 
  geom_col(position = position_dodge()) +
  scale_fill_discrete(labels = c("1" = "dataframe 1 & 3", "2" = "dataframe 2 & 4")) +
  scale_x_discrete(labels = axis_labels) +
  theme(plot.margin = unit(rep(0, 4), "pt"))

p1 / p2 +
  plot_layout(guides = 'collect')


# https://stackoverflow.com/questions/6644997/showing-data-values-on-stacked-bar-chart-in-ggplot2
Year      <- c(rep(c("2006-07", "2007-08", "2008-09", "2009-10"), each = 4))
Category  <- c(rep(c("A", "B", "C", "D"), times = 4))
Frequency <- c(168, 259, 226, 340, 216, 431, 319, 368, 423, 645, 234, 685, 166, 467, 274, 251)
Data      <- data.frame(Year, Category, Frequency)
library(ggplot2)

ggplot(Data, aes(x = Year, y = Frequency, fill = Category, label = Frequency)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))



Year      <- c("2006-07", "2007-08", "2008-09", "2009-10")
Frequency <- c(168, 259, 226, 340)
Data      <- data.frame(Year,  Frequency)
library(ggplot2)
# https://ggplot2.tidyverse.org/reference/geom_text.html
# https://ggplot2.tidyverse.org/reference/position_dodge.html
position_dodge(width = NULL, preserve = "total")
position_dodge2(
  width = NULL,
  preserve = "total",
  padding = 0.1,
  reverse = FALSE
)

Data1 = Data[order(Data$Frequency), ]
Data1
dev.off()
#annotate("text", label = "plot mpg vs. wt", x = 2, y = 15, size = 8, colour = "red")
ggplot(Data, aes(x = Year, y = Frequency, label = Frequency)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, vjust = -0.5) 

ggplot(Data1, aes(x = Year, y = Frequency, label = Frequency)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, vjust = -0.5) 





#####   label x-axie according to data frame , not default alphabeticall order
# https://stackoverflow.com/questions/25463029/ggplot2-how-to-label-the-x-ticks-using-a-column-in-my-data-frame
df <- data.frame(team=c('Mavs', 'Heat', 'Nets', 'Lakers'),
points=c(100, 122, 104, 109),
order=c(4, 1, 3, 2) )

#view data frame
df
# team points
# 1   Mavs    100
# 2   Heat    122
# 3   Nets    104
# 4 Lakers    109

str(df)

library(ggplot2)ls

#create bar plot
ggplot(df, aes(x=team, y=points)) +
  geom_col()


df$team <- reorder(df$team)
df$team <- factor(df$team)
df$team <- reorder(df$team, df$rank)

df

ggplot(data = df, aes(x = reorder(team, order), y = points)) +
  geom_bar(stat = "identity", fill="grey")

##### https://ggplot2.tidyverse.org/reference/geom_text.html
# add all kinds of text lable in different locations


# run.04.DSG_counts_plot_ranking_batch123_DSG_2023_05_15_Manugs_request.R

R.Version()
rm(list = ls())
#library(gplots)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(withr); # in the Packages withr, select the 2.4.2 version  shown as library(withr, lib.loc = "/opt/R/4.0.3/lib/R/library")
#install.packages("patchwork")
library(patchwork) #https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r
library(dplyr) #outer_join  #library(plyr) #join # 
library(forcats)  # fct_relevel

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors

#dout = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/" # din = "./res.02.DSG_counts/" 
dout = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/" # din = "./res.02.DSG_counts/" 
dir.create(dout, recursive = T)
setwd(dout); 
getwd()                 

################################################## ======== prepair the plotting ========== ##################################################
library(stringr)
library(VennDiagram)
library(dplyr)

################################################## ======== ========== ##################################################
# /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/run.02-1.DSG_counts_plot_6152_1.R
prepair_venn_plot_data = function( compounds = comp_1_compounds ) {
  for(i in 1:length(compounds) ) {
    # i=1
    # i=2
    # i=3
    # i=4
    comp    = compounds[i]     #concent = concentration[i]
    
    dPSI_temp = dPSI_bi[ , grep(comp, colnames(dPSI_bi)),  drop=FALSE ] # avoid the sinlge columne to vector
    #  dPSI_temp = dPSI_temp[ , grep(concent , colnames(  dPSI_temp)) ]
    dim(dPSI_temp)
    head(dPSI_temp)
    ## merge ##
    if( i==1 ){ 
      dPSI_selected = dPSI_temp
    }else if( i == 2 ){  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by = 0, all=T) 
    } else {  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by.x =  "Row.names", by.y = 0, all=T) 
    }
    
  }
  
  row.names(dPSI_selected) = dPSI_selected$Row.names  
  dPSI_selected$Row.names  = NULL
  
  dPSI_selected = as.matrix(dPSI_selected) #  dPSI_selected = as.numeric(dPSI_selected)
  
  head(    dPSI_selected)
  dPSI_temp = dPSI_selected[rowSums(dPSI_temp)>0,]
  dim(dPSI_temp)
  head(dPSI_temp)
  
  dPSI_temp = dPSI_temp[order(rowSums(dPSI_temp),  decreasing=T) ,]
  dPSI_temp
}





# ##############################   risdiplam in 3 iPSC_MN cells ===2019     RG7961 = risdplam
din0 = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/"
fin0 = "iPSC_MN_binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
dPSI_bi0= read.delim2( paste0(din0, fin0)  , stringsAsFactors = F, sep = "," ) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dPSI_bi0
row.names(dPSI_bi0) = dPSI_bi0$X
dPSI_bi0$X = NULL

## total number of DSG ###
y_0th        = colSums(dPSI_bi0) ##################### get the totall DS events
y_0th        = y_0th[order(names(y_0th) )]
#names(y_0th) = gsub("dPSI", "", names(y_0th))
y_0th # class(y_0th) #y_0th_rank = rank(y_0th ) # y_rank = rank(-y )


####### iPSC-MN Risdiplam 3x in iPSC-MN
##############################=========      2023-05-17 ===========================
comp_3_compounds = c("Risdiplam_3x_in_GM24468D",  "Risdiplam_3x_in_HNDS0030_01" ,  "Risdiplam_3x_in_HNDS006_01C2" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi0)

GM24468D.RG7961_3x_vector     =   rownames( dPSI_bi0[  (dPSI_bi0$GM24468D.RG7961_3x     == 1), ] ) #  RG7961_3x_vector  = Risdplam 120nM
GM24468D.RG7961_3x_vector
HNDS0030_01.RG7961_3x_vector  =   rownames( dPSI_bi0[  (dPSI_bi0$HNDS0030_01.RG7961_3x  == 1), ] ) # 
HNDS0030_01.RG7961_3x_vector
HNDS006_01C2.RG7961_3x_vector =   rownames( dPSI_bi0[  (dPSI_bi0$HNDS006_01C2.RG7961_3x == 1), ] ) # 
HNDS006_01C2.RG7961_3x_vector

getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(GM24468D.RG7961_3x_vector, HNDS0030_01.RG7961_3x_vector, HNDS006_01C2.RG7961_3x_vector   ),
  #category.names = c("Risdiplam_3x_in_GM24468D", "Risdiplam_3x_in_HNDS0030_01", "Risdiplam_3x_in_HNDS006_01C2"  ),
  category.names = c("", "", ""  ),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_risdiplam_3x_in_GM24468D_HNDS0030_01_HNDS006_01C2.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_risdiplam_3x_in_GM24468D_HNDS0030_01_HNDS006_01C2_no_name_NGN2.png'),      ###### output in the working directory!!!
  output=TRUE
)

dev.off()


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

# colnames(dPSI_bi_tmp) =  c("2184088_NGN2_3x", "2184090_NGN2_3x", "2174714_NGN2_3x", "2184741_NGN2_3x", "2174748_NGN2_3x", "2175420_NGN2_3x", 
#                            "2006152_NGN2_3x", "2186527_NGN2_3x", "2176866_NGN2_3x", "2186960_NGN2_3x", "2139701_NGN2_3x" ) 


head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)

dPSI_bi_batch2_iPSC = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25
#dPSI_bi_batch2_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25


dPSI_bi = dPSI_bi_tmp
y_2nd        = colSums(dPSI_bi) ##################### get the totall DS events
y_2nd        = y_2nd[order(names(y_2nd) )]
#names(y_2nd) = gsub("dPSI", "", names(y_2nd))
y_2nd # class(y_2nd) #y_2nd_rank = rank(y_2nd ) # y_rank = rank(-y ) #y_2nd_rank


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
dPSI_bi_batch3_iPSC = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425


colnames(dPSI_bi_tmp) =  c(  "1755497_NGN2_3x", "2006152_NGN2_3x", "2189972_NGN2_3x", "2194943_NGN2_3x", "2195127_NGN2_3x", "2195327_NGN2_3x", "2197294_NGN2_3x" ) 
dPSI_bi_batch3 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425
dPSI_bi_batch3_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425
dPSI_bi_batch3_NGN2


dPSI_bi = dPSI_bi_tmp
y_3rd        = colSums(dPSI_bi) ##################### get the totall DS events
y_3rd        = y_3rd[order(names(y_3rd) , decreasing = T)]
#names(y_3rd) = gsub("dPSI", "", names(y_3rd))
y_3rd #class(y_3rd)  #y_3rd_rank = rank(y_3rd ) # y_rank = rank(-y )

#### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels

################################################## ======== ========== ##################################################
# /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/run.02-1.DSG_counts_plot_6152_1.R
prepair_venn_plot_data = function( compounds = comp_1_compounds ) {
  for(i in 1:length(compounds) ) {
    # i=1
    # i=2
    # i=3
    # i=4
    comp    = compounds[i]     #concent = concentration[i]
    
    dPSI_temp = dPSI_bi[ , grep(comp, colnames(dPSI_bi)),  drop=FALSE ] # avoid the sinlge columne to vector
    #  dPSI_temp = dPSI_temp[ , grep(concent , colnames(  dPSI_temp)) ]
    dim(dPSI_temp)
    head(dPSI_temp)
    ## merge ##
    if( i==1 ){ 
      dPSI_selected = dPSI_temp
    }else if( i == 2 ){  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by = 0, all=T) 
    } else {  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by.x =  "Row.names", by.y = 0, all=T) 
    }
    
  }
  
  row.names(dPSI_selected) = dPSI_selected$Row.names  
  dPSI_selected$Row.names  = NULL
  
  dPSI_selected = as.matrix(dPSI_selected) #  dPSI_selected = as.numeric(dPSI_selected)
  
  head(    dPSI_selected)
  dPSI_temp = dPSI_selected[rowSums(dPSI_temp)>0,]
  dim(dPSI_temp)
  head(dPSI_temp)
  
  dPSI_temp = dPSI_temp[order(rowSums(dPSI_temp),  decreasing=T) ,]
  dPSI_temp
}


##############################=========      2023-04-25 ===========================
# =============================================== PTC518 comp_3_title     = "PTC-518-3x: iPSC vs NGN2 vs Sy5Y"
# =============================================== # 1755497 is Branaplam, 1949634 Risdiplam, 2186960
#######NGN2
comp_3_compounds = c("iPSC_2186960_3x",  "iPSC_1755497_3x" ,  "iPSC_2197294_3x" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi_batch2)
head(dPSI_bi_batch3)

iPSC_2186960_3x_vector =   rownames( dPSI_bi_batch2[  (dPSI_bi_batch2$`2186960_iPSC_3x` == 1), ] ) # 6960 from batch2
iPSC_2186960_3x_vector

iPSC_1755497_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`1755497_iPSC_3x` == 1), ] ) # Branaplam
iPSC_1755497_3x_vector

iPSC_2197294_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`2197294_iPSC_3x` == 1), ] ) # PTC518
iPSC_2197294_3x_vector

library(stringr)
library(VennDiagram)
library(dplyr)


getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(iPSC_2186960_3x_vector, iPSC_1755497_3x_vector, iPSC_2197294_3x_vector  ),
  #category.names = c("2186960", "Branaplam", "PTC518"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_no_name.png'),      ###### output in the working directory!!!
  output=TRUE
)
dev.off()


#######NGN2
##############################=========      2023-04-25 ===========================
# =============================================== PTC518 comp_3_title     = "PTC-518-3x: NGN2 vs NGN2 vs Sy5Y"
# =============================================== # 1755497 is Branaplam, 1949634 Risdiplam, 2186960
comp_3_compounds = c("NGN2_2186960_3x",  "NGN2_1755497_3x" ,  "NGN2_2197294_3x" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi_batch2_NGN2)
head(dPSI_bi_batch3_NGN2)

NGN2_2186960_3x_vector =   rownames( dPSI_bi_batch2[  (dPSI_bi_batch2$`2186960_NGN2_3x` == 1), ] ) # 6960 from batch2
NGN2_2186960_3x_vector

NGN2_1755497_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`1755497_NGN2_3x` == 1), ] ) # Branaplam
NGN2_1755497_3x_vector

NGN2_2197294_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`2197294_NGN2_3x` == 1), ] ) # PTC518
NGN2_2197294_3x_vector

getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(NGN2_2186960_3x_vector, NGN2_1755497_3x_vector, NGN2_2197294_3x_vector  ),
  #category.names = c("2186960", "Branaplam", "PTC518"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_NGN2.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_no_name_NGN2.png'),      ###### output in the working directory!!!
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



# geatmap
# https://stackoverflow.com/questions/54342274/pheatmap-display-numbers-argument
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
We can render the plot and save the result with

xx <- pheatmap(test)
#Then you can output to this to a file by opening a graphics device and re-drawing the result the way it's done in the main function

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
save_pheatmap_pdf(xx, "test.pdf")



# Venn
R.Version()
rm(list = ls())
#library(gplots)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(withr); # in the Packages withr, select the 2.4.2 version  shown as library(withr, lib.loc = "/opt/R/4.0.3/lib/R/library")
#install.packages("patchwork")
library(patchwork) #https://stackoverflow.com/questions/67858336/how-to-plot-two-grouped-barplots-vertically-with-single-x-axis-in-r
library(dplyr) #outer_join  #library(plyr) #join # 
library(forcats)  # fct_relevel

pformat = "/edgehpc/dept/compbio/users/dhuh/software/R/R_modified/pformat_whitebG.r"
jet     = colorRampPalette(c("blue","green","yellow","orange","darkred")) # make gradient of colors

#dout = "/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/" # din = "./res.02.DSG_counts/" 
dout = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/" # din = "./res.02.DSG_counts/" 
dir.create(dout, recursive = T)
setwd(dout); 
getwd()                  # [1] "/mnt/depts/dept04/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig"



################################################## ======== ========== ##################################################
library(stringr)
library(VennDiagram)
library(dplyr)

################################################## ======== ========== ##################################################
# /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/run.02-1.DSG_counts_plot_6152_1.R
prepair_venn_plot_data = function( compounds = comp_1_compounds ) {
  for(i in 1:length(compounds) ) {
    # i=1
    # i=2
    # i=3
    # i=4
    comp    = compounds[i]     #concent = concentration[i]
    
    dPSI_temp = dPSI_bi[ , grep(comp, colnames(dPSI_bi)),  drop=FALSE ] # avoid the sinlge columne to vector
    #  dPSI_temp = dPSI_temp[ , grep(concent , colnames(  dPSI_temp)) ]
    dim(dPSI_temp)
    head(dPSI_temp)
    ## merge ##
    if( i==1 ){ 
      dPSI_selected = dPSI_temp
    }else if( i == 2 ){  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by = 0, all=T) 
    } else {  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by.x =  "Row.names", by.y = 0, all=T) 
    }
    
  }
  
  row.names(dPSI_selected) = dPSI_selected$Row.names  
  dPSI_selected$Row.names  = NULL
  
  dPSI_selected = as.matrix(dPSI_selected) #  dPSI_selected = as.numeric(dPSI_selected)
  
  head(    dPSI_selected)
  dPSI_temp = dPSI_selected[rowSums(dPSI_temp)>0,]
  dim(dPSI_temp)
  head(dPSI_temp)
  
  dPSI_temp = dPSI_temp[order(rowSums(dPSI_temp),  decreasing=T) ,]
  dPSI_temp
}















# ##############################   000000----------------------------------------- iPSC Note ======PTC518 compound is BIO-2197294
# 1) venn diagram for these 3 compounds. What is the percent overlap between compounds in human cells?
# 2) excel list of OTs occurring in each compound. Which OTs occur in all 3 compounds?
# BIO-2197294 (PTC518)
# BIO-1949634 (risdiplam)
# BIO-1755497 (branaplam)


din0 = "/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DownStream_Results_batch1_2_3/analysis.03.splice_analysis_Mangus_0511_request/"
fin0 = "iPSC_MN_binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv" # fin = "binary_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M502022-04-14.csv" 
dPSI_bi0= read.delim2( paste0(din0, fin0)  , stringsAsFactors = F, sep = "," ) # dPSI_bi = read.csv( paste0(fin1) , stringsAsFactors = F , row.names = 1 , check.names = F)

dPSI_bi0
row.names(dPSI_bi0) = dPSI_bi0$X
dPSI_bi0$X = NULL

## total number of DSG ###
y_0th        = colSums(dPSI_bi0) ##################### get the totall DS events
y_0th        = y_0th[order(names(y_0th) )]
#names(y_0th) = gsub("dPSI", "", names(y_0th))
y_0th # class(y_0th) #y_0th_rank = rank(y_0th ) # y_rank = rank(-y )


####### iPSC-MN Risdiplam 3x in iPSC-MN
##############################=========      2023-05-17 ===========================
comp_3_compounds = c("Risdiplam_3x_in_GM24468D",  "Risdiplam_3x_in_HNDS0030_01" ,  "Risdiplam_3x_in_HNDS006_01C2" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi0)

GM24468D.RG7961_3x_vector =   rownames( dPSI_bi0[  (dPSI_bi0$GM24468D.RG7961_3x == 1), ] ) #  RG7961_3x_vector  = Risdplam 120nM
GM24468D.RG7961_3x_vector

HNDS0030_01.RG7961_3x_vector =   rownames( dPSI_bi0[  (dPSI_bi0$HNDS0030_01.RG7961_3x == 1), ] ) # 
HNDS0030_01.RG7961_3x_vector

HNDS006_01C2.RG7961_3x_vector =   rownames( dPSI_bi0[  (dPSI_bi0$HNDS006_01C2.RG7961_3x == 1), ] ) # 
HNDS006_01C2.RG7961_3x_vector

getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(GM24468D.RG7961_3x_vector, HNDS0030_01.RG7961_3x_vector, HNDS006_01C2.RG7961_3x_vector   ),
  category.names = c("Risdiplam_3x_in_GM24468D", "Risdiplam_3x_in_HNDS0030_01", "Risdiplam_3x_in_HNDS006_01C2"  ),
  #category.names = c("", "", ""  ),
  imagetype = "png",
  filename = paste0(out_dir, '/venn_diagramm_risdiplam_3x_in_GM24468D_HNDS0030_01_HNDS006_01C2.png'),      ###### output in the working directory!!!
  #filename = paste0(out_dir, '/venn_diagramm_risdiplam_3x_in_GM24468D_HNDS0030_01_HNDS006_01C2_no_name_NGN2.png'),      ###### output in the working directory!!!
  output=TRUE
)

dev.off()


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

# colnames(dPSI_bi_tmp) =  c("2184088_NGN2_3x", "2184090_NGN2_3x", "2174714_NGN2_3x", "2184741_NGN2_3x", "2174748_NGN2_3x", "2175420_NGN2_3x", 
#                            "2006152_NGN2_3x", "2186527_NGN2_3x", "2176866_NGN2_3x", "2186960_NGN2_3x", "2139701_NGN2_3x" ) 


head(dPSI_bi_tmp)
dim(dPSI_bi_tmp)

dPSI_bi_batch2_iPSC = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25
#dPSI_bi_batch2_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ====== 04-25


dPSI_bi = dPSI_bi_tmp
y_2nd        = colSums(dPSI_bi) ##################### get the totall DS events
y_2nd        = y_2nd[order(names(y_2nd) )]
#names(y_2nd) = gsub("dPSI", "", names(y_2nd))
y_2nd # class(y_2nd) #y_2nd_rank = rank(y_2nd ) # y_rank = rank(-y ) #y_2nd_rank


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
dPSI_bi_batch3_iPSC = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425


colnames(dPSI_bi_tmp) =  c(  "1755497_NGN2_3x", "2006152_NGN2_3x", "2189972_NGN2_3x", "2194943_NGN2_3x", "2195127_NGN2_3x", "2195327_NGN2_3x", "2197294_NGN2_3x" ) 
dPSI_bi_batch3 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425
dPSI_bi_batch3_NGN2 = dPSI_bi_tmp #### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels ======0425
dPSI_bi_batch3_NGN2


dPSI_bi = dPSI_bi_tmp
y_3rd        = colSums(dPSI_bi) ##################### get the totall DS events
y_3rd        = y_3rd[order(names(y_3rd) , decreasing = T)]
#names(y_3rd) = gsub("dPSI", "", names(y_3rd))
y_3rd #class(y_3rd)  #y_3rd_rank = rank(y_3rd ) # y_rank = rank(-y )

#### ==========+++++++++++========== 1-1  ==> 1 plot with 3 panels

################################################## ======== ========== ##################################################
# /camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/run.02-1.DSG_counts_plot_6152_1.R
prepair_venn_plot_data = function( compounds = comp_1_compounds ) {
  for(i in 1:length(compounds) ) {
    # i=1
    # i=2
    # i=3
    # i=4
    comp    = compounds[i]     #concent = concentration[i]
    
    dPSI_temp = dPSI_bi[ , grep(comp, colnames(dPSI_bi)),  drop=FALSE ] # avoid the sinlge columne to vector
    #  dPSI_temp = dPSI_temp[ , grep(concent , colnames(  dPSI_temp)) ]
    dim(dPSI_temp)
    head(dPSI_temp)
    ## merge ##
    if( i==1 ){ 
      dPSI_selected = dPSI_temp
    }else if( i == 2 ){  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by = 0, all=T) 
    } else {  
      dPSI_selected = merge(dPSI_selected, dPSI_temp , by.x =  "Row.names", by.y = 0, all=T) 
    }
    
  }
  
  row.names(dPSI_selected) = dPSI_selected$Row.names  
  dPSI_selected$Row.names  = NULL
  
  dPSI_selected = as.matrix(dPSI_selected) #  dPSI_selected = as.numeric(dPSI_selected)
  
  head(    dPSI_selected)
  dPSI_temp = dPSI_selected[rowSums(dPSI_temp)>0,]
  dim(dPSI_temp)
  head(dPSI_temp)
  
  dPSI_temp = dPSI_temp[order(rowSums(dPSI_temp),  decreasing=T) ,]
  dPSI_temp
}


##############################=========      2023-04-25 ===========================
# =============================================== PTC518 comp_3_title     = "PTC-518-3x: iPSC vs NGN2 vs Sy5Y"
# =============================================== # 1755497 is Branaplam, 1949634 Risdiplam, 2186960
#######NGN2
comp_3_compounds = c("iPSC_2186960_3x",  "iPSC_1755497_3x" ,  "iPSC_2197294_3x" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi_batch2)
head(dPSI_bi_batch3)

iPSC_2186960_3x_vector =   rownames( dPSI_bi_batch2[  (dPSI_bi_batch2$`2186960_iPSC_3x` == 1), ] ) # 6960 from batch2
iPSC_2186960_3x_vector

iPSC_1755497_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`1755497_iPSC_3x` == 1), ] ) # Branaplam
iPSC_1755497_3x_vector

iPSC_2197294_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`2197294_iPSC_3x` == 1), ] ) # PTC518
iPSC_2197294_3x_vector

library(stringr)
library(VennDiagram)
library(dplyr)


getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(iPSC_2186960_3x_vector, iPSC_1755497_3x_vector, iPSC_2197294_3x_vector  ),
  #category.names = c("2186960", "Branaplam", "PTC518"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_no_name.png'),      ###### output in the working directory!!!
  output=TRUE
)
dev.off()


#######NGN2
##############################=========      2023-04-25 ===========================
# =============================================== PTC518 comp_3_title     = "PTC-518-3x: NGN2 vs NGN2 vs Sy5Y"
# =============================================== # 1755497 is Branaplam, 1949634 Risdiplam, 2186960
comp_3_compounds = c("NGN2_2186960_3x",  "NGN2_1755497_3x" ,  "NGN2_2197294_3x" ) # latter 2 are Branaplam and PTC518 

head(dPSI_bi_batch2_NGN2)
head(dPSI_bi_batch3_NGN2)

NGN2_2186960_3x_vector =   rownames( dPSI_bi_batch2[  (dPSI_bi_batch2$`2186960_NGN2_3x` == 1), ] ) # 6960 from batch2
NGN2_2186960_3x_vector

NGN2_1755497_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`1755497_NGN2_3x` == 1), ] ) # Branaplam
NGN2_1755497_3x_vector

NGN2_2197294_3x_vector =   rownames( dPSI_bi_batch3[  (dPSI_bi_batch3$`2197294_NGN2_3x` == 1), ] ) # PTC518
NGN2_2197294_3x_vector

getwd()
out_dir=getwd()
dev.off()

venn.diagram(
  x = list(NGN2_2186960_3x_vector, NGN2_1755497_3x_vector, NGN2_2197294_3x_vector  ),
  #category.names = c("2186960", "Branaplam", "PTC518"  ),
  category.names = c("", "", ""  ),
  #category.names = c("iPSC_2197294_3x_vector", "NGN2_2197294_3x_vector",  "Sy5Y_2197294_3x_vector"),
  imagetype = "png",
  #filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_NGN2.png'),      ###### output in the working directory!!!
  filename = paste0(out_dir, '/venn_diagramm_compound_2186960_branapam_PTC518.3x_no_name_NGN2.png'),      ###### output in the working directory!!!
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


# heatmap
rm(list = ls())
library(pheatmap)
library(dplyr)
library("openxlsx")
#setwd("~/2T_Disk/Download/DATA")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/")

### 1st exp
DSG_3v3_1st = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/summarytable_with_gene-1st_experiment.csv", header = T)
DSG_3v3_1st[1:10, 1:10]
DSG_3v3_1st[1:10, ]
colnames(DSG_3v3_1st)
dim(DSG_3v3_1st)

row.names(DSG_3v3_1st) = DSG_3v3_1st$X
DSG_3v3_1st$X = NULL
dim(DSG_3v3_1st)

### 2nd exp
DSG_3v3_2nd = read.csv("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/summarytable_with_genedPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.csv", header = T)
DSG_3v3_2nd[1:10, 1:10]
DSG_3v3_2nd[1:10, ]
colnames(DSG_3v3_2nd)
dim(DSG_3v3_2nd)

row.names(DSG_3v3_2nd) = DSG_3v3_2nd$X
DSG_3v3_2nd$X = NULL
colnames(DSG_3v3_2nd)
dim(DSG_3v3_2nd)

getwd()
################# (1) data from 1st batch  === 5497 6152(1st) iPSC ===============
###== (1.1) 1755497
DSG_3v3_1st_iPSC_5497 =  DSG_3v3_1st[ , c("dPSI_TST11872_iPSC.BIO.1755497.12nM" , "dPSI_TST11872_iPSC.BIO.1755497.40nM"  ) ]    
DSG_3v3_1st_iPSC_5497 =  DSG_3v3_1st_iPSC_5497[complete.cases(DSG_3v3_1st_iPSC_5497), ] ######## complext case , mean both 3x and 10
dim(DSG_3v3_1st_iPSC_5497)
head(DSG_3v3_1st_iPSC_5497)

DSG_3v3_1st_iPSC_5497_mx  =  abs( as.matrix(DSG_3v3_1st_iPSC_5497) )
DSG_3v3_1st_iPSC_5497_df  =  as.data.frame(DSG_3v3_1st_iPSC_5497_mx)  
DSG_3v3_1st_iPSC_5497_df$Row.names = row.names(DSG_3v3_1st_iPSC_5497_df)
head(DSG_3v3_1st_iPSC_5497_df)

###== (1.2) 2006152
DSG_3v3_1st_iPSC_6152 =  DSG_3v3_1st[ , c( "dPSI_TST11872_iPSC.BIO.2006152.45nM" , "dPSI_TST11872_iPSC.BIO.2006152.150nM"  ) ]    
DSG_3v3_1st_iPSC_6152 =  DSG_3v3_1st_iPSC_6152[complete.cases(DSG_3v3_1st_iPSC_6152), ]
dim(DSG_3v3_1st_iPSC_6152)
head(DSG_3v3_1st_iPSC_6152)

DSG_3v3_1st_iPSC_6152_mx  =  abs( as.matrix(DSG_3v3_1st_iPSC_6152) )
DSG_3v3_1st_iPSC_6152_df  =  as.data.frame(DSG_3v3_1st_iPSC_6152_mx)  
DSG_3v3_1st_iPSC_6152_df$Row.names = row.names(DSG_3v3_1st_iPSC_6152_df)

################# (2) data from 2nd batch  ===  4088 4090 4748 5420 6152(2nd) 6527 6866  6960 9701  iPSC ====
###== (2.1) 4088
DSG_3v3_2nd_iPSC_4088 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.4088.3x", "dPSI_iPSC.4088.10x"    ) ]    
DSG_3v3_2nd_iPSC_4088 =  DSG_3v3_2nd_iPSC_4088[complete.cases(DSG_3v3_2nd_iPSC_4088), ]
dim(DSG_3v3_2nd_iPSC_4088)
head(DSG_3v3_2nd_iPSC_4088)

DSG_3v3_2nd_iPSC_4088_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_4088) )
DSG_3v3_2nd_iPSC_4088_df  =  as.data.frame(DSG_3v3_2nd_iPSC_4088_mx)  
DSG_3v3_2nd_iPSC_4088_df$Row.names = row.names(DSG_3v3_2nd_iPSC_4088_df)

###== (2.2) 4090
DSG_3v3_2nd_iPSC_4090 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.4090.3x", "dPSI_iPSC.4090.10x"   ) ]    
DSG_3v3_2nd_iPSC_4090 =  DSG_3v3_2nd_iPSC_4090[complete.cases(DSG_3v3_2nd_iPSC_4090), ]
dim(DSG_3v3_2nd_iPSC_4090)
head(DSG_3v3_2nd_iPSC_4090)

DSG_3v3_2nd_iPSC_4090_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_4090) )
DSG_3v3_2nd_iPSC_4090_df  =  as.data.frame(DSG_3v3_2nd_iPSC_4090_mx)  
DSG_3v3_2nd_iPSC_4090_df$Row.names = row.names(DSG_3v3_2nd_iPSC_4090_df)

###== (2.3) 4748
DSG_3v3_2nd_iPSC_4748 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.4748.3x", "dPSI_iPSC.4748.10x"    ) ]    
DSG_3v3_2nd_iPSC_4748 =  DSG_3v3_2nd_iPSC_4748[complete.cases(DSG_3v3_2nd_iPSC_4748), ]
dim(DSG_3v3_2nd_iPSC_4748)
head(DSG_3v3_2nd_iPSC_4748)

DSG_3v3_2nd_iPSC_4748_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_4748) )
DSG_3v3_2nd_iPSC_4748_df  =  as.data.frame(DSG_3v3_2nd_iPSC_4748_mx)  
DSG_3v3_2nd_iPSC_4748_df$Row.names = row.names(DSG_3v3_2nd_iPSC_4748_df)

###== (2.4) 5420 
DSG_3v3_2nd_iPSC_5420 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.5420.3x", "dPSI_iPSC.5420.10x"     ) ]    
DSG_3v3_2nd_iPSC_5420 =  DSG_3v3_2nd_iPSC_5420[complete.cases(DSG_3v3_2nd_iPSC_5420), ]
dim(DSG_3v3_2nd_iPSC_5420)
head(DSG_3v3_2nd_iPSC_5420)

DSG_3v3_2nd_iPSC_5420_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_5420) )
DSG_3v3_2nd_iPSC_5420_df  =  as.data.frame(DSG_3v3_2nd_iPSC_5420_mx)  
DSG_3v3_2nd_iPSC_5420_df$Row.names = row.names(DSG_3v3_2nd_iPSC_5420_df)

###== (2.5) 6152-2nd
DSG_3v3_2nd_iPSC_6152 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.6152.3x", "dPSI_iPSC.6152.10x" ) ]    
DSG_3v3_2nd_iPSC_6152 =  DSG_3v3_2nd_iPSC_6152[complete.cases(DSG_3v3_2nd_iPSC_6152), ]
dim(DSG_3v3_2nd_iPSC_6152)
head(DSG_3v3_2nd_iPSC_6152)

DSG_3v3_2nd_iPSC_6152_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_6152) )
DSG_3v3_2nd_iPSC_6152_df  =  as.data.frame(DSG_3v3_2nd_iPSC_6152_mx)  
DSG_3v3_2nd_iPSC_6152_df$Row.names = row.names(DSG_3v3_2nd_iPSC_6152_df)

# 4088 4090 4748 5420 6152(2nd) 6527 6866  6960 9701  ====
###== (2.6) 6527
DSG_3v3_2nd_iPSC_6527 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.6527.3x", "dPSI_iPSC.6527.10x"     ) ]    
DSG_3v3_2nd_iPSC_6527 =  DSG_3v3_2nd_iPSC_6527[complete.cases(DSG_3v3_2nd_iPSC_6527), ]
dim(DSG_3v3_2nd_iPSC_6527)
head(DSG_3v3_2nd_iPSC_6527)

DSG_3v3_2nd_iPSC_6527_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_6527) )
DSG_3v3_2nd_iPSC_6527_df  =  as.data.frame(DSG_3v3_2nd_iPSC_6527_mx)  
DSG_3v3_2nd_iPSC_6527_df$Row.names = row.names(DSG_3v3_2nd_iPSC_6527_df)

###== (2.7) 6866
DSG_3v3_2nd_iPSC_6866 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.6866.3x", "dPSI_iPSC.6866.10x"     ) ]    
DSG_3v3_2nd_iPSC_6866 =  DSG_3v3_2nd_iPSC_6866[complete.cases(DSG_3v3_2nd_iPSC_6866), ]
dim(DSG_3v3_2nd_iPSC_6866)
head(DSG_3v3_2nd_iPSC_6866)

DSG_3v3_2nd_iPSC_6866_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_6866) )
DSG_3v3_2nd_iPSC_6866_df  =  as.data.frame(DSG_3v3_2nd_iPSC_6866_mx)  
DSG_3v3_2nd_iPSC_6866_df$Row.names = row.names(DSG_3v3_2nd_iPSC_6866_df)

###== (2.8) 6960
#DSG_3v3_2nd_iPSC_6960 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.6960.3x", "dPSI_iPSC.6960.10x" ,  "dPSI_NGN2.6960.3x", "dPSI_NGN2.6960.10x"    ) ]    
DSG_3v3_2nd_iPSC_6960 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.6960.3x", "dPSI_iPSC.6960.10x"   ) ]    
DSG_3v3_2nd_iPSC_6960 =  DSG_3v3_2nd_iPSC_6960[complete.cases(DSG_3v3_2nd_iPSC_6960), ]
dim(DSG_3v3_2nd_iPSC_6960)
head(DSG_3v3_2nd_iPSC_6960)

DSG_3v3_2nd_iPSC_6960_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_6960) )
DSG_3v3_2nd_iPSC_6960_df  =  as.data.frame(DSG_3v3_2nd_iPSC_6960_mx)  
DSG_3v3_2nd_iPSC_6960_df$Row.names = row.names(DSG_3v3_2nd_iPSC_6960_df)

###== (2.9) 9701 
DSG_3v3_2nd_iPSC_9701 =  DSG_3v3_2nd[ , c( "dPSI_iPSC.9701.3x", "dPSI_iPSC.9701.10x"     ) ]    
DSG_3v3_2nd_iPSC_9701 =  DSG_3v3_2nd_iPSC_9701[complete.cases(DSG_3v3_2nd_iPSC_9701), ]
dim(DSG_3v3_2nd_iPSC_9701)
head(DSG_3v3_2nd_iPSC_9701)

DSG_3v3_2nd_iPSC_9701_mx  =  abs( as.matrix(DSG_3v3_2nd_iPSC_9701) )
DSG_3v3_2nd_iPSC_9701_df  =  as.data.frame(DSG_3v3_2nd_iPSC_9701_mx)  
DSG_3v3_2nd_iPSC_9701_df$Row.names = row.names(DSG_3v3_2nd_iPSC_9701_df)

################# (3) data from 1st batch  === 5497 6152(1st) NGN2 ===============
###== (3.1) 1755497
DSG_3v3_1st_NGN2_5497 =  DSG_3v3_1st[ , c("dPSI_TST11872_NGN2.BIO.1755497.12nM" , "dPSI_TST11872_NGN2.BIO.1755497.40nM"  ) ]    
DSG_3v3_1st_NGN2_5497 =  DSG_3v3_1st_NGN2_5497[complete.cases(DSG_3v3_1st_NGN2_5497), ]
dim(DSG_3v3_1st_NGN2_5497)
head(DSG_3v3_1st_NGN2_5497)

DSG_3v3_1st_NGN2_5497_mx  =  abs( as.matrix(DSG_3v3_1st_NGN2_5497) )
DSG_3v3_1st_NGN2_5497_df  =  as.data.frame(DSG_3v3_1st_NGN2_5497_mx)  
DSG_3v3_1st_NGN2_5497_df$Row.names = row.names(DSG_3v3_1st_NGN2_5497_df)

###== (3.2) 2006152
DSG_3v3_1st_NGN2_6152 =  DSG_3v3_1st[ , c( "dPSI_TST11872_NGN2.BIO.2006152.45nM" , "dPSI_TST11872_NGN2.BIO.2006152.150nM"  ) ]    
DSG_3v3_1st_NGN2_6152 =  DSG_3v3_1st_NGN2_6152[complete.cases(DSG_3v3_1st_NGN2_6152), ]
dim(DSG_3v3_1st_NGN2_6152)
head(DSG_3v3_1st_NGN2_6152)

DSG_3v3_1st_NGN2_6152_mx  =  abs( as.matrix(DSG_3v3_1st_NGN2_6152) )
DSG_3v3_1st_NGN2_6152_df  =  as.data.frame(DSG_3v3_1st_NGN2_6152_mx)  
DSG_3v3_1st_NGN2_6152_df$Row.names = row.names(DSG_3v3_1st_NGN2_6152_df)

################# (4) data from 2nd batch  ===  4088 4090 4748 5420 6152(2nd) 6527 6866  6960 9701  NGN2 ====
###== (4.1) 4088
DSG_3v3_2nd_NGN2_4088 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.4088.3x", "dPSI_NGN2.4088.10x"    ) ]    
DSG_3v3_2nd_NGN2_4088 =  DSG_3v3_2nd_NGN2_4088[complete.cases(DSG_3v3_2nd_NGN2_4088), ]
dim(DSG_3v3_2nd_NGN2_4088)
head(DSG_3v3_2nd_NGN2_4088)

DSG_3v3_2nd_NGN2_4088_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_4088) )
DSG_3v3_2nd_NGN2_4088_df  =  as.data.frame(DSG_3v3_2nd_NGN2_4088_mx)  
DSG_3v3_2nd_NGN2_4088_df$Row.names = row.names(DSG_3v3_2nd_NGN2_4088_df)

###== (4.2) 4090
DSG_3v3_2nd_NGN2_4090 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.4090.3x", "dPSI_NGN2.4090.10x"   ) ]    
DSG_3v3_2nd_NGN2_4090 =  DSG_3v3_2nd_NGN2_4090[complete.cases(DSG_3v3_2nd_NGN2_4090), ]
dim(DSG_3v3_2nd_NGN2_4090)
head(DSG_3v3_2nd_NGN2_4090)

DSG_3v3_2nd_NGN2_4090_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_4090) )
DSG_3v3_2nd_NGN2_4090_df  =  as.data.frame(DSG_3v3_2nd_NGN2_4090_mx)  
DSG_3v3_2nd_NGN2_4090_df$Row.names = row.names(DSG_3v3_2nd_NGN2_4090_df)
head(DSG_3v3_2nd_NGN2_4090_df)  

###== (4.3) 4748
DSG_3v3_2nd_NGN2_4748 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.4748.3x", "dPSI_NGN2.4748.10x"    ) ]    
DSG_3v3_2nd_NGN2_4748 =  DSG_3v3_2nd_NGN2_4748[complete.cases(DSG_3v3_2nd_NGN2_4748), ]
dim(DSG_3v3_2nd_NGN2_4748)
head(DSG_3v3_2nd_NGN2_4748)

DSG_3v3_2nd_NGN2_4748_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_4748) )
DSG_3v3_2nd_NGN2_4748_df  =  as.data.frame(DSG_3v3_2nd_NGN2_4748_mx)  
DSG_3v3_2nd_NGN2_4748_df$Row.names = row.names(DSG_3v3_2nd_NGN2_4748_df)

###== (4.4) 5420 
DSG_3v3_2nd_NGN2_5420 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.5420.3x", "dPSI_NGN2.5420.10x"     ) ]    
DSG_3v3_2nd_NGN2_5420 =  DSG_3v3_2nd_NGN2_5420[complete.cases(DSG_3v3_2nd_NGN2_5420), ]
dim(DSG_3v3_2nd_NGN2_5420)
head(DSG_3v3_2nd_NGN2_5420)

DSG_3v3_2nd_NGN2_5420_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_5420) )
DSG_3v3_2nd_NGN2_5420_df  =  as.data.frame(DSG_3v3_2nd_NGN2_5420_mx)  
DSG_3v3_2nd_NGN2_5420_df$Row.names = row.names(DSG_3v3_2nd_NGN2_5420_df)

###== (4.5) 6152-2nd
DSG_3v3_2nd_NGN2_6152 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.6152.3x", "dPSI_NGN2.6152.10x" ) ]    
DSG_3v3_2nd_NGN2_6152 =  DSG_3v3_2nd_NGN2_6152[complete.cases(DSG_3v3_2nd_NGN2_6152), ]
dim(DSG_3v3_2nd_NGN2_6152)
head(DSG_3v3_2nd_NGN2_6152)

DSG_3v3_2nd_NGN2_6152_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_6152) )
DSG_3v3_2nd_NGN2_6152_df  =  as.data.frame(DSG_3v3_2nd_NGN2_6152_mx)  
DSG_3v3_2nd_NGN2_6152_df$Row.names = row.names(DSG_3v3_2nd_NGN2_6152_df)

# 4088 4090 4748 5420 6152(4nd) 6527 6866  6960 9701  ====
###== (4.6) 6527
DSG_3v3_2nd_NGN2_6527 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.6527.3x", "dPSI_NGN2.6527.10x"     ) ]    
DSG_3v3_2nd_NGN2_6527 =  DSG_3v3_2nd_NGN2_6527[complete.cases(DSG_3v3_2nd_NGN2_6527), ]
dim(DSG_3v3_2nd_NGN2_6527)
head(DSG_3v3_2nd_NGN2_6527)

DSG_3v3_2nd_NGN2_6527_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_6527) )
DSG_3v3_2nd_NGN2_6527_df  =  as.data.frame(DSG_3v3_2nd_NGN2_6527_mx)  
DSG_3v3_2nd_NGN2_6527_df$Row.names = row.names(DSG_3v3_2nd_NGN2_6527_df)

###== (4.7) 6866
DSG_3v3_2nd_NGN2_6866 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.6866.3x", "dPSI_NGN2.6866.10x"     ) ]    
DSG_3v3_2nd_NGN2_6866 =  DSG_3v3_2nd_NGN2_6866[complete.cases(DSG_3v3_2nd_NGN2_6866), ]
dim(DSG_3v3_2nd_NGN2_6866)
head(DSG_3v3_2nd_NGN2_6866)

DSG_3v3_2nd_NGN2_6866_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_6866) )
DSG_3v3_2nd_NGN2_6866_df  =  as.data.frame(DSG_3v3_2nd_NGN2_6866_mx)  
DSG_3v3_2nd_NGN2_6866_df$Row.names = row.names(DSG_3v3_2nd_NGN2_6866_df)

###== (4.8) 6960
#DSG_3v3_2nd_NGN2_6960 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.6960.3x", "dPSI_NGN2.6960.10x" ,  "dPSI_NGN2.6960.3x", "dPSI_NGN2.6960.10x"    ) ]    
DSG_3v3_2nd_NGN2_6960 =  DSG_3v3_2nd[ , c(  "dPSI_NGN2.6960.3x", "dPSI_NGN2.6960.10x"    ) ]    
DSG_3v3_2nd_NGN2_6960 =  DSG_3v3_2nd_NGN2_6960[complete.cases(DSG_3v3_2nd_NGN2_6960), ]
dim(DSG_3v3_2nd_NGN2_6960)
head(DSG_3v3_2nd_NGN2_6960)

DSG_3v3_2nd_NGN2_6960_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_6960) )
DSG_3v3_2nd_NGN2_6960_df  =  as.data.frame(DSG_3v3_2nd_NGN2_6960_mx)  
DSG_3v3_2nd_NGN2_6960_df$Row.names = row.names(DSG_3v3_2nd_NGN2_6960_df)

###== (4.9) 9701 
DSG_3v3_2nd_NGN2_9701 =  DSG_3v3_2nd[ , c( "dPSI_NGN2.9701.3x", "dPSI_NGN2.9701.10x"     ) ]    
DSG_3v3_2nd_NGN2_9701 =  DSG_3v3_2nd_NGN2_9701[complete.cases(DSG_3v3_2nd_NGN2_9701), ]
dim(DSG_3v3_2nd_NGN2_9701)
head(DSG_3v3_2nd_NGN2_9701)

DSG_3v3_2nd_NGN2_9701_mx  =  abs( as.matrix(DSG_3v3_2nd_NGN2_9701) )
DSG_3v3_2nd_NGN2_9701_df  =  as.data.frame(DSG_3v3_2nd_NGN2_9701_mx)  
DSG_3v3_2nd_NGN2_9701_df$Row.names = row.names(DSG_3v3_2nd_NGN2_9701_df)



################# (5) data Transform ===============
head(DSG_3v3_1st_iPSC_5497_df)
head(DSG_3v3_1st_iPSC_6152_df)
head(DSG_3v3_2nd_iPSC_5420_df)
head(DSG_3v3_2nd_NGN2_5420_df)
head(DSG_3v3_2nd_NGN2_4748_df)
df_iPSC_list = list(DSG_3v3_1st_iPSC_5497_df[, c(1, ncol(DSG_3v3_1st_iPSC_5497_df)) ] ,
                    DSG_3v3_1st_iPSC_6152_df[, c(1, ncol(DSG_3v3_1st_iPSC_6152_df)) ] ,
                    DSG_3v3_2nd_iPSC_4088_df[, c(1, ncol(DSG_3v3_2nd_iPSC_4088_df)) ]  ,
                    DSG_3v3_2nd_iPSC_4090_df[, c(1, ncol(DSG_3v3_2nd_iPSC_4090_df)) ]  ,
                    DSG_3v3_2nd_iPSC_4748_df[, c(1, ncol(DSG_3v3_2nd_iPSC_4748_df)) ]  ,
                    DSG_3v3_2nd_iPSC_5420_df[, c(1, ncol(DSG_3v3_2nd_iPSC_5420_df)) ]  ,
                    DSG_3v3_2nd_iPSC_6152_df[, c(1, ncol(DSG_3v3_2nd_iPSC_6152_df)) ]  ,
                    DSG_3v3_2nd_iPSC_6527_df[, c(1, ncol(DSG_3v3_2nd_iPSC_6527_df)) ]  ,
                    DSG_3v3_2nd_iPSC_6866_df[, c(1, ncol(DSG_3v3_2nd_iPSC_6866_df)) ]  ,
                    DSG_3v3_2nd_iPSC_6960_df[, c(1, ncol(DSG_3v3_2nd_iPSC_6960_df)) ]  ,
                    DSG_3v3_2nd_iPSC_9701_df[, c(1, ncol(DSG_3v3_2nd_iPSC_9701_df)) ]  )
head(df_iPSC_list, 11)

DSG_3v3_iPSC = Reduce( function(x, y)  merge(x, y, by="Row.names", all=TRUE),  df_iPSC_list) # by=0 , Important

 dim(DSG_3v3_iPSC)
head(DSG_3v3_iPSC)
DSG_3v3_iPSC[1:2,1:11] 

df_NGN2_list = list(DSG_3v3_1st_NGN2_5497_df[, c(1, ncol(DSG_3v3_1st_NGN2_5497_df)) ] ,
                    DSG_3v3_1st_NGN2_6152_df[, c(1, ncol(DSG_3v3_1st_NGN2_6152_df)) ] ,
                    DSG_3v3_2nd_NGN2_4088_df[, c(1, ncol(DSG_3v3_2nd_NGN2_4088_df)) ]  ,
                    DSG_3v3_2nd_NGN2_4090_df[, c(1, ncol(DSG_3v3_2nd_NGN2_4090_df)) ]  ,
                    DSG_3v3_2nd_NGN2_4748_df[, c(1, ncol(DSG_3v3_2nd_NGN2_4748_df)) ]  ,
                    DSG_3v3_2nd_NGN2_5420_df[, c(1, ncol(DSG_3v3_2nd_NGN2_5420_df)) ]  ,
                    DSG_3v3_2nd_NGN2_6152_df[, c(1, ncol(DSG_3v3_2nd_NGN2_6152_df)) ]  ,
                    DSG_3v3_2nd_NGN2_6527_df[, c(1, ncol(DSG_3v3_2nd_NGN2_6527_df)) ]  ,
                    DSG_3v3_2nd_NGN2_6866_df[, c(1, ncol(DSG_3v3_2nd_NGN2_6866_df)) ]  ,
                    DSG_3v3_2nd_NGN2_6960_df[, c(1, ncol(DSG_3v3_2nd_NGN2_6960_df)) ]  ,
                    DSG_3v3_2nd_NGN2_9701_df[, c(1, ncol(DSG_3v3_2nd_NGN2_9701_df)) ]  )
head(df_NGN2_list, 11)

DSG_3v3_NGN2 = Reduce( function(x, y)  merge(x, y, by="Row.names", all=TRUE),  df_NGN2_list) # by=0 , Important

dim(DSG_3v3_NGN2)
head(DSG_3v3_NGN2)
DSG_3v3_NGN2[1:2,1:11] 

DSG_3v3_df         =  merge(DSG_3v3_iPSC, DSG_3v3_NGN2, by="Row.names", all=TRUE)
dim(DSG_3v3_df) 
colnames(DSG_3v3_df) 
head(DSG_3v3_df) 
tail(DSG_3v3_df) 
DSG_3v3_df_bk = DSG_3v3_df 
# BIO-1755497   # B-1755497
# BIO-2006152   # B-2006152
# BIO-6152      # B-6152
# BIO-2175420   # B-2175420
# BIO-2139701   # B-2139701
# BIO-2174748   #	B-2174748
# BIO-2176866   # B-2176866
# BIO-2186960   # B-2186960
# BIO-2186527   # B-2186527
# BIO-2184088   # B-2184088
# BIO-2184090	  # B-2184090
DSG_3v3_df = DSG_3v3_df[                                        , c(  "Row.names", 
                                                                      "dPSI_TST11872_iPSC.BIO.1755497.12nM", 
                                                                      "dPSI_TST11872_iPSC.BIO.2006152.45nM",
                                                                      "dPSI_iPSC.6152.3x" ,
                                                                      "dPSI_iPSC.5420.3x" ,
                                                                      "dPSI_iPSC.9701.3x" ,           
                                                                      "dPSI_iPSC.4748.3x" ,
                                                                      "dPSI_iPSC.6866.3x" ,
                                                                      "dPSI_iPSC.6960.3x" ,
                                                                      "dPSI_iPSC.6527.3x" ,
                                                                      "dPSI_iPSC.4088.3x" ,
                                                                      "dPSI_iPSC.4090.3x" ,
                                                                      
                                                                      "dPSI_TST11872_NGN2.BIO.1755497.12nM",
                                                                      "dPSI_TST11872_NGN2.BIO.2006152.45nM",
                                                                      "dPSI_NGN2.6152.3x" ,
                                                                      "dPSI_NGN2.5420.3x" ,
                                                                      "dPSI_NGN2.9701.3x" ,     
                                                                      "dPSI_NGN2.4748.3x" ,
                                                                      "dPSI_NGN2.6866.3x" ,
                                                                      "dPSI_NGN2.6960.3x" ,   
                                                                      "dPSI_NGN2.6527.3x" ,    
                                                                      "dPSI_NGN2.4088.3x" ,
                                                                      "dPSI_NGN2.4090.3x" ) ]

DSG_3v3_df = DSG_3v3_df[ order(DSG_3v3_df$Row.names), ]
head(DSG_3v3_df)
dim(DSG_3v3_df)

row.names(DSG_3v3_df) = DSG_3v3_df$Row.names
DSG_3v3_df$Row.names  = NULL

DSG_3v3_iPSC_df = DSG_3v3_df[, c(1:11)]
DSG_3v3_NGN2_df = DSG_3v3_df[, c(12:22)]

Total_Frenq  = apply( DSG_3v3_df,      1, function(DSG_3v3_df)       length(na.omit(DSG_3v3_df      ) ) )
iPSC_Frenq   = apply( DSG_3v3_iPSC_df, 1, function(DSG_3v3_iPSC_df)  length(na.omit(DSG_3v3_iPSC_df ) ) )
NGN2_Frenq   = apply( DSG_3v3_NGN2_df, 1, function(DSG_3v3_NGN2_df)  length(na.omit(DSG_3v3_NGN2_df ) ) )

# https://slowkow.com/notes/pheatmap-tutorial/
# colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))
# https://stackoverflow.com/questions/48176640/rns-seq-analysis-using-r-make-pheatmap-all-looks-blue
# https://www.r-bloggers.com/2010/12/r-using-rcolorbrewer-to-colour-your-figures-in-r/
# http://www.sthda.com/english/wiki/colors-in-r
library(RColorBrewer)
# colors <- brewer.pal(9, "Blues")
# colors <- colorRampPalette(brewer.pal(9,"Blues"))(100)
# colors <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
colors <- colorRampPalette(brewer.pal(9,"PuBu"))(100)
# colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
# pheatmap(mat,          col=colors,)

# ========= heatmap and clustering ========== #

head(DSG_3v3_iPSC_df )
DSG_3v3_iPSC_df_1 = DSG_3v3_iPSC_df
DSG_3v3_iPSC_df_1[is.na(DSG_3v3_iPSC_df_1)] = 0

head(DSG_3v3_iPSC_df_1 )

pheatmap(DSG_3v3_iPSC_df_1 )
pheatmap(DSG_3v3_iPSC_df_1 , color = colors )

xx = pheatmap(DSG_3v3_iPSC_df_1, color = colors  )

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(xx, "DSG_iPSC_compounds_heatmap_clustring1.pdf", 9, 24)

# install.packages("dendsort")
# library(dendsort)
# 
# sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
# 
# mat_cluster_cols <- sort_hclust( DSG_3v3_iPSC_df_1 )# mat_cluster_cols <- sort_hclust(mat_cluster_cols)
# plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

mat_cluster_cols <- hclust(dist(t(DSG_3v3_iPSC_df_1)))
plot(mat_cluster_cols, main = "DSG_iPSC_compounds_clustring_Dendrogram", xlab = "", sub = "")
# save as pdf landscape 6 X 10  
?colorRampPalette

# ========= heatmap and clustering ========== #
  
#DSG_3v3_df$Total_Frenq  = apply(DSG_3v3_df, 1, function(DSG_3v3_df) length(na.omit(DSG_3v3_df) ) )


# colnames(DSG_3v3_df) =  c(      "1st_iPSC.1755497.12nM",
#                                 "1st_iPSC.2006152.45nM",
#                                 "dPSI_iPSC.6152.3x" ,
#                                 "dPSI_iPSC.5420.3x" ,
#                                 "dPSI_iPSC.9701.3x" ,       
#                                 "dPSI_iPSC.4748.3x" ,
#                                 "dPSI_iPSC.6866.3x" ,
#                                 "dPSI_iPSC.6960.3x" ,
#                                 "dPSI_iPSC.6527.3x" ,
#                                 "dPSI_iPSC.4088.3x" ,
#                                 "dPSI_iPSC.4090.3x" ,
#                                 "1st_NGN2.1755497.12nM",
#                                 "1st_NGN2.2006152.45nM",
#                                 "dPSI_NGN2.6152.3x" ,
#                                 "dPSI_NGN2.5420.3x" ,
#                                 "dPSI_NGN2.9701.3x" ,
#                                 "dPSI_NGN2.4748.3x" ,
#                                 "dPSI_NGN2.6866.3x" ,
#                                 "dPSI_NGN2.6960.3x" ,
#                                 "dPSI_NGN2.6527.3x" ,
#                                 "dPSI_NGN2.4088.3x" ,
#                                 "dPSI_NGN2.4090.3x" )
                                #"Total_Frenq") 

#DSG_3v3_df <- DSG_3v3_df%>% select(Total_Frenq, everything())

DSG_3v3_df_with_freq <-  cbind( iPSC_Frenq ,NGN2_Frenq,Total_Frenq , DSG_3v3_df)

head(DSG_3v3_df_with_freq)

write.csv( DSG_3v3_df_with_freq ,file = paste0("DSG_3v3_df_with_freq_11_11", ".csv"), quote = F, row.names = T)

# df[order(df$var1), ]
# #sort data frame by multiple columns alphabetically 
# df[with(df, order(var1, var2)), ]
