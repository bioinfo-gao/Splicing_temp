
#index
FILES=/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/bam/point125bam/*1755497-3x-?.Aligned.sortedByCoord.out.bam

for file in $FILES ; do 
      #samtools view -bo  /camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/bam/point125bam/$(basename $file ) -s 125.125 $file 
      samtools index $file  $outfile ======
done

# subset
# https://stackoverflow.com/questions/8512462/looping-through-all-files-in-a-directory
# for file in *; do 
#     if [ -f "$file" ]; then 
#         echo "$file" 
#     fi 
# done

# https://bioinformatics.stackexchange.com/questions/3565/subset-smaller-bam-to-contain-several-thousand-rows-from-multiple-chromosomes
#samtools view -bo subset.bam -s 123.4 alignments.bam chr1 chr2 # select 40% (the .4 part) of the reads (123 i
#ls /camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/bam/*1755497-10x-?.Aligned.sortedByCoord.out.bam
#!/bin/bash
FILES=/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/bam/*1755497-10x-?.Aligned.sortedByCoord.out.bam

for files in $FILES ; do 
#    if [ -f "$file" ]; then 
      #echo "$file" 
      ls $file 
#  fi 
done
# 
samtools view -bo subset.bam -s 123.5 alignments.bam

#setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/code/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/code_and_genome_bam/script")
#link_content = read.csv("full_record_02_25.txt", header = F) 
link_content = read.csv("full_2022_02_28_total.txt", header = F) # link_content = read.csv("2022_02_26_total.txt", header = F) #link_content = read.csv("full_record_02_25_bam_selected.txt", header = F) 
link_content
str(link_content)
dim(link_content)

file_list = link_content[grepl('sortedByCoord.out.bam$', link_content$V1),  ] # ending file_list = link_content[grepl('.bam$', link_content$V1),  ] # ending
file_list
write.csv(file_list, file = "full_2_28_total_bam.txt", quote=F, sep = ",", row.names = F, col.names = F) #=======<<<<<<<<<<


for(i in 1:length(file_list) ){
  #i=2
  fileConn <- file(paste0( i, "_bam_path.bash")) #   fileConn <- file(paste0( i, "_bam_path.bash"))
  #fileConn <- file(paste0( i, "_bam_path.txt")) #   fileConn <- file(paste0( i, "_bam_path.bash"))
  writeLines(file_list[i] , fileConn) # -i download all files in the files 
  #writeLines(paste0( "nohup wget  ", file_list[i] , "  &" ), fileConn) # -i download all files in the files 
  #fileConn <- file("output.txt")
  #writeLines(file_list[i], file = paste0( i, "_bam_path.txt"), quote=F, sep = "\t", row.names = F, col.names = F)
  # https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r 
  # another good method is use sink()
}


setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Christina_RASLseq/")
getwd()
# rank                                                      Exon rank in transcript    structure
genes = c("PDE7A",  "POMT2",  "EVC",  "DENND1A",  "WDR27",  "KCNT2",  "DLGAP4",   "DHFR",  
          "PITPNB",  "XRN2",  "PDXDC1",  "FHOD3",  "KCNT2",  "CIP2A",  
          "TENT2",  "HLTF",  "HTR3A",  "GNAQ",  "KDSR",  "L3MBTL2",  "SLC7A6",
          "ADAMTS19",  "RDX",  "DLG5",  "BTBD10",  "FBL",  "ASNS",  "AGPS",  "MLLT10",  "RAD21",  "ELP4",  "VPS41")

## ----useEnsembl-----------------------------------------------------------------------------------
library(biomaRt)
listEnsembl()

mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
mart
## ----listDatasets---------------------------------------------------------------------------------
listDatasets( mart )
listAttributes( mart) # ====<<<< # ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#getAtt <- c('external_gene_name', 'ensembl_gene_id',  'ensembl_transcript_id', 'ensembl_transcript_id_version', 'chromosome_name', 'start_position', 'end_position', 'ensembl_exon_id', 'strand','exon_chrom_start','exon_chrom_end',  'rank')
getAtt <- c('external_gene_name', 'ensembl_gene_id','chromosome_name','strand')


elocs = list()

for ( i in 1:length(genes)) { # for ( i in 14:15) {for ( i in 16:25) { # for ( i in 14:15) {
  #i=1   #
  #i = 11 
  i 
  print(i)
  elocs[[i]] = getBM(attributes=getAtt, 
                     filters="hgnc_symbol", 
                     value=genes[i],
                     mart=mart )
  print( paste0( "=========================  (      ", i,  "      ) =========================  ") ) 
  print( paste0( "=========================  (  ", genes[i],  "  ) =========================  ") ) 
  print( elocs[[i]] )
  i = i+1
}
# i =1 PDE7A
# 53              PDE7A ENSG00000205268       ENST00000519231             ENST00000519231.5               8       65714334     65842322     -1         65779720       65779803 ENSE00003654459    5
# 54              PDE7A ENSG00000205268       ENST00000519231             ENST00000519231.5               8       65714334     65842322     -1         65747705 ----- 65747803 ENSE00002125114    6
# 55              PDE7A ENSG00000205268       ENST00000523253             ENST00000523253.1               8       65714334     65842322     -1         65782783       65782843 ENSE00003615863    3
# i 2 POMT2
# 15              POMT2 ENSG00000009830       ENST00000261534             ENST00000261534.9              14       77274956     77320883     -1         77288762       77288831 ENSE00003644918   11
# 16              POMT2 ENSG00000009830       ENST00000261534             ENST00000261534.9              14       77274956     77320883     -1         77286744       77286822 ENSE00003503606   12
# == 3
# 5                 EVC ENSG00000072840       ENST00000264956            ENST00000264956.11               4        5711201      5814305      1          5733351        5733435 ENSE00001008086    5
# 6                 EVC ENSG00000072840       ENST00000264956            ENST00000264956.11               4        5711201      5814305      1          5741716        5741814 ENSE00001193057    6
# ====4 
# 12            DENND1A ENSG00000119522       ENST00000473039             ENST00000473039.5               9      123379654    123930152     -1        123452276      123452347 ENSE00003485358   12
# 34            DENND1A ENSG00000119522       ENST00000373624             ENST00000373624.6               9      123379654    123930152     -1        123454739      123454779 ENSE00003509546   16
=== 5
19              WDR27 ENSG00000184465       ENST00000448612             ENST00000448612.6               6      169457212    169702067     -1        169662304      169662424 ENSE00003659062    9
20              WDR27 ENSG00000184465       ENST00000448612             ENST00000448612.6               6      169457212    169702067     -1        169660663      169660766 ENSE00003480055   10
21              WDR27 ENSG00000184465       ENST00000448612             ENST00000448612.6               6      169457212    169702067     -1        169659451      169659518 ENSE00003685937   11
=== 6
19              KCNT2 ENSG00000162687       ENST00000367433             ENST00000367433.9               1      196225779    196609225     -1        196326717      196326889 ENSE00001315063   19
20              KCNT2 ENSG00000162687       ENST00000367433             ENST00000367433.9               1      196225779    196609225     -1        196315892      196316026 ENSE00003474978   20
=== 7
20             DLGAP4 ENSG00000080845       ENST00000339266            ENST00000339266.10              20       36306336     36528637      1         36496705       36497066 ENSE00003684659    8
21             DLGAP4 ENSG00000080845       ENST00000339266            ENST00000339266.10              20       36306336     36528637      1         36500199       36500611 ENSE00003488351   10
25             DLGAP4 ENSG00000080845       ENST00000339266            ENST00000339266.10              20       36306336     36528637      1         36499588       36499676 ENSE00003471293    9
=== 8
10               DHFR ENSG00000228716       ENST00000505337             ENST00000505337.5               5       80626226     80654983     -1         80633877       80633992 ENSE00003465725    5
12               DHFR ENSG00000228716       ENST00000505337             ENST00000505337.5               5       80626226     80654983     -1         80629039       80629165 ENSE00002082987    6
  
=== 9 
7              PITPNB ENSG00000180957       ENST00000335272            ENST00000335272.10              22       27851669     27920134 ENSE00003693284     -1         27894555       27894638    7
8              PITPNB ENSG00000180957       ENST00000335272            ENST00000335272.10              22       27851669     27920134 ENSE00003623733     -1         27873738       27873815    8
=== 10
16               XRN2 ENSG00000088930       ENST00000377191             ENST00000377191.5              20       21303331     21389825      1         21344090       21344208 ENSE00000660879   16
17               XRN2 ENSG00000088930       ENST00000377191             ENST00000377191.5              20       21303331     21389825      1         21346415       21346550 ENSE00000564389   17
# === 11 ????? Latest Assembly
# chr16:14,974,591-15,153,218
# (GRCh38/hg38)
# Size: 178,628 bases Orientation: Plus strand
# 
# Previous Assembly
# chr16:15,068,592-15,233,267
# (GRCh37/hg19 by NCBI Gene)
# Size: 164,676 bases Orientation: Plus strand
# 
# chr16:15,068,448-15,233,196
# (GRCh37/hg19 by Ensembl)
# Size: 164,749 bases Orientation: Plus strand
=== 12
18              FHOD3 ENSG00000134775       ENST00000257209             ENST00000257209.8              18       36297713     36780220      1         36740656       36740838 ENSE00003722930   18
19              FHOD3 ENSG00000134775       ENST00000257209             ENST00000257209.8              18       36297713     36780220      1         36742737       36742856 ENSE00000916123   19
===
11              CIP2A ENSG00000163507       ENST00000295746            ENST00000295746.13               3      108549864    108589644     -1        108566497      108566638 ENSE00003505942   11
12              CIP2A ENSG00000163507       ENST00000295746            ENST00000295746.13               3      108549864    108589644     -1        108565355      108565454 ENSE00003595116   12
===
12              TENT2 ENSG00000164329       ENST00000423041             ENST00000423041.6               5       79612120     79688246      1         79641105       79641196 ENSE00001083575    7
13              TENT2 ENSG00000164329       ENST00000423041             ENST00000423041.6               5       79612120     79688246      1         79642844       79642910 ENSE00001643516    8
=== 16
37               HLTF ENSG00000071794       ENST00000310053            ENST00000310053.10               3      149030127    149086554     -1        149055303      149055400 ENSE00003570796   14
38               HLTF ENSG00000071794       ENST00000310053            ENST00000310053.10               3      149030127    149086554     -1        149050232      149050375 ENSE00000779435   15  
=== 17
24              HTR3A ENSG00000166736       ENST00000375498             ENST00000375498.6              11      113975075    113990313      1        113986825      113987046 ENSE00003521883    8
25              HTR3A ENSG00000166736       ENST00000375498             ENST00000375498.6              11      113975075    113990313      1        113975075      113975392 ENSE00001106099    1
26              HTR3A ENSG00000166736       ENST00000375498             ENST00000375498.6              11      113975075    113990313      1        113989465      113990310 ENSE00002086203    9  
=== 18
8                GNAQ ENSG00000156052       ENST00000411677             ENST00000411677.1               9       77716097     78031811     -1         77922161       77922345 ENSE00001773857    2
9                GNAQ ENSG00000156052       ENST00000411677             ENST00000411677.1               9       77716097     78031811     -1         77815616       77815770 ENSE00003536290    3
=== 19
7                KDSR ENSG00000119537       ENST00000645214             ENST00000645214.2              18       63327726     63367228     -1         63344410       63344493 ENSE00003524868    7
8                KDSR ENSG00000119537       ENST00000645214             ENST00000645214.2              18       63327726     63367228     -1         63338800       63338883 ENSE00003569694    8
=== 20
16            L3MBTL2 ENSG00000100395       ENST00000489136             ENST00000489136.5              22       41205282     41231271      1         41216139       41216262 ENSE00003466954    5 **
23            L3MBTL2 ENSG00000100395       ENST00000466589             ENST00000466589.5              22       41205282     41231271      1         41217123       41217202 ENSE00003638297    5 **
24            L3MBTL2 ENSG00000100395       ENST00000466589             ENST00000466589.5              22       41205282     41231271      1         41219419       41219536 ENSE00003600505    6 **
=== 21
26             SLC7A6 ENSG00000103064       ENST00000219343            ENST00000219343.11              16       68264516     68301823      1         68274691       68275249 ENSE00001310370    3
27             SLC7A6 ENSG00000103064       ENST00000219343            ENST00000219343.11              16       68264516     68301823      1         68287746       68287871 ENSE00003637599    4
=== 22
22           ADAMTS19 ENSG00000145808       ENST00000274487             ENST00000274487.9               5      129460281    129738683      1        129684120      129684273 ENSE00001172271   18 *
33           ADAMTS19 ENSG00000145808       ENST00000509467             ENST00000509467.1               5      129460281    129738683      1        129688100      129688214 ENSE00002058109    1
34           ADAMTS19 ENSG00000145808       ENST00000509467             ENST00000509467.1               5      129460281    129738683      1        129694720      129694855 ENSE00003469389    2
=== 23
61                RDX ENSG00000137710       ENST00000647231             ENST00000647231.1              11      109864295    110296712     -1        110258106      110258189 ENSE00003572576    6
60                RDX ENSG00000137710       ENST00000647231             ENST00000647231.1              11      109864295    110296712     -1        110263960      110264234 ENSE00003658656    5
59                RDX ENSG00000137710       ENST00000647231             ENST00000647231.1              11      109864295    110296712     -1        110264779      110264874 ENSE00003588483    4
58                RDX ENSG00000137710       ENST00000647231             ENST00000647231.1              11      109864295    110296712     -1        110272536      110272619 ENSE00003664757    3
57                RDX ENSG00000137710       ENST00000647231             ENST00000647231.1              11      109864295    110296712     -1        110279681      110279756 ENSE00003482846    2

=== 24
8                DLG5 ENSG00000274429       ENST00000632221             ENST00000632221.1 CHR_HSCHR10_1_CTG4       77790791     77926526     -1         77816551       77816701 ENSE00003746703    8
9                DLG5 ENSG00000274429       ENST00000632221             ENST00000632221.1 CHR_HSCHR10_1_CTG4       77790791     77926526     -1         77812215       77812377 ENSE00003742338    9
=== 25
20             BTBD10 ENSG00000148925       ENST00000530907             ENST00000530907.5              11       13388008     13463297     -1         13419460       13419745 ENSE00003612602    3
21             BTBD10 ENSG00000148925       ENST00000530907             ENST00000530907.5              11       13388008     13463297     -1         13417158       13417260 ENSE00003495010    4
=== 26
17                FBL ENSG00000280548       ENST00000630901             ENST00000630901.2 CHR_HG2021_PATCH       39834458     39846379     -1         39837711       39837843 ENSE00003768231    2
18                FBL ENSG00000280548       ENST00000630901             ENST00000630901.2 CHR_HG2021_PATCH       39834458     39846379     -1         39836556       39836668 ENSE00003764043    3
=== 27
31               ASNS ENSG00000070669       ENST00000437628             ENST00000437628.5               7       97851677     97872542     -1         97855353       97855459 ENSE00000707047    7
32               ASNS ENSG00000070669       ENST00000437628             ENST00000437628.5               7       97851677     97872542     -1         97854580       97854680 ENSE00000707046    8
=== 28
2                AGPS ENSG00000018510       ENST00000680705             ENST00000680705.1               2      177392746    177567024      1        177420269      177420358 ENSE00003913213    2
3                AGPS ENSG00000018510       ENST00000680705             ENST00000680705.1               2      177392746    177567024      1        177434327      177434417 ENSE00003914986    3
=== 29
22             MLLT10 ENSG00000078403       ENST00000631589             ENST00000631589.1              10       21524646     21743630      1         21727856       21727928 ENSE00001803822   16
23             MLLT10 ENSG00000078403       ENST00000631589             ENST00000631589.1              10       21524646     21743630      1         21730900       21731054 ENSE00002719058   17
=== 30
14              RAD21 ENSG00000164754       ENST00000687358             ENST00000687358.1               8      116845934    116874776     -1        116861841      116861940 ENSE00003464321    4
15              RAD21 ENSG00000164754       ENST00000687358             ENST00000687358.1               8      116845934    116874776     -1        116858352      116858458 ENSE00003579379    5
=== 31
45               ELP4 ENSG00000109911       ENST00000350638            ENST00000350638.10              11       31509755     31790324      1         31627110       31627194 ENSE00003616769    6
46               ELP4 ENSG00000109911       ENST00000350638            ENST00000350638.10              11       31509755     31790324      1         31632217       31632405 ENSE00003580122    7
=== 32
7               VPS41 ENSG00000006715       ENST00000310301             ENST00000310301.9               7       38722974     38932394     -1         38817817       38817882 ENSE00000460377    7
8               VPS41 ENSG00000006715       ENST00000310301             ENST00000310301.9               7       38722974     38932394     -1         38796745       38796864 ENSE00003787385    8
=== 33
  
# elocs_order = elocs [ with(elocs, order(external_gene_name, exon_chrom_start,  exon_chrom_end, ensembl_exon_id, rank) ),  ]  # unique_elocs_order = unique_elocs [ with(unique_elocs, order(external_gene_name, rank) ),  ] 
# elocs_order 
# elocs_order[!duplicated(elocs_order), ]
# unique_elocs_order  <- elocs_order[!duplicated(elocs_order), ]
# unique_elocs_order 


#unique_elocs = unique(elocs)

  
#   https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
#   # Sort by vector name [z] then [x]
#   dataframe[
#     with(dataframe, order(z, x)),
#   ]
# Similarly, to sort by multiple columns based on column index, add additional arguments to order() with differing indices:
#   
#   # Sort by column index [1] then [3]
#   dataframe[
#     order( dataframe[,1], dataframe[,3] ),
#   ]  


library(leafviz)
library(stringr)


workding_dirs = list.dirs("/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run", recursive=F) # no argument pattern= "TST11872", also list.files

for ( dir in workding_dirs ) {
   #rm(list=setdiff(ls(), "x")) #https://stackoverflow.com/questions/6190051/how-can-i-remove-all-objects-but-one-from-the-workspace-in-r
   raw_Rdata_file  = paste0(dir, "/leafcutter/leafviz.Rdata"           )
   copy_Rdata_file = paste0(dir, "/leafcutter/CopyOfRaw_leafviz.Rdata" )
   print(raw_Rdata_file)
   # file.copy(raw_Rdata_file, copy_Rdata_file )
   load(raw_Rdata_file)
    dim(counts)
   head(counts)
   temp_intron_data =  row.names(counts)
   row.names(counts) = str_sub(temp_intron_data , 1 , nchar(temp_intron_data)-2 ) #stringr chrX:1390290:1391899:clu_2:-   remove the trailing : and -, troublesome
   save.image(raw_Rdata_file)
}

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-12nM/leafcutter/leafviz.Rdata")
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1755497-40nM/leafcutter/leafviz.Rdata")
 
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-3uM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-1949634-10uM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-45nM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2006152-150nM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-261nM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2059811-870nM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-3uM/leafcutter/leafviz.Rdata"  )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060573-10uM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-1000nM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2060884-3700nM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-1680nM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2070692-5600nM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-300nM/leafcutter/leafviz.Rdata" )
# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2135644-1uM/leafcutter/leafviz.Rdata" )

# leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-300nM/leafcutter/leafviz.Rdata"  )
leafviz( "/camhpc/home/zgao1/NGS_projects/TST11872/dnanexus/20220204181515_zhen.gao/Result/Analysis_2_2022_02_try_run/TST11872_iPSC-BIO-2136770-1000nM/leafcutter/leafviz.Rdata" )




fastqc_dir="/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/fastqc"

echo "source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/mryals/.conda/envs/multiqc;cd $fastqc_dir; multiqc $fastqc_dir" | qsub -N multiqc -q cpu.q -l h_rt=272:00:00 -l h_vmem=256G -pe thread 12


etwd("/edgehpc/dept/compbio/users/zgao1/Procedure/Spelicing/Old_Before_Edge/") # all RNAseq
setwd("/edgehpc/dept/compbio/users/zgao1/project_RNAseq/TST12145/") #
setwd("/edgehpc/dept/compbio/users/zgao1/project_RNAseq/TST12086/") # all RNAseq
setwd("/edgehpc/dept/compbio/users/zgao1/project_RNAseq/TST11872/") # all RNAseq
setwd("/edgehpc/dept/compbio/users/zgao1/project_RNAseq/TST11955/") # all RNAseq

setwd(    "/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DNANenxs_Master/")
setwd(    "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DNANenxs_Master/")
"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/From_DNANenxs_Master_3vs3_09-25/"
"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/From_DNANenxs_Master_3vs3_09-25/"

"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/DSG_Result/splicing_offtargets_rerun_4vs4_09-25/"
"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/DSG_Result/splicing_offtargets_rerun_4vs4_09-25/"

"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/splicing_offtargets_3vs3_09-25/"
"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/splicing_offtargets_3vs3_09-25/"

"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921054434_Zhen.Gao_GMfibro_3vs3/DSG_Result/splicing_offtargets_3vs3_09-27_abs_dPSI"
"/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230921051846_Zhen.Gao_ShSy5Y_3vs3/DSG_Result/splicing_offtargets_3vs3_09-27_abs_dPSI/"


setwd("/mnt/depts/dept04/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/GM_Full_DEG_Result_reanalysis_09_29_3vs3/DEG_plots")

setwd(    "/edgehpc/dept/compbio/projects/TST12188/splicing_events/")
setwd(    "/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/GM_Full_DEG_Result_reanalysis_09_07")
setwd("/mnt/depts/dept04/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/Sh_Full_DEG_Result_reanalysis_09_07")

setwd("/mnt/depts/dept04/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/DSG_Result/From_DNANenxs_Master_analysis.02.splicing_offtargets_09_12")
setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/DSG_Result/From_DNANenxs_Master_analysis.02.splicing_offtargets_09_12/")

setwd("/edgehpc/dept/compbio/projects/TST12188/")      # ===> Christina SMN2 SM selectivity profiling- round 1
setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/Analysis_Code/DSG/2023-10-05-compounds-correlations/")      # ===> Christina SMN2 SM selectivity profiling- round 1

setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/")      # ===> Christina SMN2 SM selectivity profiling- round 1
setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230823221304_Zhen.Gao_ShSy5Y_GRCh38/bam/")      # ===> Christina SMN2 SM selectivity profiling- round 1
setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/20230824054737_Zhen.Gao_GM_noSMN1/bam/")      # ===> Christina SMN2 SM selectivity profiling- round 1
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/") # Joon
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master_analysis.02.splicing_offtargets_3vs3/") # Joon
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/2023_All_data_code_RNA3.06/")
setwd("/edgehpc/dept/compbio/projects/TST11872/dnanexus/20220204181515_zhen.gao/")
setwd("/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/")
setwd("/edgehpc/dept/compbio/genomes/human/gencode.v34/fasta/")



setwd("/edgehpc/dept/compbio/projects/TST11354/motif/motif_rmats_SE_nearExon_dPSI0.1_pm5bp_matchedbG/") # run.04.from_bed_to_sequence.bash
#cp -r /edgehpc/dept/compbio/projects/TST11354/motif/  /edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Dann_motif # copy and reName at the same time
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Dann_motif/motif_rmats_SE_nearExon_dPSI0.1_pm5bp_matchedbG")

# Joon Lee 05/11/2023
# to Magnus Pfaffenbach Andreas Weihofen;Colin Choi
# RE: RNAseq 
setwd("/edgehpc/dept/compbio/projects/TST11451") # iPSC-MN /edgehpc/dept/compbio/projects/TST11354/ /mnt/depts/dept04/compbio/users/zgao1/project_splicig/Email_Dann.txt
setwd("/edgehpc/dept/compbio/projects/TST11451/master_table/res_comb_HNDS0030_01-RG7961_120nM_"  )    #EC50 = 40nM 3x = 120 nm accroding to Joon Email 
setwd("/edgehpc/dept/compbio/projects/TST11451/master_table/res_comb_HNDS006_01C2-RG7961_120nM_" )
setwd("/edgehpc/dept/compbio/projects/TST11451/master_table/res_comb_GM24468D-RG7961_120nM_"     )



setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/DSG_Result/From_DNANenxs_Master/")
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/Local_Call_DSG_code/")
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/SJ/")
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/20230719034037/")

setwd("/edgehpc/dept/compbio/projects/TST12188/dnanexus/")
setwd("/edgehpc/dept/compbio/projects/TST12145/dnanexus/20230622062158_Zhen.Gao/")
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/")
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/2023_All_data_code_RNA3.06/")
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_1755497_10x/")

setwd("/edgehpc/dept/compbio/projects/TST11872/dnanexus/20220204181515_zhen.gao/")
setwd("/edgehpc/dept/compbio/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/")

setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/EA_DEG_automatic_files/")

setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig/ggplot2_exploring_interim/")
setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_DEG_for_Magnus_06_01/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig/")

setwd("/mnt/depts/dept04/compbio/users/zgao1/Procedure/DNANexus_procedures//")           #
setwd("/mnt/depts/dept04/compbio/users/zgao1/Procedure/MulitiOmics/")           #
setwd("/mnt/depts/dept04/compbio/users/zgao1/Procedure/multiqc/")           #
setwd("/mnt/depts/dept04/compbio/users/zgao1/Procedure/Spelicing/")           #
setwd("/mnt/depts/dept04/compbio/users/zgao1/project_RNAseq/TST12086/")           #

setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/rank_based")
# setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/DEG/code/The_1st_DEG_big_table.R")

setwd("camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DSG_ranking_and_join" )
#   $  ls
# 3rd_catherine1.R  3rd_catherine4.R  3rd_catherine.R  3rd_RNAseq.R  overlap_1st_v_2nd_1.R  overlap_1st_v_2nd_2_for_catherine3.R

setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_dPSI_a_sig/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/code_downstream/") # ===========>>>> run.04.DSG_counts_plot_ranking_batch123_DSG_count_Magnus_Dann_final_5.R 


setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/Result/iPSC_2189972_3x//")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/2023_All_data_code_RNA3.06/")

setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/code/TST12086/2023_modified_code_used/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/Result/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/genome_bam")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/code/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/bam/")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/code_and_genome_bam/script")
setwd("/camhpc/ngs/projects/TST12086/dnanexus/20230219201433_Zhen.Gao/genome_bam/bad_temp")

"/edgehpc/dept/compbio/projects/RASL_seq_BSSI/programs/dev/TST12036_downsample10/711_TST12036_downsample10_RASLseq_gene_QC_EC50_indiv_FC_comparison.html" # Matthew Ryals demonstation of RASL-seq

setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230305034844_Zhen.Gao") # setwd("/edgehpc/dept/compbio/projects/TST12086/dnanexus/20230218080338_Zhen.Gao/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST12806/full_record/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST12806/full_record/") # /edgehpc/dept/compbio/users/zgao1/project_splicing/TST12806/All_1st_vs_2nd_vs_3rd/DSG/sample_size1.txt  SAM
setwd("/camhpc/ngs/projects/TST12086/dnanexus/")


setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/2022-07-28-ZG_Edge/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/Splicing-Detection-3vs3-10-07/")
setwd("/mnt/depts/dept04/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/full_overlap/") #DSG
#setwd("/mnt/depts/dept04/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/total_back_11-10/back-2022-11-4-overlap-iPSC-NGN2/6152/") #DSG
setwd("/mnt/depts/dept04/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/DEG")           #DEG


setwd("/camhpc/ngs/projects/TST11872/dnanexus/Suyun_Branaplam/DSG") # Shuyun
#   "/edgehpc/dept/compbio/projects/TST11955/dnanexus/Suyun_Branaplam/code/Email_Dann.txt
#   
setwd("/edgehpc/dept/compbio/projects/TST11955/dnanexus/Suyun_Branaplam/Jessica_code/")


setwd("/edgehpc/dept/compbio/projects/TST11451/") # iPSC-MN from Dann
setwd("/edgehpc/dept/compbio/projects/TST11354/") # fibrobalst from Dann
#setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/2022-03-02-ZG_verified/Explore")

setwd("/camhpc/ngs/projects/TST11000/dnanexus/20210709173527_fergal.casey")
setwd("/camhpc/ngs/projects/TST11797/dnanexus/20210809135550_maria.zavodszky/EA20210812_0")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/EA20230123_0")

setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/EA20230123_0") # RNAseq MultiOmics
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/code_downstream/") # Fbinary_numDSG_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.png
#https://rstudio-prod.hpc.biogen.com/s/e3dc3cfc78a31963b422b/file_show?path=%2Fcamhpc%2Fngs%2Fprojects%2FTST11955%2Fdnanexus%2F20220723054602_Zhen.Gao%2Fcode%2FAll_code_3vs3%2Fcode_downstream%2Fbinary_numDSG_dPSIrm0.3_dPSIlc0.25_padj0.05_EdPSI0.3_PdPSI0.9_M50.png
setwd("/edgehpc/dept/compbio/users/dhuh/SMASM/20181210_TST11354_sma_fibroblast_compound_treated/majiq/res.03.deltapsi_voila_thrdpsi_0.3/")

setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/heatmap_for_catherine_excel/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_splicing/TST11955/All_1st_vs_2nd/overlap_1st_vs_2nd/full_overlap")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/NO_threshold/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/")

########## ========== OLD from CAMHPC
setwd("/camhpc/home/zgao1/Coding_system_tool_utilities/conda/")
setwd("/camhpc/home/zgao1")
setwd("/camhpc/home/zgao1/Working_code")
setwd("~/camhpc_zgao1/Working_code")

# Edge
setwd("~/CompBio_zgao1_LK/project_RNAseq/SMA_RASL-Seq/standard_code/")
setwd("/edgehpc/dept/compbio/users/zgao1/project_RNAseq/SMA_RASL-Seq/standard_code/")
setwd("/edgehpc/dept/compbio_old/projects/splice_pipeline_BSSI/data/annotation")

# run make-master_table in CAMHPC
setwd("/camhpc/ngs/projects/TST11955/dnanexus/SMA_RASL-Seq/standard_code/")

setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_1st_vs_2nd/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_3vs3_overlap_4vs8//")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/BACK-contain_R4_bad_samples/Re_Analyisis-4vs8-10-12/Splicing-Detection-4vs8-10-12/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/BACK-contain_R4_bad_samples/Re_Analyisis-4vs8-10-12/code-downstram-10-12/analysis.02.splicing_offtargets/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/BACK-contain_R4_bad_samples/analysis.03.splice_analysis/")

out_dir="/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/Result/S4088_4090_4v8/iPSC_BIO-4088-10x/" # (1)

setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/BACK-contain_R4_bad_samples/code-downstram-10-12/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/All_code_3vs3/analysis.02.splicing_offtargets_default_threshold/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/code/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/20220723054602_Zhen.Gao/Result/")
setwd("/camhpc/ngs/projects/TST11955/dnanexus/sample_sheet/")


setwd("/camhpc/ngs/projects/TST11872")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/2022-03-02-ZG_verified/Explore")
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/") # 2022-04-05 ZG√
setwd("/camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/analysis.02.splicing_offtargets/res.02.DSG_counts/") # 2022-04-05 ZG√


setwd("/mnt/depts/dept04/compbio/users/zgao1/")#setwd.R")
setwd("/mnt/depts/dept04/compbio/users/zgao1/projects_link/TST11955")#setwd.R")
setwd("/mnt/depts/dept04/compbio/users/zgao1/projects_link/TST11872")#setwd.R")
setwd("/edgehpc/dept/compbio/users/zgao1/")#setwd.R")
setwd("/edgehpc/dept/compbio/users/zgao1/coding_system_tool_utilities/Environment_and_Setting_code/")#setwd.R")
# ln -s /camhpc/ngs/projects/ projects # no mkdir needed for projects
 

setwd("/edgehpc/dept/compbio/users/zgao1")
setwd("/edgehpc/dept/compbio/users/SPLICE/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs") # =========== Not permission


########## ========== Dann # 1) Fibroblast: this one does not have the “master table”, but separate results are at:
setwd("/edgehpc/dept/compbio/users/dhuh/SMASM/20181210_TST11354_sma_fibroblast_compound_treated/")
  majiq/res.03.deltapsi_voila_thrdpsi_0.3/
  rmats/
  Leafcutter/results_m1_p0.001/
Compound names:
  RG7916 = Risdiplam
  LMI070 = Branaplam
  907272 = analogous compound from a chemical library
Also, DSG lists I complied can be found at (if it helps):
  /edgehpc/dept/compbio/users/dhuh/SMASM/20181210_TST11354_sma_fibroblast_compound_treated/analysis_02_compare_combine_lc_rm_mj/DSG_list.xlsx   ,     tab “DSG_per_dose_compound”
########## ========== Dann # 2) iPSC-MN:
setwd("/edgehpc/dept/compbio/users/dhuh/SMASM/20190913_TST11451_sma_iPSCMN_compound_treated/master_table/")
  res_comb_GM24468D-RG7961*/master_all.txt    # “GM..” is the name of the cell line used.
  res_comb_HNDS0030_01-RG7961*/master_all.txt
  res_comb_HNDS006_01C2-RG7961*/master_all.txt
########## ========== Dann # (0) NGN2
setwd("/camhpc/home/dhuh/project_RNAseq/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs/")
#For the codes that Dan ran, you can refer to 
setwd("/home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing/")
#      /home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/run.00.download_from_DNAnexus_mehools_code.qsub

# LOCAL
dir.create("/Users/zgao1/code/ggsashimi/examples/") # local MAC mounted files =========================================
dir.create("/Volumes/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/bam") # local mounted files =========================================
"/Volumes/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/bam/TST11872_iPSC-BIO-1755497-1-12nM.genome.sorted.bam"

"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|-*_master_table.csv" # this is where results from all three algorithms are combined.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- leafcutter"   #  leafcutter results. Data from “leafcutter_ds_res_cluster_significance.txt” and “leafcutter_ds_res_effect_sizes.txt” are used for the master table.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- majiq_v1"   # majiq result. Data from “cutoff_*_reference_alternative.deltapsi.tsv” are used for the master table.
"/camhpc/dept/compbio/project/splice_pipeline_BSSI/results/dev/|-BIO-*/|- rmats"        # rmats result. Data from *.MATS.JCEC.txt” are used for the master table.

###      =========> JoySplicing -> /home/jpoulo/SPLICE/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs

# Splice analysis example codes:
#   For (1), please refer to 
# “/camhpc/dept/compbio/project/splice_pipeline_BSSI/programs/dev/example_codes/” # , but I haven’t checked it myself.
# 
# For the codes that I ran, you can refer to 
# “/home/dhuh/project_RNAseq/TST11742_RNAseq_profiling_of_Skyhawk_SMA_SM_offtarget/analysis.02.splicing/”,
# 
# From “run.01.*” to  “run.06.*”
# 
# https://platform.dnanexus.com/projects/FPXy9Pj0G3b983Zv4x939bbz/monitor/analysis/G7yfpKQ0G3b9XFVY2KXY9Xfy
# https://platform.dnanexus.com/projects/FPXy9Pj0G3b983Zv4x939bbz/data/analyses/TST11872
# 
# 
# Hi team,
# I just found out that the cambridge-hpc is mounted on edge cluster. So we can access all the files from Cambridge-hpc from edge easily, although the connection is a slow one.
# See below
# I ssh into edge cluster:
# snegi1@BMD-C02ZC3U4LVDQ ~ % ssh snegi1@edge.hpc.biogen.com
# Then I can access our ngs projects using the regular path:
# 
# [snegi1@edge-hpc-log-101 ~]$ cd /camhpc/ngs/projects/TST11887

# total 305168
# drwxr-xr-x. 7 zgao1 zgao1         178 May 27 23:20 Coding_system_tool_utilities
# drwxr-xr-x. 2 zgao1 zgao1          72 Feb  8 13:27 Documents
# drwxr-xr-x. 2 zgao1 1000515        66 Feb 22 19:51 Genome_Ref
# drwxr-xr-x. 3 zgao1 zgao1          68 Feb  8 16:15 jupyter-notebook-dir
# -rw-r--r--. 1 zgao1 zgao1    76422664 Apr  3 21:15 leafviz_test1.RData
# -rw-r--r--. 1 zgao1 zgao1   123372999 Apr  3 22:11 leafviz_test2.RData
# -rw-r--r--. 1 zgao1 zgao1    76499696 Apr  3 23:28 leafviz_test3.RData
# drwxr-xr-x. 5 zgao1 zgao1        4853 Mar  3 15:44 log_master_table
# drwxr-xr-x. 2 zgao1 zgao1        1910 Apr  2 01:32 log-splicing-final-worked-records
# drwxr-xr-x. 2 zgao1 1000515        55 Mar 15 15:26 NGS_projects
# drwxr-xr-x. 2 zgao1 zgao1           0 Mar  2 15:56 output_dir
# lrwxrwxrwx. 1 zgao1 zgao1          67 Mar 29 01:15 path -> /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/path
# drwxr-xr-x. 2 zgao1 zgao1          81 May 30 11:23 Python_Code
# drwxr-xr-x. 3 zgao1 zgao1          45 Feb  7 21:43 R
# drwxr-xr-x. 3 zgao1 zgao1          45 Jan 24 19:10 R_Code
# drwxr-xr-x. 5 zgao1 zgao1         701 Jun 24 11:46 Working_code
# drwxr-xr-x. 5 zgao1 zgao1         559 Jun  4 16:37 working_record
# drwxrwxr-x. 3 zgao1 zgao1         408 Mar 29 01:18 zgao1@edge.hpc.biogen.com:CompBio_zgao1_LK
# lrwxrwxrwx. 1 zgao1 zgao1         124 Mar 29 00:13 ZG_AS_code -> /camhpc/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/code/example_TST11742_HNDS0030-RGX71083-10x/ZG_spelicng_code/
#   
# cd 
#   lrwxrwxrwx. 1 zgao1 zgao1   80 Mar 15 15:26 JoySplicing -> /home/jpoulo/SPLICE/TST11726_RNA_profiling_of_NGN2_treated_with_RNA_targeted_SMs
# lrwxrwxrwx. 1 zgao1 1000515 30 Feb 21 16:21 TST11872 -> /camhpc/ngs/projects/TST11872/

# Tigran Babayan (EXTERNAL)
# Dera Colleagues, 
# HPC team is ready to run the next step as part of COMPBIO DATA migration to EDGE, which should enable a new EDGE directory to be used:
#   All files from /edgehpc/dept/compbio will be copied to the read-only folder /edgehpc/dept/compbio_old
# Except Active projects:
#   /edgehpc/dept/compbio/ngs/TST12039  (will be moved to /edgehpc/dept/compbio/projects/TST12039)
# /edgehpc/dept/compbio/projects/GSE154659_Renthal_MouseDRG
# /edgehpc/dept/compbio/projects/TST11964_analysis
# /edgehpc/dept/compbio/projects/TST11882_igv
# /edgehpc/dept/compbio/projects/TST11978
# /edgehpc/dept/compbio/projects/AD_early_sig
# 
# The new directories structure will be implemented in /edgehpc/dept/compbio
# /edgehpc/dept/compbio
# /users/<username>
#   /data 
# /projects
# /reference
# /scripts
# /scratch
# /public
# Changes will be executed in two runs:
#   Non-active folders will be moved on 27th Sep
# Active project folders will be fixed on 28th Sep from 4 am to 8 am, when it could be temporary unavailable  
# If you have any concern regarding the proposed plan please reach out to HPC team so we can address that.
# 
# P.S Please take into account that /edgehpc/dept/compbio_old will be read only, but if you need something from that folder you will be able to copy that to a new directory structure or request HPC team to move necessary data to the proper directory.
# P.P.S. When it is covered we will run active projects migration from CAMHPC to Edge, but that will be announced separately!
#   
#   
#   Thanks for your partition,
# Tigran

# Thoma Carlile /edgehpc/dept/compbio/ngs/ONT/?
Tigran Babayan (EXTERNAL)
By default, it will be moved to 

/edgehpc/dept/compbio_old/ngs/ONT/ and will be read-only.

But if it is an active project we could move that to 

/edgehpc/dept/compbio/projects/ngs/ONT/ (and it will mean some downtime from 4 am to 8 am 28th Sep)

Hi Baohong,
Sorry I missed this email. That path works for me. @Andrei could someone on your team move the data?
Thanks,
Thomas
Baohong Zhang
Thomas,Will df -h

/edgehpc/dept/compbio/instruments         work?
  

Regards,

Tigran


/mnt/depts/dept04/compbio = edgehpc/dept/compbio/projects
/mnt/depts/dept04/compbio zgao1@edge-hpc-log-102 $ 
  $  tree -L 1
.
├── data
├── edge_condaEnv
├── edge_tools
├── genomes
├── gtxai
├── gtxau
├── human_genetics
├── jc_bfd
├── Omicsoft
├── projects
├── public
├── reference
├── scratch -> /scratch/compbio
├── scripts
└── users

/edgehpc/dept/compbio_old/projects/splice_pipeline_BSSI/data/annotation


df -h | sort -n
/mnt/depts/dept04/compbio zgao1@edge-hpc-cpu-109 $ 
camhpcisixnfs.biogen.com:/home                         910T  710T  201T  78% /camhpc/home
camhpcisixnfs.biogen.com:/ifs/HPC/CAMHPCISI/scratch    190T  4.5T  186T   3% /camhpc/scratch
camhpcisixnfs.biogen.com:/instruments                  145T   75T   71T  52% /camhpc/instruments
c
camhpcisixnfs.biogen.com:/ngs                          500T  464T   37T  93% /camhpc/ngs
camhpcisixnfs.biogen.com:/pkg                          8.0T  4.2T  3.9T  53% /camhpc/pkg
camhpcisixnfs.biogen.com:/project                      3.6P  3.0P  544T  85% /camhpc/projectcd 

camhpcisixnfs.biogen.com:/proteomics                   170T  154T   17T  91% /camhpc/proteomics
camhpcisixnfs.biogen.com:/pubdata                       65T   63T  3.0T  96% /camhpc/pubdata
camhpcisixnfs.biogen.com:/sbgrid                       3.6P  3.0P  544T  85% /programs

dept04                                                 1.9P  1.7P  193T  90% /mnt/depts/dept04



$  df -h | sort -n /edgehpc zgao1@edge-hpc-cpu-109 $ 
admin                                                  1.0T   14G 1011G   2% /edgehpc/admin
apps                                                   3.0T  2.2T  916G  71% /edgehpc/apps

dept03                                                 470T  435T   36T  93% /mnt/depts/dept03
dept04                                                 1.9P  1.7P  193T  90% /mnt/depts/dept04

hpc-test                                               1.0T   14M  1.0T   1% /edgehpc/hpctest
modulefiles                                             50G   14M   50G   1% /edgehpc/modulefiles
pubdata                                                 15T  5.7T  9.4T  38% /edgehpc/pubdata
scheduler                                               50G  8.0K   50G   1% /edgehpc/scheduler
scratch                                                100T   83T   18T  83% /scratch


$PATH
-bash: /edgehpc/dept/compbio/edge_tools/RNAsequest:/cm/shared/apps/slurm/current/sbin:/cm/shared/apps/slurm/current/bin:/edgehpc/apps/gb/spack/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:
  /home/zgao1/.local/bin:
  /home/zgao1/bin: 
/edgehpc/dept/compbio/edge_tools/RNAsequest
/edgehpc/dept/compbio/edge_tools/RNAsequest
