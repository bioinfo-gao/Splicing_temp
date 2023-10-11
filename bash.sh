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
