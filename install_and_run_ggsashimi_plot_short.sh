###################### ========= MAC ============= ########################### https://github.com/guigolab/ggsashimi/issues/25. chr4:3212709-3213958
# copy from server /camhpc/homez/gao1/Working_code/ggsashimi/
# run in local mac 
# ../ggsashimi.py -b $input_single_bam   -c chr4:3211709-3214958 -g $anno_file -M 2 -C 3 -O 3 --shrink --alpha 0.25 --base-size=16 --ann-height=3 --height=3 --width=18 -P palette3.txt -o HTT
# ../ggsashimi.py -b HTT_bams_compound1.tsv   -c chr4:3211709-3214240 -g $anno_file -M 2 -C 3 -O 3 --shrink --alpha 0.25 --base-size=16 --ann-height=3 --height=3 --width=18 -P palette3.txt -o HTT_compound1

anno_file="/Users/zgao1/Library/CloudStorage/OneDrive-Biogen/Code/ggsashimi/examples-TST11872_Mac2023/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf"

../ggsashimi.py -b HTT_bams_compound111.txt   -c chr4:3211709-3214200 -g $anno_file -M 2 -C 3 -O 3 --shrink --alpha 0.2 --base-size=16 --ann-height=1 --height=3 --width=18 -P palette3.txt -o HTT_compound1


# ## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
# #./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file     -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt # 1
# ../ggsashimi.py -b input_bams.tsv    -c chr10:27040584-27048100 -g annotation.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt
# ../ggsashimi.py -b input_bams.tsv    -c chr10:27040584-27048100                   -M 10 -C 3 -O 3          --alpha 1    --base-size=16                --height=3 --width=18 --fix-y-scale -F tiff -R 350  -A mean # Mark mean only
# ## Example #2. Mean coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors, fixed y-scale
#   -b BAM,          --bam           BAM            Individual bam file or file with a list of bam files. In the case of a list of files the format is tsv:1col: id for bam file, 2col: path of bam file, 3+col: additional columns
#   -c COORDINATES,  --coordinates   COORDINATES    Genomic region. Format: chr:start-end.   Remember that bam coordinates are 0-based
#   -g GTF,          --gtf           GTF            Gtf file with annotation (only exons is enough)
# # -o OUT_PREFIX,   --out-prefix OUT_PREFIX                   Prefix for plot file name [default=sashimi]
#   -M MIN_COVERAGE, --min-coverage  MIN_COVERAGE   Minimum number of reads supporting a junction to be drawn [default=1]
#   -C COLOR_FACTOR, --color-factor  COLOR_FACTOR   Index of column with color   levels (1-based)
#   -O OVERLAY,      --overlay       OVERLAY        Index of column with overlay levels (1-based)
#                    --shrink                       Shrink the junctions by a factor for nicer display         
#                    --alpha         ALPHA          Transparency level for density histogram [default=0.5]
#                    --base-size     BASE_SIZE      Base font size of the plot in pch [default=14]
#                    --ann-height    ANN_HEIGHT     Height of annotation plot in inches [default=1.5]
#                    --height        HEIGHT         Height of the individual signal plot in inches [default=2]
#                    --width         WIDTH          Width of the plot in inches [default=10]
#   -P PALETTE,      --palette       PALETTE        Color palette file. tsv file with >=1 columns, where the color is the first column. Both R color names and hexadecimal values are valid
#                    --fix-y-scale                  Fix y-scale across individual signal plots                    [default=False]
#   -F OUT_FORMAT,   --out-format    OUT_FORMAT     Output file format: <pdf> <svg> <png> <jpeg> <tiff>             [default=pdf]
#   -A AGGR,         --aggr          AGGR           Aggregate function for overlay: <mean> <median>  <mean_j> <median_j>. Use mean_j | median_j to keep  density overlay but aggregate junction counts   [default=no aggregation]
#   
#   
#   
#   -L LABELS, --labels LABELS Index of column with labels (1-based) [default=1]
#   -R OUT_RESOLUTION, --out-resolution OUT_RESOLUTION    Output file resolution in PPI (pixels per inch). Applies only to raster output formats [default=300]
#   

###################### ========= MAC ============= ########################### https://github.com/guigolab/ggsashimi/issues/25
# succss record in local MAC
mkdir -p ~/code/ # == /User/zgao1/code/

git clone https://github.com/guigolab/ggsashimi.git
conda install pysam
zgao1@ zhens-mbp $  conda activate py27 # changed in 2023
# (py2) ~ :  
zgao1@ zhens-mbp $  conda install -c bioconda pysam # changed 2023
(py2) 
cd ~/code/ggsashimi :  
chmod a+x ggsashimi.py

../sashimi-plot.py -b input_bams.tsv -c chr10:27040584-27048100 -g annotation.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt

input_single_bam="/Volumes/ngs/projects/TST11872/dnanexus/20220204181515_zhen.gao/bam/TST11872_iPSC-BIO-1755497-1-12nM.genome.sorted.bam"
anno_file="/Volumes/dept/compbio/project/splice_pipeline_BSSI/data/annotation/Human.GRCh38.v34.l1_5.ERCC/Human.GRCh38.v34.l1_5.ERCC.transcript.gtf"
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4   --height=3  --width=18 -P palette.txt # 1
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4   --height=18 --width=18 -P palette.txt # 2
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=18  --height=5  --width=18 -P palette.txt #3
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=14 --ann-height=1.5 --height=2  --width=10 -P palette.txt #4
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=14 --ann-height=1.5 --height=2  --width=28 -P palette.txt #5
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25 -P palette.txt #6
./ggsashimi.py -b $input_single_bam -c chr10:27040584-27048100 -g $anno_file -M 10 -C 3 -O 3 --shrink --alpha 0.25  -P palette.txt -F svg # 7

Error in loadNamespace(name) : there is no package called ‘svglite’
Calls: ggsave ... tryCatch -> tryCatchList -> tryCatchOne -> <Anonymous>
Execution halted

R install.packages("svglite", dependencies=T)
ERROR: dependency ‘systemfonts’ is not available for package ‘svglite’
* removing ‘/opt/anaconda3/envs/py2/lib/R/library/svglite’

conda install -c conda-forge r-svglite

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Individual bam file or file with a list of bam files. In the case of a list of files the format is tsv:1col: id for bam file, 2col: path of bam file, 3+col: additional columns
  -c COORDINATES, --coordinates COORDINATES     Genomic region. Format: chr:start-end. Remember that            bam coordinates are 0-based
  -o OUT_PREFIX, --out-prefix OUT_PREFIX        Prefix for plot file name [default=sashimi]
  -S OUT_STRAND, --out-strand OUT_STRAND        Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=both]
* -M MIN_COVERAGE, --min-coverage MIN_COVERAGE  Minimum number of reads supporting a junction to be drawn [default=1]
  -j JUNCTIONS_BED, --junctions-bed JUNCTIONS_BED
                        Junction BED file name [default=no junction file]
* -g GTF, --gtf GTF     Gtf file with annotation (only exons is enough)
  -s STRAND, --strand STRAND  Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=NONE]
  --shrink              Shrink the junctions by a factor for nicer display                               [default=False]
* -O OVERLAY, --overlay OVERLAY
                        Index of column with overlay levels (1-based)
  -A AGGR, --aggr AGGR  Aggregate function for overlay: <mean> <median>  <mean_j> <median_j>. Use mean_j | median_j to keep
                        density overlay but aggregate junction counts                         [default=no aggregation]
  -C COLOR_FACTOR, --color-factor COLOR_FACTOR
                        Index of column with color levels (1-based)
  --alpha ALPHA         Transparency level for density histogram [default=0.5]
  -P PALETTE, --palette PALETTE
                        Color palette file. tsv file with >=1 columns, where
                        the color is the first column. Both R color names and
                        hexadecimal values are valid
  -L LABELS, --labels LABELS
                        Index of column with labels (1-based) [default=1]
  --fix-y-scale         Fix y-scale across individual signal plots                    [default=False]
  --height HEIGHT           Height of the individual signal plot in inches [default=2]
  --ann-height ANN_HEIGHT   Height of annotation plot in inches [default=1.5]
  --width WIDTH             Width of the plot in inches [default=10]
  --base-size BASE_SIZE     Base font size of the plot in pch [default=14]
  -F OUT_FORMAT, --out-format OUT_FORMAT        Output file format: <pdf> <svg> <png> <jpeg> <tiff>             [default=pdf]
  -R OUT_RESOLUTION, --out-resolution OUT_RESOLUTION    Output file resolution in PPI (pixels per inch). Applies only to raster output formats [default=300]
  --debug-info          Show several system information useful for debugging purposes [default=None]
  --version             show program's version number and exit'



## ====================
##  ggsashimi examples
## ====================

## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
../ggsashimi.py -b input_bams.tsv -c chr10:27040584-27048100 -g annotation.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt

## Example #2. Mean coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors, fixed y-scale
../ggsashimi.py -b input_bams.tsv -c chr10:27040584-27048100 -M 10 -C 3 -O 3 -A mean --alpha 1 -F tiff -R 350 --base-size=16 --height=3 --width=18 --fix-y-scale

