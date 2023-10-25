
import argparse
import datetime
import os
import sys
import time
import random
import re 
import subprocess
import pandas as pd

parser = argparse.ArgumentParser(description = "Splicing analysis pipeline in argparse format.")

#Make bam file list arguments
parser.add_argument('--alt_bam_file', default=[], nargs='+', help='a list full file paths of indexed, sorted bam files to undergo splicing analysis in the alternative group.  Separate each bam file using a space after calling the option.')
parser.add_argument('--ref_bam_file', default=[], nargs='+', help='a list of full file paths of indexed, sorted bam files to undergo splicing analysis in the reference group. Separate each bam file using a space after calling the option.')

#Leave out sequencing options
# parser.add_argument('STARIndexFolder',help='The full path to the STAR index folder prepared using the alignment settings for your alignment run.')
# parser.add_argument('genome_file',help='The full path to the reference genome FASTA file to align to.')

#Add store true options to run each tool separately
parser.add_argument('--rmats', help="Include this flag to run the rMATs analysis (v4.0.2)", default=False, action='store_true')
parser.add_argument('--leafcutter', help="Include this flag to run the leafcutter analysis (uses regtools for bam2junc)", default=False, action='store_true')
parser.add_argument('--majiq_v1', help="Include this flag to run the majiq analysis (v1.1.4)", default=False, action='store_true')
parser.add_argument('--majiq_v2', help="Include this flag to run the majiq analysis (v2) (Note: not available)", default=False, action='store_true')

#We will need these options
parser.add_argument('--read_type',help="The sequencing type. Choose either 'single' or 'paired'. [default: paired]", default='paired')
parser.add_argument('--read_length',help='The sequencing read length from the sequencing run.',type=int)
parser.add_argument('--anno_file', help='The full path to the gene structure annotation gtf file.')
parser.add_argument('--gff3_file',help='The full path to the general feature format file.')
parser.add_argument('--leafcutter_anno_file',help='The full path to the leafcutter annotation file.')
parser.add_argument('--leafcutter_strandedness',help='Strand specificity of library.  0 = unstranded, 1 = first-strand/RF, 2 = second-strand/FR.' )
parser.add_argument('--genome_file',help='The full path to the reference genome FASTA file to align to (used to create majiq config file).') #Needed for majiq v1, not needed for majiq v2

#Add extra arguments for rMATS
#None with v4.0.2

#Add extra arguments for leafcutter
#Add group file argument - add as full file path, call it 'covariate_file'. Define the format.  
# First column is (sorted) bam file name (just the name, no path needed) and will need to have a column name 'bamfiles'
# Second columns is ref/alt
# Function will read this in, subset it to a file containing the bam file names and a column of ref/alt assignments
parser.add_argument('--leafcutter_covar', help="Include this flag to include a leafcutter covariate file", default=True, action='store_true') #DHmod from default=False
parser.add_argument('--leafcutter_covar_file', help="The full path for a whitespace-delimited, subject-level covariate metadata file used for leafcutter.  Needs at least a column called 'bamfiles' with the bam file names from the experiment.  Can add subsequent column(s) of covariates.")

#Covariate file will not be used for this project.  The function to generate a dummy covariate file was added and then we could turn the covariate file back on at another point in time

#Add extra arguments for majiq_v1
parser.add_argument('--majiq_cutoffs', default=[], nargs='+', help='the dPSI cutoffs to test.  Separate each with a space after calling the option.')

#Make an out directory argument
parser.add_argument('--out_dir', help='The full path for the output files.')

#Set up the parsed arguments
args = parser.parse_args()

###### Run rmats #####
if args.rmats:
  #starting with bam files, rmats requires two files called "b1.txt" and "b2.txt" config files
  #They are comma-separated lists of bam files (full paths)
  def rmatsFile(out_dir,ref_bamfiles, alt_bamfiles):
    #Set the output directory
    os.chdir(out_dir)
    
    #make a directory to hold rmats files
    if not os.path.exists("rmats"):
      os.makedirs("rmats")
    
    #make the string needed by rmats
    ref_bamstring = ','.join(ref_bamfiles)
    alt_bamstring = ','.join(alt_bamfiles)
    
    #Seems to need to create a new file
    b1_file = open(os.path.join(out_dir, "rmats", "b1.txt"), 'w')
    b2_file = open(os.path.join(out_dir, "rmats", "b2.txt"), 'w')
    
    #Now write the strings out
    b1_file.write(ref_bamstring)
    b1_file.close()
    b2_file.write(alt_bamstring)
    b2_file.close()
    
  #Run the actual command to make those files out
  rmatsFile(args.out_dir, args.ref_bam_file, args.alt_bam_file)
  
  #Run the rMATS analysis
  rmats = subprocess.run("source /etc/profile.d/modules_bash.sh; module purge; module add rmats/4.0.2;  module load python/2.7.11; python /camhpc/pkg/rmats/4.0.2/centos6/rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --nthread 8".format(os.path.join(args.out_dir, "rmats", "b1.txt"), os.path.join(args.out_dir, "rmats", "b2.txt"), args.anno_file, os.path.join(args.out_dir, "rmats"), args.read_type, args.read_length), shell=True)

###### Run leafcutter #######
if args.leafcutter:
  #Placeholder function in case formatting is needed
  def leafcutterFile(out_dir):
    os.chdir(out_dir)
    if not os.path.exists("leafcutter"):
      os.makedirs("leafcutter")
  
  leafcutterFile(args.out_dir)

  #bam2junc 
  #Assume already run samtools index on bam files
  #Assume RF strandedness for stranded library prep
  juncfilelist = "juncfiles.txt"
  
  for bamfile in args.ref_bam_file+args.alt_bam_file:
    #Get the file name
    bam_filename = os.path.basename(bamfile)
    juncname = [bam_filename, ".junc"]
    #Need to set -s as required 0 = unstranded, 1 = first-strand/RF, 2 = second-strand/FR
    bam2junc = subprocess.run("/home/lhou1/bin/regtools junctions extract -a 8 -m 50 -s {} -M 500000 {} -o {}; echo {} >> {}".format(args.leafcutter_strandedness, bamfile, os.path.join(args.out_dir, "leafcutter", "".join(juncname)), os.path.join(args.out_dir, "leafcutter", "".join(juncname)), os.path.join(args.out_dir, "leafcutter", juncfilelist)), shell=True)
  
  #junc2clust
  clustname = "junc2clust" #Don't use os.path.join(args.out_dir, "leafcutter", clustname) in the -o below, just use clustname
  #Need intermediate step for removing junctions regtools annotates with a "?"
  base_ref=[os.path.basename(x) for x in args.ref_bam_file]
  base_alt=[os.path.basename(x) for x in args.alt_bam_file]
  junc_files = [x + ".junc" for x in base_ref+base_alt]
  for juncfile in junc_files:
    junc_clean = subprocess.run("grep -v '?' {} > {}; mv {} {}".format(os.path.join(args.out_dir, "leafcutter", juncfile), os.path.join(args.out_dir, "leafcutter", "tmp"), os.path.join(args.out_dir, "leafcutter", "tmp"), os.path.join(args.out_dir, "leafcutter", juncfile)), shell=True)
  #Run junc2clust after cleaning
  leafcutterpath1 = os.path.join('/home/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/scripts','leafcutter_cluster_regtools.py')
  junc2clust = subprocess.run("source /etc/profile.d/modules_bash.sh; module purge; module add anaconda; source activate /home/emarshal/.conda/envs/leafcutter; cd {}; python {} -j {} -m 50 -o {} -l 500000 --checkchrom True".format(os.path.join(args.out_dir, "leafcutter"),leafcutterpath1, os.path.join(args.out_dir, "leafcutter", juncfilelist), clustname), shell=True)
  rmjuncfilelist = subprocess.run("rm {}".format(os.path.join(args.out_dir, "leafcutter", juncfilelist)), shell=True) #Remove juncfiles.txt or it will grow as code is rerun
  
  #differential splicing (ds) - two parts - config file (groups_file.txt) and the ds script
  
  if args.leafcutter_covar:
  #Create groups_file.txt in args.out_dir from args.leafcutter_covar_file
    def create_leafcutter_groups_file(out_dir, ref_bamfiles, alt_bamfiles, covar_file):
      #Read covar file in
      covar = pd.read_csv(covar_file, delim_whitespace=True)
      #Check the column bamfiles is present
      if ~covar.columns.isin(['bamfiles']).any():
        raise ValueError('The covariate file must have a column name called bamfiles.')
      #Get the bam file basenames from ref and alt
      my_ref_bam_list = [os.path.basename(x) for x in ref_bamfiles]
      my_alt_bam_list = [os.path.basename(x) for x in alt_bamfiles]
      #Filter covar for my_bam_list
      ref_covar_filt = covar[covar['bamfiles'].isin(my_ref_bam_list)]
      alt_covar_filt = covar[covar['bamfiles'].isin(my_alt_bam_list)]
      #Make the group columns for each
      ref_covar_filt['group'] = 'reference'
      alt_covar_filt['group'] = 'alternative'
      #concat and arrange columns such that 1st column = bamfiles, 2nd column = group
      covar_format = pd.concat([ref_covar_filt, alt_covar_filt])
      my_columns = ['bamfiles', 'group']
      covar_format = covar_format[ my_columns + [ col for col in covar_format.columns if col not in my_columns ] ]
      #Write this file out, it is the group_file.txt (no rownames, no header)
      covar_format.to_csv(os.path.join(out_dir,"groups_file.txt"), sep='\t', index=False, header=False)
    #Write the groups_file
    create_leafcutter_groups_file(out_dir=os.path.join(args.out_dir,"leafcutter"), ref_bamfiles=args.ref_bam_file, alt_bamfiles=args.alt_bam_file, covar_file=args.leafcutter_covar_file)
      
  if not args.leafcutter_covar:
    #Create a dummy version of groups_file.txt
    def create_dummy_leafcutter_groups_file(out_dir, ref_bamfiles, alt_bamfiles):
      #Get the bam file basenames from ref and alt
      my_ref_bam_list = [os.path.basename(x) for x in ref_bamfiles]
      my_alt_bam_list = [os.path.basename(x) for x in alt_bamfiles]
      #Append second 'column' of groups
      ref_dummy = [my_ref + ' reference' for my_ref in my_ref_bam_list]
      alt_dummy = [my_alt + ' alternative' for my_alt in my_alt_bam_list]
      #Combine lists
      dummy_groups = ref_dummy + alt_dummy
      #Write to a file (set out_dir)
      dummy_file = open(os.path.join(out_dir, "groups_file.txt"), 'w')
      for element in dummy_groups:
        dummy_file.write(element + "\n")
      dummy_file.close()
    create_dummy_leafcutter_groups_file(out_dir=os.path.join(args.out_dir,"leafcutter"), ref_bamfiles=args.ref_bam_file, alt_bamfiles=args.alt_bam_file)
  
  #Run the differential splicing script
  # assume the mapping file, call it groups_file.txt, has been rendered before running the pipeline
  # assume the exon file has been rendered from the gtf file, this is the file path args.leafcutter_anno_file from args (use Rscript gtf_to_exons.R on the gtf file).
  # i is min. samples per intron
  # g is min. samples per group
  groupsfile = "groups_file.txt"
  dsname = "junc2clust_perind_numers.counts.gz"
  #leafcutterpath2 = os.path.join('/home/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/scripts','leafcutter_ds.R')
  #Use my leafcutter script that accounts for the regtools strandedness formatting
  leafcutterpath2 = os.path.join('/home/mryals/splicing','leafcutter_ds_mry.R')
  leafcutter_ds_outname = "leafcutter_ds_res"
  ds = subprocess.run("source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; source activate /home/emarshal/.conda/envs/leafcutter; cd {}; Rscript {} --num_threads 8 --exon_file {} -i 2 -g 2 {} {} -o {}".format(os.path.join(args.out_dir, "leafcutter"),leafcutterpath2, args.leafcutter_anno_file, os.path.join(args.out_dir, "leafcutter", dsname), os.path.join(args.out_dir, "leafcutter", groupsfile), os.path.join(args.out_dir, "leafcutter", leafcutter_ds_outname)), shell=True)
  
  #ds result visualizations
  #Makes a pdf of the differentially spliced clusters at FDR 5%
  #Requires the junc2clust, significance file, groups file, and an exon annotation file (preset here)
  #exon_file = os.path.join('/home/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/leafcutter/data','gencode19_exons.txt.gz')
  #sig_file = "leafcutter_cluster_significance.txt"
  #plots_script = os.path.join('/home/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/leafcutter/scripts','ds_plots.R')
  #ds_plots = subprocess.run("source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; source activate /home/emarshal/.conda/envs/leafcutter; cd {}; Rscript {} -e {} {} {} {} -f 0.05".format(os.path.join(args.out_dir, "leafcutter"), plots_script, exon_file, os.path.join(args.out_dir, "leafcutter", dsname), os.path.join(args.out_dir, "leafcutter", groupsfile), os.path.join(args.out_dir, "leafcutter", sig_file))
  
  #ds shiny app visualizations 
  #The gtf annotation file has to be parsed for the app.  
  #Requires the junc2clust, the significance file, the effect sizes file, and the annotation code at minimum
  #Prep results creates an Rdata file with annotations that can be used when compiling the full leafcutter result table
  #Use my leafcutter script that accounts for the regtools strandedness formatting
  #prep_script = os.path.join('/home/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/leafviz','prepare_results.R')
  prep_script = os.path.join('/home/mryals/splicing/','prepare_results_mry.R')
  prep_name = "leafviz.Rdata"
  prep_viz = subprocess.run("source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; source activate /home/emarshal/.conda/envs/leafcutter; cd {}; Rscript {} -m {} -o {} -f 1 {} {} {} {}".format(os.path.join(args.out_dir,"leafcutter"),prep_script,os.path.join(args.out_dir,"leafcutter",groupsfile),os.path.join(args.out_dir,"leafcutter",prep_name),os.path.join(args.out_dir,"leafcutter",dsname),os.path.join(args.out_dir,"leafcutter",leafcutter_ds_outname+"_cluster_significance.txt"),os.path.join(args.out_dir,"leafcutter",leafcutter_ds_outname+"_effect_sizes.txt"), os.path.join("/camhpc/dept/compbio/project/splice_pipeline_BSSI/mry_annot","leafcutterviz","leafcutterviz") ), shell=True)

###### Run majiq v1 #######
if args.majiq_v1:
  #Make an output subfolder for majiq
  def majiqFile(out_dir):
    os.chdir(out_dir)
    if not os.path.exists("majiq_v1"):
      os.makedirs("majiq_v1")
  
  majiqFile(args.out_dir)
  
  #Run the majiq_build - two steps - make a config file, config.txt and run majiq build script
  def create_majiq_config_file_v1(out_dir, ref_bamfiles, alt_bamfiles, geno_file, read_length):
  
    #Create the file 
    out_file = open(os.path.join(out_dir,'majiq_v1','config.txt'), 'w')
    #Get the bam file basenames from ref and alt and strip the .bam (required by majiq)
    my_ref_bam_list = [item.replace('.bam','') for item in [os.path.basename(x) for x in ref_bamfiles]]
    my_alt_bam_list = [item.replace('.bam','') for item in [os.path.basename(x) for x in alt_bamfiles]]
    #create the lines that go into config.txt
    ref_line = 'reference=' + ','.join(my_ref_bam_list)
    alt_line = 'alternative=' + ','.join(my_alt_bam_list)
  
    #Write out the lines of the info part of config.txt
    out_file.write('[info]\n')
    out_file.write('readlen=' + str(read_length) + '\n')
    out_file.write('samdir=' + os.path.dirname(ref_bamfiles[1]) + '\n') #This is different than when using v2
    out_file.write('genome=hg19\n') #this is based on the SMN1 genome fasta file, update to hg38 if needed
    #out_file.write('strandedness=forward\n') #This was not included in the wrapper, defaults to 'None'
    out_file.write('type=strand-specific\n') #required for v1 (with negative strand as reference) if using strand-specific RNASeq data
    out_file.write('genome_path=' + geno_file + '\n') #required for v1 QC 
    #Write out the lines of the experiments part of config.txt
    out_file.write('[experiments]\n')
    out_file.write(ref_line + '\n')
    out_file.write(alt_line + '\n')
    out_file.close()
  #Create the config file
  create_majiq_config_file_v1(out_dir = args.out_dir, ref_bamfiles = args.ref_bam_file, alt_bamfiles = args.alt_bam_file, geno_file=args.genome_file, read_length=args.read_length)
  #Run majiq_build
  #Assume gff3 file is already made from the existing gtf file, it will go to args.gff3_file
  majiq_build = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4; majiq build {} -c {} -o {} -j 8'.format(args.gff3_file, os.path.join(args.out_dir,"majiq_v1","config.txt"), os.path.join(args.out_dir,"majiq_v1")), shell=True)
  
  #Run majiq_deltapsi
  #Need the list of reference and control majiq files from majiq build
  ref_majiq1 = [item.replace('.bam','.majiq') for item in [os.path.basename(x) for x in args.ref_bam_file]]
  ref_majiq2 = [os.path.join(args.out_dir,'majiq_v1',ref) for ref in ref_majiq1]
  ref_majiq2 = ' '.join(ref_majiq2) #additional formatting step #1
  ref_majiq2 = '\"' + ref_majiq2 + '\"' #additional formatting step #2
  alt_majiq1 = [item.replace('.bam','.majiq') for item in [os.path.basename(x) for x in args.alt_bam_file]]
  alt_majiq2 = [os.path.join(args.out_dir,'majiq_v1',alt) for alt in alt_majiq1]
  alt_majiq2 = ' '.join(alt_majiq2) #additional formatting step #1
  alt_majiq2 = '\"' + alt_majiq2 + '\"' #additional formatting step #2
  majiq_deltapsi = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4;grp1={};grp2={}; majiq deltapsi -grp1 $grp1 -grp2 $grp2 -o {} --name reference alternative -j 8'.format(ref_majiq2, alt_majiq2, os.path.join(args.out_dir, "majiq_v1")), shell=True)
  
  #Run majiq_voila deltapsi - the result visualizations and tsv output are bundled into a single v1 command
  #Set the loop of jobs to run for the args.majiq_cutoffs
  #For the purposes of this pipeline, we can disable the html reports and only get the tsv file
  majiq_deltapsi_file=os.path.join(args.out_dir, "majiq_v1", "reference_alternative.deltapsi.voila")
  splicegraph_file = os.path.join(args.out_dir, "majiq_v1", "splicegraph.sql")
  for cutoff in args.majiq_cutoffs:
    voila_deltapsi = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4; voila deltapsi {} -s {} -o {} -j 8 --non-changing-threshold {} --threshold {} --disable-html --show-all'.format(majiq_deltapsi_file, splicegraph_file, os.path.join(args.out_dir,"majiq_v1"), cutoff, cutoff), shell=True)
    #The output file for a single run needs to be renamed with a cutoff in the filename
    out_file = os.path.join(args.out_dir, "majiq_v1", "reference_alternative.deltapsi.tsv")
    out_file_rename = os.path.join(args.out_dir, "majiq_v1", "cutoff_" + str(cutoff) + "_" + "reference_alternative.deltapsi.tsv")
    voila_deltapsi_rename = subprocess.run('mv {} {}'.format(out_file, out_file_rename), shell=True)
    
###### Run majiq v2 #######
if args.majiq_v2:
  #Make an output subfolder for majiq
  def majiqFile(out_dir):
    os.chdir(out_dir)
    if not os.path.exists("majiq_v2"):
      os.makedirs("majiq_v2")
  
  majiqFile(args.out_dir)
  
  #Run the majiq_build - two steps - make a config file, config.txt and run majiq build script
  #Create config file
  #Assume hg38 genomy assembly
  #Assume stranded library and that all bam files used the same strandedness
  #Assume the bam files have been sorted and indexed
  #Assume the bam files are located in a common directory
  #Assume strandedness in library preparation
  #due to the unusual formatting, this function borrows heavily from the existing function in wrapper_majiq.py
  def create_majiq_config_file_v2(out_dir, ref_bamfiles, alt_bamfiles, geno_file, read_length):
  
    #Create the file 
    out_file = open(os.path.join(out_dir,'majiq','config.txt'), 'w')
    #Get the bam file basenames from ref and alt and strip the .bam (required by majiq)
    my_ref_bam_list = [item.replace('.bam','') for item in [os.path.basename(x) for x in ref_bamfiles]]
    my_alt_bam_list = [item.replace('.bam','') for item in [os.path.basename(x) for x in alt_bamfiles]]
    #create the lines that go into config.txt
    ref_line = 'reference=' + ','.join(my_ref_bam_list)
    alt_line = 'alternative=' + ','.join(my_alt_bam_list)
  
    #Write out the lines of the info part of config.txt
    out_file.write('[info]\n')
    out_file.write('readlen=' + str(read_length) + '\n')
    out_file.write('bamdirs=' + os.path.dirname(ref_bamfiles[1]) + '\n')
    out_file.write('genome=hg38\n')
    out_file.write('strandedness=forward\n') #This was not included in the wrapper, defaults to 'None'
    #out_file.write('genome_path=' + geno_file + '\n') #Not needed for v2
    #Write out the lines of the experiments part of config.txt
    out_file.write('[experiments]\n')
    out_file.write(ref_line + '\n')
    out_file.write(alt_line + '\n')
    out_file.close()
  
  create_majiq_config_file_v2(out_dir = args.out_dir, ref_bamfiles = args.ref_bam_file, alt_bamfiles = args.alt_bam_file, geno_file=args.genome_file, read_length=args.read_length)
  
  #Run majiq_build
  #Assume gff3 file is already made from the existing gtf file, it will go to args.gff3_file
  majiq_build = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4; majiq build {} -c {} -o {} -j 8'.format(args.gff3_file, os.path.join(args.out_dir,"majiq","config.txt"), os.path.join(args.out_dir,"majiq")), shell=True)
  
  #Run majiq_deltapsi
  #Need the list of reference and control majiq files
  ref_majiq1 = [item.replace('.bam','.majiq') for item in [os.path.basename(x) for x in args.ref_bam_file]]
  ref_majiq2 = [os.path.join(args.out_dir,'majiq',ref) for ref in ref_majiq1]
  alt_majiq1 = [item.replace('.bam','.majiq') for item in [os.path.basename(x) for x in args.alt_bam_file]]
  alt_majiq2 = [os.path.join(args.out_dir,'majiq',alt) for alt in alt_majiq1]
  majiq_deltapsi = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4; majiq deltapsi -grp1 {} -grp2 {} -o {} --name reference alternative -j 8'.format(ref_majiq2, alt_majiq2, os.path.join(args.out_dir, "majiq")), shell=True)
  
  #Run majiq_voila deltapsi - maybe skip this until we know how we're planning to present the results
  
  #Run majiq_voila to create a human-readable tsv file
  deltapsi_file="reference_alternative.deltapsi.voila"
  majiq_voila_tsv = subprocess.run('source /etc/profile.d/modules_bash.sh; module purge; module load anaconda; unset PYTHONPATH; source activate /home/emarshal/.conda/envs/majiq-1.1.4; voila tsv {} {} -f {}'.format(os.path.join(args.out_dir,"majiq","splicegraph.sql"), os.path.join(args.out_dir, "majiq", deltapsi_file),  os.path.join(args.out_dir, "majiq", deltapsi_file+'.tsv') ), shell=True)


#Here can either pipe in an optparse result summary code or can stop here and then run the optparse script separately
