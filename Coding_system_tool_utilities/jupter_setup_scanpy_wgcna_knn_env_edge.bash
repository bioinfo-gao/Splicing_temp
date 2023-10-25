module load anaconda3/4.9.2
conda activate /edgehpc/dept/compbio/users/zgao1/conda_gao/r_main
jupyter lab --no-browser --ip=`/usr/bin/hostname -I | /usr/bin/grep -oE "10.106.192.[0-9]{1,3}"`
