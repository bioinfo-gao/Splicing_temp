currdir=$(pwd)

for din in $currdir/TST*; do
	if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist, 
       		 bash $din/argparse_wrapper_qsub_leafcutter.bash # then run the leafcutter job
	fi
done
