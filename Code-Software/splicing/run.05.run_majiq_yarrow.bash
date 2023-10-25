currdir=$(pwd)

for din in $currdir/TST*; do
	if [ ! -d ${din}/majiq_v1 ]; then
        	bash $din/argparse_wrapper_qsub_majiq.bash
	fi		
done
