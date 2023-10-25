currdir=$(pwd)
for din in $currdir/TST*excrep4; do
	echo $din
        bash $din/argparse_wrapper_qsub_rmats.bash
done
