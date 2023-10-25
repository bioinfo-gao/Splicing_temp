currdir=$(pwd)

for din in $currdir/TST*; do
	bash $din/optparse_qsub_make_master_table.sh
done
