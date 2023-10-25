template=./template_script/optparse_qsub_make_master_table.sh
	
currdir=$(pwd)
for din in $currdir/TST*; do
	doutname=${din##*/}
	cp $template $din
	sed -i "s/INDIR/$doutname/g" $doutname/optparse_qsub_make_master_table.sh
done

