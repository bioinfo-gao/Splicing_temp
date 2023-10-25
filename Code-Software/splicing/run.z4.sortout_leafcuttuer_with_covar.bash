currdir=$(pwd)

:<<"end"
for din in $currdir/TST*; do
	if [ -d ${din}/leafcutter ]; then # if leafcutter directory does not exist, 
		newname=${din}/leafcutter_prev_nocovar
		mv ${din}/leafcutter $newname
	fi
done


for din in ./take03/TST*/leafcutter; do
	lin=(${din//\// })
	dout=${lin[2]}
	cp -r $din $dout
done

end


for fin in $currdir/TST*/*master_table.csv; do
	fout=${fin//.csv/_nocovar.csv}
	mv $fin $fout
done

