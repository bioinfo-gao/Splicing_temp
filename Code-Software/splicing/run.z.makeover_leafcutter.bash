for din in TST*; do
	if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist,
	  sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_leafcutter.bash
        fi
done

#
#for din in TST*excrep4; do
#        if [ ! -d ${din}/leafcutter ]; then # if leafcutter directory does not exist,
#          sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_leafcutter.bash
#		echo 
#        fi
#done

