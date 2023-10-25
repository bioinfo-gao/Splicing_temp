for din in TST*; do
	if [ ! -d ${din}/majiq_v1 ]; then # if leafcutter directory does not exist,
	  #sed -i "s/cpu.q/cpu.q@ginseng02/g" $din/argparse_wrapper_qsub_majiq.bash
	  sed -i "s/cpu.q@ginseng03/cpu.q@ginseng04/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/cpu.q@ginseng01/cpu.q/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/-V cpu.q@ginseng20/cpu.q/g" $din/argparse_wrapper_qsub_majiq.bash
	  #sed -i "s/thread 12/thread 12 -V/g"  $din/argparse_wrapper_qsub_majiq.bash
        fi
done


