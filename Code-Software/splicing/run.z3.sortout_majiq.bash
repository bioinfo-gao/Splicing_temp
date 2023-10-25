for din in TST*; do
	if compgen -G "$din/majiq_v1/cutoff_0.1*" > /dev/null; then
#            ls $din/majiq_v1/cutoff_0.1*
		echo "ex"
        else 
            rm -r $din/majiq_v1/
#		echo $din	
        fi
done
