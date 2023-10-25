for din in TST*; do
	dout=${din}_excrep4
	mkdir $dout

	cp $din/*.bash $dout

	sed -i 's/.bam)/.bam | grep -v rep4 )/g' $dout/*.bash
	sed -i "s/$din/$dout/g" $dout/*.bash
done

