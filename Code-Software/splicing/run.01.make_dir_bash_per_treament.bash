
template=./template_script/*.bash

while IFS="" read -r line; do
	lin=($line)
       	alt=${lin[0]} # alt prefix
       	ref=${lin[1]} # ref prefix
	
	temp1=${alt//X-7/X7}
	temp2=${temp1//C-1/C1}
	doutname=${temp2/S-0/S0} #output directory names

	mkdir -p $doutname # make the output directory per treatment
	
	cp $template $doutname # copy the template scripts
	sed -i "s/OUTDIR/$doutname/g" $doutname/*.bash
	sed -i "s/ALTBAMFILE/$alt/g" $doutname/*.bash # modify the scripts per treatment
	sed -i "s/REFBAMFILE/$ref/g" $doutname/*.bash

done < ./config/config.txt
