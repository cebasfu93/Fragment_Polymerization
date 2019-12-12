#!/bin/bash

NM="scripts/NanoModeler/NanoModeler.py"

for l in 16
do
	echo Starting with L${l}
	mkdir -p MIDFILES/L${l}
	mkdir -p OUTPUT/L${l}
	MID="MIDFILES/L${l}"
	OUT="OUTPUT/L${l}"
	
	echo Polymerizing L${l}
	python scripts/Polymerize.py INPUT/L${l}.in
        sed -i s"/THI/L${l} /" ${MID}/L${l}.mol2
	
	echo Looking for the tail of L${l}
	python scripts/find_tail.py ${MID}/L${l}.mol2 > ${MID}/L${l}_tail.txt
	cola="$(head -1 ${MID}/L${l}_tail.txt)"
	

	sed -i s"/LLL/..\/..\/MIDFILES\/L${l}\/L${l}.mol2/" ${NM}
        sed -i s"/CCC/${cola}/" ${NM}
	
	cd scripts/NanoModeler/
	echo Running NanoModeler with L${l}
        python NanoModeler.py > ../../${MID}/NP${l}.log
        cp tmp*/NP.top ../../${MID}/NP${l}.top
        cp tmp*/NP.gro ../../${OUT}/NP${l}.gro
        rm -r tmp*
        cd ../../
	
	sed -i s"/..\/..\/MIDFILES\/L${l}\/L${l}.mol2/LLL/" ${NM}
        sed -i s"/LIG1_C=${cola}/LIG1_C=CCC/" ${NM}
	
	echo Separating parameters files for L${l}
	awk "/\[ atomtypes \]/{f=2;nextnext} /\[ moleculetype \]/{f=0} f" ${MID}/NP${l}.top > ${OUT}/NP${l}.atomtypes
        sed "1,/\[ moleculetype \]/d" ${MID}/NP${l}.top > ${OUT}/NP${l}.itp
        sed -i "1i\[ moleculetype \]" ${OUT}/NP${l}.itp
        sed -i '/\[ system \]/,$d' ${OUT}/NP${l}.itp
        sed -i "3s/.*/NP${l}\t\t3/" ${OUT}/NP${l}.itp
	
done
