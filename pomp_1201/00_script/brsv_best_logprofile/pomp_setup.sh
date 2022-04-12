#!/bin/bash

scirptf=$(find . -name "hhs*_loglik.R")
#find . -name "*_script.R"
#echo $scirptf


for i in $scirptf
do
	
	runname=$(sed "s/_loglik.R//g;s#./hhs##g" <<< $i)
	echo $runname
	#sed -i -e '10s/.*/inc_data_add %.>%/' $i
	#sed -i -e 's#co-infect#co_infect#g' $i
	#sed -i -e 's#365/30#365/180#g' $i
	#sed -i -e 's#w1=1#w1=2#g' $i
	#sed -i -e 's#w2=1#w2=2#g' $i
	cat ./sub.sh | sed "s#jobname#$runname#" > pomp_sub.sh
	echo "R CMD BATCH $i" >> pomp_sub.sh
	sbatch ./pomp_sub.sh
	echo "$i is submitted"
	
done	





