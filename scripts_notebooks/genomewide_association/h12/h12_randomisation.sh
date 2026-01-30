#!/bin/bash

pops=(
	NonPBO
	PBO
)
chroms=(
	2RL
	3RL
	X
)

for chrom in ${chroms[@]}
do
	for pop in ${pops[@]}
	do
		python h12.py $pop $chrom 1000 &
	done
	#wait
done


