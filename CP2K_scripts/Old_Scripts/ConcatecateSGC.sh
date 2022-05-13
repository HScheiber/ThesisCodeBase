#!/bin/bash

path=`pwd`
n=( $(ls -d */) )
for i in "${n[@]}" # loop through first layer of directories
do 
	m=( $(ls -d $i*) )
	
	ITER=0
	for j in "${m[@]}" # loop through second layer of directories
	do
		efile=( $(ls -f $path/$j | grep 'frc-1.xyz') ) # Find energy file
		
		if [ $ITER -eq 0 ]
		then
			outfile=${i%/*}"-frc.xyz" # Make output file
			cat $path/$j/$efile > $path/$i$outfile
		else
			cat $path/$j/$efile >> $path/$i$outfile # concatecate into output file
		fi
		
		ITER=$(expr $ITER + 1) #iterate by 1
	done
done