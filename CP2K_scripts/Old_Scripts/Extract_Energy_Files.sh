#!/bin/bash

path=`pwd`
n=( $(ls -d */) )
for i in "${n[@]}" # loop through first layer of directories
do 
	m=( $(ls -d $i*) )
	
	for j in "${m[@]}" # loop through second layer of directories
	do
		efile=( $(ls -f $path/$j | grep '.ener') ) # Find energy file
		
		if [ -z ${efile+x} ]
		then # if no file exists
			: # continue
		else # if file exists
			mkdir -p $path/"out/"$j
			cp $path/$j/$efile $path/"out/"$j/$efile
		fi
	done
done