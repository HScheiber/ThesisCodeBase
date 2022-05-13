#!/bin/bash

filename='PES.txt'
path=`pwd`

touch $path/$filename

n=( $(ls -d */) )
for i in "${n[@]}" # loop through first layer of directories
do 
		efile=( $(ls -f $path/$i | grep '.out-001') ) # Find energy file

		if [ -z ${efile+x} ]
		then # if no file exists
			echo -n $i | sed 's/[^0-9|.]*//g' >> $path/$filename #Input lattice parameter
			echo -e '\t' >> $path/$filename # Put a tab
			: # continue
		else # if file exists
			 
			echo -n $i | sed 's/[^0-9|.]*//g' >> $path/$filename #Input lattice parameter
			echo -e -n '\t' >> $path/$filename # Put a tab
			grep 'total all-electron energy =' $path/$i/$efile | sed 's/[^0-9|.|-]*//g' | cut -c 2- >> $path/$filename # total energy
		fi
done
