#!/bin/bash


mkdir -p frames

for((i=0;i<=3;i++))
do
	i=$(printf "%05d\n" "$i")
	frame="${i}frame.eps"
	fc="filec=\"network\/${i}inter.txt\""
	fp="filep=\"sample\/${i}absolute.txt\""
	echo $fp
	echo $fc
	sed  "s/fp/$fp/g" plotsampleVOID.gnu > tmp.gnu
	sed  "s/fc/$fc/g" tmp.gnu > plotsample.gnu 

	geps plotsample.gnu
	mv plotsample.eps $frame
	mv $frame frames


done
