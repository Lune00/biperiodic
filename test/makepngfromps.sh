#!/bin/bash


mkdir -p framespng

for((i=0;i<=10;i++))
do
	i=$(printf "%05d\n" "$i")
	fileps="analyse/frames/${i}frame.ps"
	filepng="${i}frame.png"
	echo "$fileps -> $filepng"
	ps2pdf -dPDFSETTINGS=/prepress -dEPSCrop $fileps tmp.pdf
	gs -sDEVICE=pngalpha -r3000 -o $filepng tmp.pdf
	mv $filepng framespng
done

