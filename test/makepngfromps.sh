#!/bin/bash


mkdir -p framespng

psfiles="analyse/frames/*.ps"
#for file in "${folderps}/*.ps"
for file in $psfiles
#for((i=0;i<=10;i++))
do
	[ -e "$file" ] || continue
	echo $file
	filepng="`echo "$file" | cut -d'.' -f1`.png"
	echo $filepng
#	i=$(printf "%05d\n" "$i")
#	fileps="analyse/frames/${i}frame.ps"
#	filepng="${i}frame.png"
#	echo "$fileps -> $filepng"
	ps2pdf -dPDFSETTINGS=/prepress -dEPSCrop $file tmp.pdf
	gs -sDEVICE=pngalpha -r2000 -o $filepng tmp.pdf
	mv $filepng framespng
done

