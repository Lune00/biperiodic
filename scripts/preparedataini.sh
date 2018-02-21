#/bin/bash
a=$(ls -l data/sample/ | tail -n 1) ; fileS=$(echo "$a" | cut -d' ' -f 9) ; 
a=$(ls -l data/cell/ | tail -n 1) ; fileC=$(echo "$a" | cut -d' ' -f 9) ; 
a=$(ls -l data/network/ | tail -n 1) ; fileN=$(echo "$a" | cut -d' ' -f 9) ; 

echo "$fileS $fileC $fileN"

fileS="data/sample/$fileS"
fileC="data/cell/$fileC"
fileN="data/network/$fileN"

echo "$fileS $fileC $fileN"

mkdir -p data_ini
mkdir -p  data_ini/sample
mkdir -p data_ini/cell
mkdir -p  data_ini/network

cp $fileS data_ini/sample/00000reduced.txt
cp $fileC data_ini/cell/00000cell.txt
cp $fileN data_ini/network/00000inter.txt

