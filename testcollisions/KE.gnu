set terminal epslatex color colortext solid size 12.5cm,7.5cm standalone 
#background rgb "#00003B"
set size ratio 1/1.61
set key spacing 1. width 1. bottom
set output 'KE.tex'

set palette rgb 33,13,10

unset colorbox

set ylabel '\Large Energy'
set xlabel '\Large $t$' offset 2,0

#set xr[5900:6000]
#set yr [0:8]
#set xr[275:325]

plot 'analyse/energy.txt' u 1:2 w lp lw 2 title 'Kinetic',\
     '' u 1:4 w lp lw 2 title 'Elastic',\
     '' u 1:3 w lp lw 2 title 'Rotation',\
     '' u 1:5 w lp lw 2 title 'Total'




