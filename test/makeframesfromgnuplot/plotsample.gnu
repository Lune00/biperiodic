set terminal epslatex color colortext solid size 12.5cm,7.5cm standalone background rgb "#00003B"

set output "plotsample.tex"

unset border
unset xtics
unset ytics

set size ratio -1

#Ca va changer aussi non?
#Il va falloir "creer des particules images"
Lx = 6 
Ly = 6

trans_p = 0.1

filep="sample/00003absolute.txt"
filec="network/00003inter.txt"

plot  filep u 3:4:2 w circle lc rgb "#4CB3FF" fs transparent solid 1. noborder  notitle,\
      filep u ($3+Lx):4:2 w circle lc rgb "#00003b" fs transparent solid  trans_p notitle,\
      filep u ($3-Lx):4:2 w circle lc rgb "#00003b" fs transparent solid trans_p notitle,\
      filep u ($3):($4+Ly):2 w circle lc rgb "#00003b" fs transparent solid trans_p notitle,\
      filep u ($3):($4-Ly):2 w circle lc rgb "#00003b" fs transparent solid trans_p notitle,\
      filep u ($3-Lx):($4-Ly):2 w circle lc rgb "#00003b" fs transparent solid trans_p notitle,\
      filep u ($3+Lx):($4-Ly):2 w circle lc rgb "#00003b" fs transparent solid trans_p notitle,\
      filec u 1:2:3:4 w vectors lw 2 lc rgb "yellow" notitle
