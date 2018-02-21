while(1){set multiplot layout 4,4;
plot 'analyse/SProfile.txt' u 3:2:1 w p palette pt 7 ps 2 ;
plot 'analyse/SProfile.txt' u 4:2:1 w p palette pt 7 ps 2;
plot 'analyse/SProfile.txt' u 5:2:1 w p palette pt 7 ps 2;
plot 'analyse/energy.txt' u 1:2;
plot 'analyse/stress_int.txt' u (-$3/$5) w lp lw 2;
plot 'analyse/compacity.txt' u 1:2 w lp lw 2;
plot 'analyse/interpenetration.txt' u 1:2 w lp lw 2;
pause 2}
