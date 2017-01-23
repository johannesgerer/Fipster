export GP=/home/data/promotion/Programmierung/gnuplot;
cd $1
gnuplot -e "sel=\"$2\"" plot1.sh  $GP/loop3.gp
