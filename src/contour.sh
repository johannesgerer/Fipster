export GP=/home/data/promotion/Programmierung/gnuplot;
cd $1
source ./plot1.sh

source "$GP/combine.sh"

gnuplot -e "sel=\"contour\";strike=$3;hedge=$4" plot1.sh  $GP/loop3.gp
