#!/bin/bash
minfr=$2
fsb=$3
beams=$4
g2la=$5
parallel=$6
out=$7

mkdir out-$7



for i in $8
do

while [ "$(pidof $1 | wc -w)" -ge 8 ]; do 
  sleep 10
done

   mkdir out-$7/BR$i-inst
   ./$1 ../data/BR$i.txt 0 99 $minfr $fsb $beams $g2la $parallel 150 out-$7/BR$i-inst/inst > out-$7/BR$i &
done


