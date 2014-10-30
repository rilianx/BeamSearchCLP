#!/bin/bash
minfr=$2
fsb=$3
out=$4

mkdir out-$4



for i in $5
do

while [ "$(pidof $1 | wc -w)" -ge 8 ]; do 
  sleep 10
done

   ./$1 ../data/BR$i.txt 0 99 $minfr $fsb 150 no > out-$4/BR$i &
done


