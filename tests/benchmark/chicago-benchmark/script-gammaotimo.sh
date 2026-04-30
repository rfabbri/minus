#!/bin/bash
# Author: Juliana Ventura
# See GEMINI.md in benchmarksdir below

#Script para rodar os primeiros 140 problemas aleatórios gerados a partir dos dados sintéticos
#MAIO 2025

set -x
minusdir=/home/juliana/Documents/doutorado/minus   # TODO: infer this automatically eg from cmake
synthdata=$minusdir/scripts/synthdata/synthdata
minus=$minusdir/bin/minus-chicago
benchmarkdir=$minusdir/tests/benchmark
configurations=$benchmarkdir/FRePTcorrigido.txt

for i in $(seq 1 2)
do 
    mkdir frpt-P$i
    fr1=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f1)
    fr2=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f2)
    fr3=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f3)
    pt1=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f4)
    pt2=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f5)
    pt3=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f6)
    for j in $(seq 1 5) # number of runs
    do 
        script -c "$synthdata $fr1 $fr2 $fr3 $pt1 $pt2 $pt3 0 1 |$minus -i -gt> out.txt" temp.txt

#	cat out.txt | cut -d '(' -f 2| cut -d ")" -f 1| sed -n '1,114p'>> frpt-P$i/vec_gamma-P$i.txt
       cat temp.txt | cut -d' ' -f 3 | sed -n '50,361p' >> frpt-P$i/sol_step-P$i.txt
       grep "isvalid)" temp.txt | cut -d')' -f2 >> frpt-P$i/realsol-P$i.txt
       echo '313,10'>> frpt-P$i/realsol-P$i.txt
       solution=$(grep index temp.txt | cut -d : -f 2)
       if [ $solution -eq  ]
           then solution=313
       fi
       echo $solution>> frpt-P$i/gt_sol-P$i.txt
       cat temp.txt | cut -d' ' -f 2-5| sed -n '2,3p' >> frpt-P$i/frpt-$i.txt
    done
done

