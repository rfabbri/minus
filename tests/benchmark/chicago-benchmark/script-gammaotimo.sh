#!/bin/bash

#Script para rodar os primeiros 140 problemas aleatórios gerados a partir dos dados sintéticos
#MAIO 2025

minusdir=/home/juliana/Documents/doutorado/minus
synthdata=$minusdir/scripts/synthdata/synthdata
minus=$minusdir/bin/minus-chicago

for i in $(seq 8 8)
do 
    mkdir frpt-P$i
    fr1=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f1)
    fr2=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f2)
    fr3=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f3)
    pt1=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f4)
    pt2=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f5)
    pt3=$(cat ../FRePTcorrigido.txt | sed -n "$(($i))p" |cut -d' ' -f6)
    for j in $(seq 1 5)
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

