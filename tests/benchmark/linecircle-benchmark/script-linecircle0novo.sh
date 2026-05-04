#!/bin/bash

#Script para rodar os primeiros 140 problemas aleatórios (CirLiNovo.txt) no minus-linecircle
#ABR 2026

mkdir licir0novo
for i in $(seq 1 140)
do 
    #Parametros
    cir1=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f1)
    cir1i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f2)
    cir2=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f3)
    cir2i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f4)
    cir3=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f5)
    cir3i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f6)
    li1=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f7)
    li1i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f8)
    li2=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f9)
    li2i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f10)
    li3=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f11)
    li3i=$(cat CirLiNovo.txt | sed -n "$(($i))p" |cut -d' ' -f12)
    #Pontos gt solution
    pt1x=$(cat Pts_linecircleNovo.txt | sed -n "$(($i))p" |cut -d' ' -f1)
    pt1y=$(cat Pts_linecircleNovo.txt | sed -n "$(($i))p" |cut -d' ' -f2)
    pt2x=$(cat Pts_linecircleNovo.txt | sed -n "$(($i))p" |cut -d' ' -f3)
    pt2y=$(cat Pts_linecircleNovo.txt | sed -n "$(($i))p" |cut -d' ' -f4)
    
        script -c "echo $cir1 0 $cir2 0 $cir3 0 $li1 0 $li2 0 $li3 0 $pt1x 0 $pt1y 0 $pt2x 0 $pt2y 0|./../../../../../minus/bin/minus-linecircle -i -gt>> licir0novo/out0novo.txt" temp.txt
	
	echo $cir1 0 $cir2 0 $cir3 0 $li1 0 $li2 0 $li3 0 $pt1x 0 $pt1y 0 $pt2x 0 $pt2y 0 >> licir0novo/configANDgt.txt
	sed -n 's/.*LOG\ 0//p' temp.txt >> licir0novo/sol_step0novo.txt  #num steps pt1
	sed -n 's/.*LOG\ 1//p' temp.txt >> licir0novo/sol_step0novo.txt   #num steps pt2
        grep -q "no valid solution" temp.txt; echo $? >> licir0novo/found_sol0novo.txt    # retorna 1 se o MINUS acha solucao e 0 caso contrario
        grep -q "ground-truth not found" temp.txt; echo $? >> licir0novo/gt_sol0novo.txt    # retorna 1 se a solucao do MINUS= gt solution e 0 caso contrario
done

