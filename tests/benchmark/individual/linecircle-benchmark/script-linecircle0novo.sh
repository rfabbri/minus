#!/bin/bash

#Script para rodar os primeiros 140 problemas aleatórios ($configurations) no minus-linecircle
#ABR 2026

set -x
set -e
minusdir=/Users/rfabbri/cprg/vxlprg/lemsvpe/minus # TODO: infer this automatically eg from cmake
synthdata=$minusdir/scripts/synthdata/synthdata
minus=$minusdir/bin/minus-chicago
benchmarkdir=$minusdir/tests/benchmark/individual/linecircle-benchmark
configurations=$benchmarkdir/inputs/configuration-specs/*
expname=licir0novo

rm -rf $exp
mkdir  $exp
cd $exp
#mkdir licir0novo

if [[ "`uname`" != Linux ]]; then
  MYOS="OSX"
  mysed=gsed
else
  MYOS="Linux"
  mysed=sed
fi


for i in $(seq 1 140); do  # each configuration
  #Parametros
  cir1=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f1)
  cir1i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f2)
  cir2=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f3)
  cir2i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f4)
  cir3=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f5)
  cir3i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f6)
  li1=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f7)
  li1i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f8)
  li2=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f9)
  li2i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f10)
  li3=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f11)
  li3i=$(cat $configurations | $mysed -n "$(($i))p" |cut -d' ' -f12)
  #Pontos gt solution
  pt1x=$(cat Pts_linecircleNovo.txt | $mysed -n "$(($i))p" |cut -d' ' -f1)
  pt1y=$(cat Pts_linecircleNovo.txt | $mysed -n "$(($i))p" |cut -d' ' -f2)
  pt2x=$(cat Pts_linecircleNovo.txt | $mysed -n "$(($i))p" |cut -d' ' -f3)
  pt2y=$(cat Pts_linecircleNovo.txt | $mysed -n "$(($i))p" |cut -d' ' -f4)
    
  script -c "echo $cir1 0 $cir2 0 $cir3 0 $li1 0 $li2 0 $li3 0 $pt1x 0 $pt1y 0 $pt2x 0 $pt2y 0|./../../../../../minus/bin/minus-linecircle -i -gt>> licir0novo/out0novo.txt" temp.txt
	
	echo $cir1 0 $cir2 0 $cir3 0 $li1 0 $li2 0 $li3 0 $pt1x 0 $pt1y 0 $pt2x 0 $pt2y 0 >> configANDgt.txt
	$mysed -n 's/.*LOG\ 0//p' temp.txt >> sol_step0novo.txt  #num steps pt1
	$mysed -n 's/.*LOG\ 1//p' temp.txt >> sol_step0novo.txt   #num steps pt2
  grep -q "no valid solution" temp.txt; echo $? >> found_sol0novo.txt    # retorna 1 se o MINUS acha solucao e 0 caso contrario
  grep -q "ground-truth not found" temp.txt; echo $? >> gt_sol0novo.txt    # retorna 1 se a solucao do MINUS= gt solution e 0 caso contrario
done

