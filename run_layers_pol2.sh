#!/bin/bash

qefile='data/mppc14160.dat'
eff=0.82
pixel=3
nlarr=(1 2 4 8 10 20 30 40 50 75 100)
datfn=${0%.sh}.dat
logfn=${0%.sh}.log

rm -f $logfn $datfn

for nl in ${nlarr[@]}; do
    echo -e "\nЧисло слоев $nl"
    ./farichres -q -e $eff -p $pixel -n 1.05 -T 40 -D 200 -N $nl -mpol2 $qefile | tee -a $logfn | fgrep 'Ошибка угла на трек'
    err=$(tail -3 $logfn | sed -nr 's/Ошибка угла на трек.*: +([0-9.]+)\b.*$/\1/p')
    echo "$nl    $err" >> $datfn
done


