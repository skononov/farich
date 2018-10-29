#!/bin/bash

qefile='data/mppc14160.dat'
eff=0.82
pixel=3
nlmin=1
nlmax=10

rm -f run.log run.dat

for ((nl=$nlmin; nl<=nlmax; nl++)); do
    echo -e "\nЧисло слоев $nl"
    ./farichres -q -N $nl -e $eff -p $pixel -n 1.05 -T 40 -D 200 -m $qefile | tee -a run.log | fgrep 'Ошибка угла на трек'
    err=$(tail -3 run.log | sed -nr 's/Ошибка угла на трек.*: +([0-9.]+)\b.*$/\1/p')
    echo "$nl    $err" >> run.dat
done


