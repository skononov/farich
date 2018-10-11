#!/bin/bash

qefile='data/mppc14160.dat'
eff=0.88
pixel=3
nlmin=1
nlmax=10

echo '' > run.log
echo '' > run.dat

for ((nl=$nlmin; nl<=nlmax; nl++)); do
    echo -e "\nЧисло слоев $nl"
    ./farichres -q -N $nl -e $eff -p $pixel -n 1.05 -T 40 -m $qefile | tee -a run.log | fgrep 'Ошибка угла на трек:'
    err=$(tail -2 run.log | sed -nr 's/Ошибка угла на трек: +([0-9.]+) .*/\1/p')
    echo "$nl    $err" >> run.dat
done


