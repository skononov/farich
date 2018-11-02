#!/bin/bash

qefile='data/mppc14160.dat'
eff=0.82
pixel=3
nlarr=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)

logfn=${0%.sh}.log

rm -f $logfn run_layers_fast.dat run_layers_nt.dat run_layers_pol2.dat

for nl in ${nlarr[@]}; do
    echo -e "\n================ Число слоев $nl ================"

    echo "==== Быстрая оптимизация ===="
    ./farichres -q -N $nl -e $eff -p $pixel -n 1.05 -T 40 -D 200 -o "" $qefile | tee -a $logfn | fgrep 'Ошибка угла на трек'
    err=$(tail -3 $logfn | sed -nr 's/Ошибка угла на трек.*: +([0-9.]+)\b.*$/\1/p')
    echo "$nl    $err" >> run_layers_fast.dat

    if [ $nl -le 10 ]; then
       echo "==== NT-оптимизация ===="
       ./farichres -q -N $nl -e $eff -p $pixel -n 1.05 -T 40 -D 200 -mnt -o "" $qefile | tee -a $logfn | fgrep 'Ошибка угла на трек'
       err=$(tail -3 $logfn | sed -nr 's/Ошибка угла на трек.*: +([0-9.]+)\b.*$/\1/p')
       echo "$nl    $err" >> run_layers_nt.dat
    fi

    echo "==== Pol2-оптимизация ===="
    ./farichres -q -N $nl -e $eff -p $pixel -n 1.05 -T 40 -D 200 -mpol2 -o "" $qefile | tee -a $logfn | fgrep 'Ошибка угла на трек'
    err=$(tail -3 $logfn | sed -nr 's/Ошибка угла на трек.*: +([0-9.]+)\b.*$/\1/p')
    echo "$nl    $err" >> run_layers_pol2.dat
done

python plot_layers.py run_layers_fast.dat run_layers_nt.dat run_layers_pol2.dat

