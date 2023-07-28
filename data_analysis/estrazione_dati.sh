#!/bin/bash

counter=0
max=600
step=60
lattice=16


for ((n = 0; n <=$max; n= n + $step))
do
        if [ ${#n} -eq 1 ]
        then
                N="  $n"
                num="$n"
        fi
        if [ ${#n} -eq 2 ]
        then
                N=" $n"
                num="$n"
        fi
        if [ ${#n} -eq 3 ]
        then
                N="$n"
                num="$n"
        fi


        grep -h "n = $N" runQ_${lattice}.log runQ_${lattice}_josh.log runQ_${lattice}_marca.log runQ_${lattice}_josh.log runQ_${lattice}_simo.log runQ_${lattice}_more.log runQ_${lattice}_more_extra.log > log_data_n$num.txt

        grep " Q = " log_data_n$num.txt | awk {'print $15'} | tr "\\\ " " " > log_Q_n$num.txt

        counter=$((counter+1))
done

rm log_data_n*.txt

echo -e "\nDATA EXTRACTED FROM\n\n runQ_${lattice}.log \n runQ_${lattice}_josh.log \n runQ_${lattice}_morelato.log \n runQ_${lattice}_fede.log \n\n$counter FILES log_Q WRITTEN"


gcc -std=c99 -Wall -o histo_Q histo_Q.c

./histo_Q $max $step

echo -e "\n$counter FILES histo_Q CREATED BY histo_Q.c \n\n "

echo -e "\n THE VALUES FOR <Q^2> FOR A LATTICE ${lattice}^4 ARE:\n "


gcc -std=c99 -Wall Q_square.c -o Q_square -lm

./Q_square $max $step
