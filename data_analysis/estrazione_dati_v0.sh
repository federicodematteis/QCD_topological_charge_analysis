#!/bin/bash

#file n=0

grep -h "n =   0" runQ_12.log > log_data_n0.txt 

grep " Q = " log_data_n0.txt | awk {'print $15'} | tr "\\\ " " " > log_Q_n0.txt


counter=1

#files n=40,80

for ((n = 40; n <=80; n= n + 40))
do

	grep -h "n =  $n" runQ_12.log > log_data_n$n.txt 

	grep " Q = " log_data_n$n.txt | awk {'print $15'} | tr "\\\ " " " > log_Q_n$n.txt
	
	counter=$((counter+1))
done

#files n=120,...,400

for ((n = 120; n <=400; n= n + 40))
do

	grep -h "n = $n" runQ_12.log > log_data_n$n.txt 

	grep " Q = " log_data_n$n.txt | awk {'print $15'} | tr "\\\ " " " > log_Q_n$n.txt
	
	counter=$((counter+1))
done

rm log_data_n*.txt

echo -e "\nDATA EXTRACTED FROM runQ_12.log\n\n$counter FILES log_Q WRITTEN\n "

gcc -std=c99 -Wall -o histo_Q histo_Q.c 

./histo_Q
