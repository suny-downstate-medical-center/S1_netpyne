#!/bin/bash

for i in {275..500}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i &
   sleep 30
done

for i in {501..501}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i 
done

for i in {502..700}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i &
   sleep 30
done

for i in {701..701}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i 
done

for i in {702..1034}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i &
	sleep 30
done
