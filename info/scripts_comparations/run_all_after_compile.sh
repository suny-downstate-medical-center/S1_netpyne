#!/bin/bash

for i in {0..1034}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_input_from_S1_4steps.py $i &
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_Stochkv_gk-gion_input_from_S1_4steps.py $i &
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_V_gStockKv_SynActives.py $i &
   sleep 1
done

