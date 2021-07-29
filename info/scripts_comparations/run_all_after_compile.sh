#!/bin/bash

for i in {1018..1034}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_somatic_efel.py $i 
done

