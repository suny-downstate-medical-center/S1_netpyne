#!/bin/bash

for i in {0..5}; do 
   echo $i | time -f "Elapsed Time = %E, CPU Usage = %P, CPU Time= %S, Job = $i" -o time.txt --append ipython compare_net.py $i &
done

