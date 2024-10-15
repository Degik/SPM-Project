#!/bin/bash

n_values=(500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 13000 13500 14000 14500 15000)
workers=20
output_file="results_wavefront.txt"

echo "Wavefront Performance Tests" > $output_file
echo "=============================" >> $output_file
echo "" >> $output_file

for n in "${n_values[@]}"
do
  echo "Testing with N=$n and Workers=$workers" >> $output_file
  echo "--------------------------------------" >> $output_file
  
  echo "Running ./wavefront_pf..." >> $output_file
  ./wavefront_pf $n $workers >> $output_file 2>&1
  echo "" >> $output_file

  echo "Running ./wavefront_pf_cache..." >> $output_file
  ./wavefront_pf_cache $n $workers >> $output_file 2>&1
  echo "" >> $output_file

  echo "Running ./wavefront_farm..." >> $output_file
  ./wavefront_farm $n $workers >> $output_file 2>&1
  echo "" >> $output_file

  echo "--------------------------------------" >> $output_file
  echo "" >> $output_file
done

echo "Test completed. Results saved in $output_file"
