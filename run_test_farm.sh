#!/bin/bash

n_values_10_workers=(500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000)
n_values_25_workers=(10000 11000 12000 13000 14000 15000 16000 17000)

workers_10=10
workers_25=25

output_file_10_workers="results_wavefront_10_workers.txt"
output_file_25_workers="results_wavefront_25_workers.txt"

echo "Wavefront Performance Tests (10 Workers)" > $output_file_10_workers
echo "======================================" >> $output_file_10_workers
echo "" >> $output_file_10_workers

for n in "${n_values_10_workers[@]}"
do
  echo "Testing with N=$n and Workers=$workers_10" >> $output_file_10_workers
  echo "--------------------------------------" >> $output_file_10_workers
  
  echo "Running ./wavefront_farm..." >> $output_file_10_workers
  ./wavefront_farm $n $workers_10 >> $output_file_10_workers 2>&1
  echo "" >> $output_file_10_workers

  echo "--------------------------------------" >> $output_file_10_workers
  echo "" >> $output_file_10_workers
done

echo "Test completed for 10 workers. Results saved in $output_file_10_workers"

echo "Wavefront Performance Tests (25 Workers)" > $output_file_25_workers
echo "======================================" >> $output_file_25_workers
echo "" >> $output_file_25_workers

for n in "${n_values_25_workers[@]}"
do
  echo "Testing with N=$n and Workers=$workers_25" >> $output_file_25_workers
  echo "--------------------------------------" >> $output_file_25_workers
  
  echo "Running ./wavefront_farm..." >> $output_file_25_workers
  ./wavefront_farm $n $workers_25 >> $output_file_25_workers 2>&1
  echo "" >> $output_file_25_workers

  echo "--------------------------------------" >> $output_file_25_workers
  echo "" >> $output_file_25_workers
done

echo "Test completed for 25 workers. Results saved in $output_file_25_workers"