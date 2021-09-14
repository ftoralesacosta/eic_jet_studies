#!/bin/bash
rm DeltaR_eJet_Distributions_MinpT
make DeltaR_eJet_Distributions_MinpT
array_fields=(1.4 3.0)
array_files=("1p4T_combined.root" "3T_combined.root")
array_thresholds=(0.0 0.1 0.3 0.5 0.7 1.0)

for B in "${!array_files[@]}";
do
  echo "${array_files[B]}" "${array_fields[B]}" 
  for pT in "${array_thresholds[@]}" #Min pT Thresholds in GeV 
  do
    # ./eJet_Distributions $B $pT
    command="./DeltaR_eJet_Distributions_MinpT ${array_files[B]} ${array_fields[B]} ${pT}"
    echo $command
    $command &
  done
done
