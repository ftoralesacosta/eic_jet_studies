#!/bin/bash

rm DeltaR_eJet_Distributions
make DeltaR_eJet_Distributions

./DeltaR_eJet_Distributions 3T_combined.root 3.0 &
./DeltaR_eJet_Distributions 1p4T_combined.root 1.4 &
