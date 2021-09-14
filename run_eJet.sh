#!/bin/bash

rm eJet_Distributions
make eJet_Distributions

./eJet_Distributions 3T_combined.root 3.0 &
./eJet_Distributions 1p4T_combined.root 1.4 &
