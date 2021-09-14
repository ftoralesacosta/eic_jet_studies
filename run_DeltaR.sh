#!/bin/bash

rm DeltaR_histo_maker
make DeltaR_histo_maker

./DeltaR_histo_maker 0 1.4 1p4T_combined.root &
./DeltaR_histo_maker 0 3.0 3T_combined.root &
