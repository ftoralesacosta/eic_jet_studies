#!/bin/bash

rm histo_maker
make histo_maker

./histo_maker 0 1.4 1p4T_combined.root &
./histo_maker 0 3.0 3T_combined.root &
