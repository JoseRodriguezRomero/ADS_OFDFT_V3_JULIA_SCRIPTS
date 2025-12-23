#!/bin/bash

bash PerformEnergyScan.sh H H $1
bash PerformEnergyScan.sh O O $1
bash PerformEnergyScan.sh C C $1
bash PerformEnergyScan.sh N N $1

bash PerformEnergyScan.sh H C $1
bash PerformEnergyScan.sh H O $1

bash PerformEnergyScan.sh C N $1
bash PerformEnergyScan.sh C O $1

bash PerformEnergyScan.sh N O $1
