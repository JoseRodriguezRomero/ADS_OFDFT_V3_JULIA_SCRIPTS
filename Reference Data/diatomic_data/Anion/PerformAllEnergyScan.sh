#!/bin/bash

sh PerformEnergyScan.sh H H
sh PerformEnergyScan.sh O O
sh PerformEnergyScan.sh C C
sh PerformEnergyScan.sh N N

sh PerformEnergyScan.sh H C
sh PerformEnergyScan.sh H O

sh PerformEnergyScan.sh C N
sh PerformEnergyScan.sh C O

sh PerformEnergyScan.sh N O
