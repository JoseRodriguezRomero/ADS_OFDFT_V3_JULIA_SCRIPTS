#!/bin/bash

sh PerformVibrationalCalc.sh H H $1
sh PerformVibrationalCalc.sh O O $1
sh PerformVibrationalCalc.sh C C $1
sh PerformVibrationalCalc.sh N N $1

sh PerformVibrationalCalc.sh H C $1
sh PerformVibrationalCalc.sh H O $1

sh PerformVibrationalCalc.sh C N $1
sh PerformVibrationalCalc.sh C O $1

sh PerformVibrationalCalc.sh N O $1
