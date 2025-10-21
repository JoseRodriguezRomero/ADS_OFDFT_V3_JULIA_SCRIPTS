#!/bin/bash

python ParseEnergyScan.py H H
python ParseEnergyScan.py O O
python ParseEnergyScan.py C C
python ParseEnergyScan.py N N

python ParseEnergyScan.py H C
python ParseEnergyScan.py H O

python ParseEnergyScan.py C N
python ParseEnergyScan.py C O

python ParseEnergyScan.py N O
