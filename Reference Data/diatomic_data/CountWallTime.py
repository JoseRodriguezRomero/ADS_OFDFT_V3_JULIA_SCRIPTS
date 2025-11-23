import os
import math
from pathlib import Path


def Round(value, digits = 3):
    factor = 10 ** digits
    return math.ceil(value * factor) / factor


def TimeToSeconds(t):
    h, m, s = t.split(":")
    return int(h) * 3600 + int(m) * 60 + float(s)


def GetWallTimeFromFile(file_lines):
    for line in file_lines:
        if "Psi4 wall time for execution:" in line:
            line_splitted = line.split()
            return TimeToSeconds(line_splitted[-1])
            
    return 0.0


def GetEnergyScansWallTime(ion_state):
    files = list(Path(ion_state + "/EnergyScans/").glob("*.out"))
    wall_time = 0.0
    
    for file_name in files:
        file = open(file_name,"r")
        wall_time = wall_time + GetWallTimeFromFile(file.readlines())
        
    return wall_time


n_time = GetEnergyScansWallTime("Neutral")
c_time = GetEnergyScansWallTime("Cation")
a_time = GetEnergyScansWallTime("Anion")
t_time = n_time + c_time + a_time

print("Neutral scans took: " + str(Round(n_time / 3600.0)) + " Hours")
print(" Cation scans took: " + str(Round(c_time / 3600.0)) + " Hours")
print("  Anion scans took: " + str(Round(a_time / 3600.0)) + " Hours")
print(" ")
print("   Total wall time: " + str(Round(t_time / 3600.0)) + " Hours")
