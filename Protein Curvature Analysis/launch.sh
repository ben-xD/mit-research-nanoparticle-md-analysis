#!/bin/bash
#SBATCH -J Contact3
#SBATCH -o contact3.out
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p regular-cpu

bash Contact3_Frame_Calc.bash

