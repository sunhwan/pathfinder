#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -o out
#PBS -e err
#PBS -l walltime=12:00:00
#PBS -q roux
#PBS -N anmpathway

cd $PBS_O_WORKDIR

python taskmanager.py {{ uuid }}

