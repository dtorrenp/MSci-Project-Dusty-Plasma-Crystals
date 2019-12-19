#!/bin/sh
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=24:mem=16000mb:ompthreads=24
#PBS -N test
module load intel-suite
module load gcc/8.2.0

##FINITE
##cd $HOME/WORK/FINITE
##g++ -std=c++11 $HOME/WORK/FINITE/MSci_project_HPC.cpp -o $HOME/WORK/FINITE/MSci_project_HPC
##exec $HOME/WORK/FINITE/MSci_project_HPC

##PERIODIC
cd $HOME/WORK/PERIODIC
g++ -std=c++11 $HOME/WORK/PERIODIC/MSci_project_periodic_HPC.cpp -o $HOME/WORK/PERIODIC/MSci_project_periodic_HPC
exec $HOME/WORK/PERIODIC/MSci_project_periodic_HPC
