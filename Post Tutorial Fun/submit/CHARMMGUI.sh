#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=1:00:00
#PBS -N openmm_GUI
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
cd CHARMM-gui/charmm-gui-0404029429/openmm/
python -u openmm_run.py -i step4_equilibration.inp -t toppar.str -p step3_input.psf -c step3_input.crd -b sysinfo.dat -orst step4_equilibration.rst -odcd step4_equilibration.dcd > step4_equilibration.out
conda deactivate
