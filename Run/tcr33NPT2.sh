#!/bin/bash
# tcr33NPT2.sh: 20 ns production NPT of tcr33 on Rockfish
# run sbatch --dependency=afterok:7295307 tcr33NPT2.sh on Rockfish. 
#SBATCH -J tcr35f.10       # job name
#SBATCH -o %j.out         # output and error file name (%j expands to jobID)
#SBATCH -p a100           # Queue name: defq or a100
#SBATCH -N 1
#SBATCH --cpus-per-gpu=1
#SBATCH --gres=gpu:1
#SBATCH -x gpu01          # Exclude the V100 node, gpu01
#SBATCH -A mcb140208_gpu
#SBATCH -t 36:00:00       # run time (hh:mm:ss) - 48 hour (max)
# Load environment variables
module load amber/20

name=tcr33
# Starting eq* (increase by 4)
i=10
((j=i+1))

echo "Starting MD with $name.eq$i.rst"
# 5 ns production NPT on Rockfish with restart saving every 1 ns.
cat<<eof>mdin7
A NPT simulation for common production-level simulations
 &cntrl
    imin=0,        ! No minimization
    irest=1,       ! This is a restart of an old MD simulation
    ntx=5,         ! So our inpcrd file has velocities
    !! Temperature control
    ntt=3,         ! Langevin dynamics
    gamma_ln=1.0,  ! Friction coefficient (ps^-1)
    temp0=310,   ! Target temperature
    !! Potential energy control
    cut=12.0,      ! nonbonded cutoff, in Angstroms
    fswitch=10.0,  ! Force-based switching
    !! MD settings
    nstlim=2500000, ! 2.5M steps, 5 ns total
    dt=0.002,      ! time step (ps)
    !! SHAKE
    ntc=2,         ! Constrain bonds containing hydrogen
    ntf=2,         ! Do not calculate forces of bonds containing hydrogen
    !! Control how often information is printed
    ntpr=5000,     ! Print energies every 5000 steps, 10 ps
    ntwx=5000,     ! Print coordinates every 5000 steps, 10 ps, to the trajectory
    ntwr=10000,    ! Print a restart file every 10K steps
    ntxo=2,        ! Write NetCDF format
    ioutfm=1,      ! Write NetCDF format 
    ntwprt=16877   ! save protein atoms and first two CHL1 in trajectory
    !! Wrap coordinates when printing them to the same unit cell
    iwrap=1,
    !! Constant pressure control
    barostat=2,    ! MC barostat
    ntp=3,         ! 1=isotropic, 2=anisotropic, 3=semi-isotropic w/ surften
    pres0=1.0,     ! Target external pressure, in bar
    !! Constant surface tension (needed for semi-isotropic scaling). Uncomment
    !! for this feature. csurften must be nonzero if ntp=3 above
    csurften=3,    ! Interfaces in 1=yz plane, 2=xz plane, 3=xy plane
    gamma_ten=0.0, ! Surface tension (dyne/cm). 0 gives pure semi-iso scaling
    ninterface=2,  ! Number of interfaces (2 for bilayer)
    !! Set water atom/residue names for SETTLE recognition
    watnam='WAT',  ! Water residues are named WAT
    owtnm='O',     ! Water oxygens are named O
 /
 &ewald
    vdwmeth = 0,
 /
eof


pmemd.cuda -O -i mdin7 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf >& /dev/null
((i=i+1))
((j=i+1))
pmemd.cuda -O -i mdin7 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf >& /dev/null
((i=i+1))
((j=i+1))
pmemd.cuda -O -i mdin7 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf >& /dev/null
((i=i+1))
((j=i+1))
pmemd.cuda -O -i mdin7 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf >& /dev/null


