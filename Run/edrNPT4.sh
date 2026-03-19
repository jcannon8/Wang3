#!/bin/bash
# edrNPT4.sh: Heat and two 500 ps NPT of edr models on Expanse
#SBATCH -J edr7          # job name
#SBATCH -o %j.out        # output and error file name (%j expands to jobID)
#SBATCH -p gpu-shared
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=93G
#SBATCH --account=umc108
#SBATCH --no-requeue
#SBATCH -t 48:00:00         # run time (hh:mm:ss) - 48 hour (max)
# Load environment variables
module purge
module load gpu
module load slurm
module load openmpi			
module load amber

# Range of models to process (25 per job)
for ((i=7; i<=7;i++)); do
# Start with edm.*.min.rst, output edr.*.eq[1-3].rst, use edm.*.top.
# Heat with protein and lipid restraints, 
# MD with Nck-binding restraints, edr.RST.
name=../edm/edm.$i
oname=edr.$i
# Position-restrained 200 ps NVT dynamic heating to 310 K
cat << eof >mdin3
NVT position-restrained 200 ps MD heating to 310 K
 &cntrl
  imin=0, irest=0, ntx=1,
  ntt=3, gamma_ln=1.0, nmropt=1, ! Ramp temperature 0 to 310
  cut=12.0, fswitch=10.0, ntb=1,
  ntc=2, ntf=2,                  ! SHAKE
  ntr=1,                         ! Protein and lipid restraints
  nstlim=100000, dt=0.002,       ! 200 ps, 2 fs time step
  ntpr=1000, ntwx=5000, ntwr=10000, ! print energy, coords, restart
  ntxo=2, ioutfm=1,              ! Save in NetCDF format
  iwrap=0,
  watnam='WAT', owtnm='O', 
 /
 &ewald
    vdwmeth = 0,
 /
 &wt type='TEMP0', istep1=0, istep2=100000, value1=0.0, value2=310.0, /
 &wt type = 'END', 
 /
 LISTOUT=POUT
 &wt type='END', &end
 /
Protein position restraint
10.0
RES 1 83 84 157 158 213
END
Membrane posres
2.5
FIND
O3 * * CHL1
P * * DOPC
P * * POPC
P * * POPI
P * * DPPC
P * * POPS
P * * DOPE
SEARCH
RES 1 82232
END
END
eof
#
pmemd.cuda -O -i mdin3 -p $name.top -c $name.min.rst \
-ref $name.min.rst \
-r $oname.eq1.rst -o $oname.eq1.out -x $oname.cdf -inf $oname.inf
# Omitted the ">& /dev/null" 
rm $oname.cdf
########################################################### 
# 500 ps production NPT 
cat<<eof>mdin5
500 ps NPT 310 K with edr.RST restraints, 10 ps frames
 &cntrl
  imin=0, irest=1, ntx=5, 
  ntt=3, gamma_ln=1.0, temp0=310.0, ig=-1,
  cut=12.0, fswitch=10.0,
  nstlim=250000, dt=0.002,
  ntc=2, ntf=2,     ! SHAKE
  nmropt=1,         ! Use Nck-binding restraints
  ntpr=5000, ntwr=5000, ntwx=5000,
  iwrap=1,          ! Keep atoms in periodic box
  ntxo=2, ioutfm=1, ! Save in NetCDF format 
  ntwprt=3405       ! Save only 3405 protein atoms
  barostat=2,       ! MC barostat
  ntp=3,            ! semi-isotropic scaling
  pres0=1.0,        ! 1 bar pressure
  csurften=3,       ! Interfaces in xy plane
  gamma_ten=0.0,    ! No surface tension for pure semi-iso scaling
  ninterface=2,     ! Number of interfaces (2 for bilayer)
  ! Set water atom/residue names for SETTLE recognition
  watnam='WAT',     ! Water residues are named WAT
  owtnm='O',        ! Water oxygens are named O
 /
 &ewald
    vdwmeth = 0,
 /
  &wt
    type='END'
 /
DISANG=edr.RST
LISTIN=POUT
LISTOUT=POUT
END
eof
# Run two 500 ps NPT
pmemd.cuda -O -i mdin5 -p $name.top -c $oname.eq1.rst \
-r $oname.eq2.rst -o $oname.eq2.out -x $oname.eq2.cdf -inf $oname.inf 
#
pmemd.cuda -O -i mdin5 -p $name.top -c $oname.eq2.rst \
-r $oname.eq3.rst -o $oname.eq3.out -x $oname.eq3.cdf -inf $oname.inf 
# Use edrNPT5.sh to continue

done
