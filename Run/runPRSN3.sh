#!/bin/bash
# runPRSN3.sh: Run 100 ns of explicit solvent 2jxb3 model at 310 K.
# This uses Amber programs on Eire.
# Only use first of 20 models in 2jxb
sed -n '1,1091p' 2jxb.pdb >2jxb3.pdb
name=2jxb3
# Step 1: minimize
cat<<eof >leapin
# Use recommended force fields
source leaprc.protein.ff14SB
source leaprc.water.tip3p
nck=loadpdb $name.pdb
solvateOct nck TIP3PBOX 12.0 iso
# Add counterions with solvent replacement. 
# 3 Cl- will neutralize charge, 26 extra bring NaCl to 150 mM.
# There are 9744 waters.
addIons nck Na+ 26 Cl- 29
saveamberparm nck $name.top $name.crd
# Just to check
savepdb nck $name.out.pdb
quit
eof
tleap -f leapin
# Last protein atom 1444.
# Step 1: Minimize solvent
cat << eof >mdin1
Initial 1000 step solvent minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=9, ntb=1,
  ntr=1,
  ntpr=100,
 /
Hold the protein fixed
500.0
RES 1 1444
END
END
eof
mpirun -n 8 sander.MPI -O -i mdin1 -p $name.top -c $name.crd \
-ref $name.crd \
-r $name.min.rst -o $name.min.out -inf $name.inf 
# Step 2: Minimize everything
cat << eof >mdin2
Total 1000 step  minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=9, ntb=1,
  ntpr=100,
 /
END
eof
mpirun -n 8 sander.MPI -i mdin2 -p $name.top -c $name.min.rst \
-r $name.min2.rst -o $name.min2.out -inf $name.inf 
# Step 3: Position-restrained 20 ps NVT dynamic heating
cat << eof >mdin3
NVT position-restrained 20 ps MD heating
 &cntrl
  imin=0, irest=0, ntx=1,
  cut=9, ntb=1,
  ntr=1,
  ntc=2, ntf=2,
  ntt=1, tempi=0.0,
  nstlim = 10000, dt = 0.002,
  ntpr=100,
  nmropt=1, restraint_wt=10.0, restraintmask=':1-1444',
 &end
 &wt type='TEMP0', istep1=0, istep2=10000, value1=0.0, value2=310.0, /
 &wt type = 'END', /
 LISTOUT=POUT
&wt type='END', &end
 /
eof
pmemd.cuda -O -i mdin3 -p $name.top -c $name.min2.rst \
-ref $name.min2.rst \
-r $name.eq1.rst -o $name.eq1.out -inf $name.inf >& /dev/null 
# 10 ns NPT unrestrained
cat << eof >mdin4
Equilibration for 10 ns at 310K, save every 1 ps.
 &cntrl
  imin=0, irest=1, ntx=7,
  ntb=2, pres0=1.0, ntp=1, taup=2.0,
  cut=9, ntr=0,
  ntc=2, ntf=2,
  tempi=310.0, temp0=310.0,
  ntt=1,
  nstlim=5000000, dt=0.002, ioutfm=1, ntwprt=1444
  ntpr=500, ntwx=500,
 /
END
eof
# Run 100 ns MD in 10 ns intervals
for ((i=1;i<=10;i++)); do
((j=i+1))
pmemd.cuda -O -i mdin4 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf \
-inf $name.inf </dev/null
done

