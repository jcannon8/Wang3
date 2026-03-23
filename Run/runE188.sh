#!/bin/bash
# runE188.sh: Run 100 ns of explicit solvent Lck bound to CD3ε Tyr188 ITAM at 310 K.
# Derived from runLck.sh
cd lck
# Script anEdt22.tcl surveyed the edt trajectories and found frames with
# low RMSD with two gauche- angles.
#sort -nk3 edtLck4.dat|head -n20
# 35  1941  1.61   -64.4   -65.7  2.29   179.4    90.1
# sort -nk6 edtLck4.dat|head -n20 
# 84 13851  3.55   -77.5    99.3  1.47   -50.5   -75.5
#
# Use vmd -dispdev none
# Lck substrate atomselection
set lckpep [mol new /home/jcannon/lck/Lck-Peptide.pdb]
set lck [atomselect $lckpep "chain A"]
set lcksub [atomselect $lckpep all]
set lckH [atomselect $lckpep "chain A and mass>2"]; # Lck heavy atoms
# Choose residues for fitting
set t 4; # The 9 residues around pTyr for fitting.
set w 7; # How far away from pTyr are Lck collisions allowed? 
set c [expr 285- $t]; # first Lck substrate residue to fit
set d [expr 285+ $t]; # last Lck substrate residue to fit
set pSel [atomselect $lckpep "name CA and chain B and resid $c to $d"]

# Get CD3e ITAM 188 and ITAM 199 coordinates fit to Lck-bound peptide.
set mod 35
set tcr [mol new /t/tcr16/edt.prot.psf]
# The NPT MD started with edt.*.eq2.
for {set i 2} {$i<=53} {incr i} {
  set mdFile /t/tcr16/edt/edt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
#
set n 188
set a [expr $n-7]
set b [expr $n+7]
set e [expr $n-4]
set g [expr $n+4]
# Only save heavy atoms of E188
set itam188 [atomselect $tcr "chain E and resid $a to $b and mass>2" frame 1941]
set iSel [atomselect $tcr "chain E and name CA and resid $e to $g" frame 1941]
# fit ITAM to substrate and move ITAM
set M [measure fit $iSel $pSel]
$itam188 move $M
set rms [measure rmsd $pSel $iSel]
$itam188 writepdb E188.pdb
###
# Now, ITAM 199
set mod 84
set tcr [mol new /t/tcr16/edt.prot.psf]
# The NPT MD started with edt.*.eq2.
for {set i 2} {$i<=53} {incr i} {
  set mdFile /t/tcr16/edt/edt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
#
set n 199
set a [expr $n-7]
set b [expr $n+7]
set e [expr $n-4]
set g [expr $n+4]
set itam199 [atomselect $tcr "chain E and resid $a to $b and mass>2" frame 13851]
set iSel [atomselect $tcr "chain E and name CA and resid $e to $g" frame 13851]
# fit ITAM to substrate and move ITAM
set M [measure fit $iSel $pSel]
$itam199 move $M
set rms [measure rmsd $pSel $iSel]
$itam199 writepdb E199.pdb
###########################################################
# Build E188 and run MD with Amber programs.
name=E188
# lck2.pdb from runLck.sh
# Fix atom name in E188.pdb
sed 's/CD  ILE/CD1 ILE/' E188.pdb > E188a.pdb
# Step 1: build system
cat<<eof >leapin
# Use recommended force fields
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.DNA.OL15
loadamberparams atp.frcmod
loadamberprep atp2.prep
lck=loadpdb lck2.pdb
pep=loadpdb E188a.pdb
atp=loadpdb ATP.pdb
lckPep=combine {lck pep atp}
check lckPep
# ^That gave charge of -7.
solvateOct lckPep TIP3PBOX 12.0 iso
# ^That added 13,036 waters
# Add counterions with solvent replacement. 
# Number of ions = 0.01867 * (Desired Concentration) * (Number of waters)
# =37
# 7 Na+ will neutralize charge, 44 extra bring NaCl to 150 mM.
addIons lckPep Na+ 44 Cl- 37
saveamberparm lckPep $name.top $name.crd
# Just to check
savepdb lckPep $name.out.pdb
quit
eof
tleap -f leapin
# Order: Lck atoms 1-4401, resid 1-273;
# E188 atoms 4402-4641, resid 274-288;
# ATP atoms 4642-4684, resid 289
# Last protein atom 4641.
# Step 2: Minimize solvent
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
RES 1 4641
END
END
eof
pmemd.cuda -O -i mdin1 -p $name.top -c $name.crd \
-ref $name.crd \
-r $name.min.rst -o $name.min.out -inf $name.inf 
# Step 3: Minimize everything
cat << eof >mdin2
Total 1000 step  minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=9, ntb=1,
  ntpr=100,
 /
END
eof
pmemd.cuda -i mdin2 -p $name.top -c $name.min.rst \
-r $name.min2.rst -o $name.min2.out -inf $name.inf 
# Step 4: Position-restrained 20 ps NVT dynamic heating
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
  nmropt=1, restraint_wt=10.0, restraintmask=':1-4641',
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
# Step 5: 10 ns NPT unrestrained
cat << eof >mdin4
Equilibration for 10 ns at 310K, save every 10 ps.
 &cntrl
  imin=0, irest=1, ntx=7,
  ntb=2, pres0=1.0, ntp=1, taup=2.0,
  cut=9, ntr=0,
  ntc=2, ntf=2,
  tempi=310.0, temp0=310.0,
  ntt=1,
  nstlim=5000000, dt=0.002, ioutfm=1, ntwprt=4641
  ntpr=5000, ntwx=5000,
 /
END
eof
# Run 100 ns NPT MD in 10 ns intervals
for ((i=1;i<=10;i++)); do
((j=i+1))
pmemd.cuda -O -i mdin4 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf \
-inf $name.inf </dev/null
done


