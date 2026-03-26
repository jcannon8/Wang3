# makeTcr38d.tcl: Add CD3 zeta dimer to tcr38.* to make tcr38.*.**.
# Use output from makeTcr38c.tcl to choose which zeta3 frames to add.
# Derived from makeTcr37b.tcl
# Use vmd -dispdev none 
cd /t/tcr38b
set mod 11;  # output tcr38.[0-8] model number
set aRef 2; # chain a zeta reference (0-8), actually use only 0-2
set bRef 2; # chain b zeta reference (0-8), actually use only 0-2
#  0  3212  6680 |  0  0 | 0 0 0 0 0 0 0 0 0
set aFrame 3212; # zeta3 frame to add to chain A
set bFrame 6680; # zeta3 frame to add to chain B
set inMod tcr38.$mod; # input tcr38 parent model to add zeta chains
set outMod tcr38.$mod.$aRef$bRef; # output model
# Use 211.2 ns, eq45a autoimaged frame for tcr38.0, tcr38.3
# use 200 ns, eq49a autoimaged frame for tcr38.2.??
set eq eq45
set eqa eq45a
set tcr [mol new /t/tcr38/$inMod.psf]
mol addfile /t/tcr38/$inMod.$eqa.rst type netcdf waitfor all
####################
# Add frame aFrame zeta3 CT to chain a. Similar to AddTail15b.tcl
set cd3a [mol new /home/jcannon/tcr15/zeta.top type parm7]
mol addfile /home/jcannon/tcr15/zeta3.c$aRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/zeta3.c$aRef.eq3.cdf type netcdf waitfor all
# Cytoplasmic tail resid 6-112 of zeta3. New addition to tcr15 are residues 2-112.
# CD3 epsilon: 6JXR chain a resid 22 to 57 in tcr15
# Fit between tcr and zeta3 uses LYS PHE SER ARG residues in each.
# Align resid 2-5 of cd3a and resid 54 to 57 of tcr15
set End [atomselect $tcr "chain A and name CA and resid 54 to 57"]
set Start [atomselect $cd3a "resid 2 to 5 and name CA" frame $aFrame]
set cd3CD [atomselect $cd3a all frame $aFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move CD3 CT to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3aN [atomselect $tcr "chain A and resid 22 to 55"]
set cd3aC [atomselect $cd3a "resid 4 to 112" frame $aFrame]
$cd3aN set chain a
$cd3aC set chain a
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $cd3a "resid 4 to 112"] get index] {
  set r [[atomselect $cd3a "index $i"] get resid]
  [atomselect $cd3a "index $i"] set resid [expr $r + 52]
}
# Fuse chain using bash commands.
$cd3aN writepdb cd3aN.pdb
$cd3aC writepdb cd3aC.pdb
set out aChain.pdb 
exec cat cd3aN.pdb > $out 
exec cat cd3aC.pdb >> $out
exec echo "TER" >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
# Remove terminal carboxy oxygens, psfgen will put them back with guesscoord
exec sed -i /OT1/d $out; exec sed -i /OT2/d $out; 
file delete cd3aN.pdb; file delete cd3aC.pdb;
###########################################################
# Add frame bFrame zeta3 CT to chain b
set cd3b [mol new /home/jcannon/tcr15/zeta.top type parm7]
mol addfile /home/jcannon/tcr15/zeta3.c$bRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/zeta3.c$bRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $cd3b get numframes]
# Cytoplasmic tail resid 4-112 of zeta3. New addition to tcr15 are residues 2-112.
# CD3 epsilon: 6JXR chain b resid 24 to 55 in tcr15
# Fit between tcr and zeta3 uses VAL LYS PHE residues in each.
# In retrospect, zeta3 should have started with one more N-terminal residue so
# that its start was ARG VAL LYS PHE, instead of VAL LYS PHE.
# Increase to "CA C N", the atoms aligned because it is short.
# Align resid 1-3 of cd3b and resid 53 to 55 of tcr15
set End [atomselect $tcr "chain B and name CA C N and resid 53 to 55"]
set Start [atomselect $cd3b "name CA C N and resid 1 to 3" frame $bFrame]
set cd3CD [atomselect $cd3b all frame $bFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move CD3 CT to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3bN [atomselect $tcr "chain B and resid 24 to 54"]
set cd3bC [atomselect $cd3b "resid 3 to 112" frame $bFrame]
$cd3bN set chain b
$cd3bC set chain b
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $cd3b "resid 3 to 112"] get index] {
  set r [[atomselect $cd3b "index $i"] get resid]
  [atomselect $cd3b "index $i"] set resid [expr $r + 52]
}
# Fuse chain using bash commands.
$cd3bN writepdb cd3bN.pdb
$cd3bC writepdb cd3bC.pdb
set out bChain.pdb 
exec cat cd3bN.pdb > $out 
exec cat cd3bC.pdb >> $out
exec echo "TER" >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
# Remove terminal carboxy oxygens, psfgen will put them back with guesscoord
exec sed -i /OT1/d $out; exec sed -i /OT2/d $out; 
file delete cd3bN.pdb; file delete cd3bC.pdb
###########################################################
# Save all the other parts
[atomselect $tcr "chain D"] writepdb dChain.pdb
[atomselect $tcr "chain E"] writepdb eChain.pdb
[atomselect $tcr "chain F"] writepdb fChain.pdb
[atomselect $tcr "chain G"] writepdb gChain.pdb
[atomselect $tcr "chain M"] writepdb mChain.pdb
[atomselect $tcr "chain N"] writepdb nChain.pdb
# Save S1 structural cholesterol and S2 DOPC. 
set S1 [atomselect $tcr "chain S and resname CHL1"] 
$S1 set segid S1; $S1 writepdb S1.pdb
set S2 [atomselect $tcr "chain S and resname DOPC"]
$S2 set segid S2; $S2 writepdb S2.pdb
[atomselect $tcr "chain J"] writepdb mem.pdb
###########################################################
# Remove cytoplasmic water that collides with new chain A and B CT.
set water [atomselect $tcr water]
# Set beta=1 for water to keep.
$water set beta 1
# Test aChain collision with water
set aChain [mol new aChain.pdb]
set aProt [atomselect $aChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $aProt $water]
set wats [lindex $con 1]
# Use atomselection "same <keyword> as <selection>"
[atomselect $tcr "same residue as index $wats"] set beta 2
[atomselect $tcr "beta 2"] num; # 1014 waters removed
# Test bChain collision with water
set bChain [mol new bChain.pdb]
set bProt [atomselect $bChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $bProt $water]
set wats [lindex $con 1]
# Use atomselection "same <keyword> as <selection>"
[atomselect $tcr "same residue as index $wats"] set beta 3
[atomselect $tcr "beta 3"] num; # 1074 waters removed
# Saved waters still marked beta=1 above.
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr "water and beta 1"]
$saveWat set chain z
set Res [lsort -unique -integer [$saveWat get residue]]
set last [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$last} {incr n} {
  set j [expr ($n-1) * $chunk +1]
  set k [expr $j + $chunk -1]
  puts "$n $j $k"
  set toSave [atomselect $tcr "water and beta 1 and residue $j to $k"]
  # Need to adjust residue# so that they go from 1 to $chunk in each file.
  # That ensures they will be unique within each psfgen segment.
  unset -nocomplain z
  set nRes [expr [$toSave num] / 3]
  for {set i 1} {$i<=$nRes} {incr i} {lappend z $i $i $i}
  $toSave set resid $z
  $toSave writepdb water$n.pdb
}
# water1.pdb to water[$n-1].pdb saved
###########################################################
# Adjust ions
set ions [atomselect $tcr "name SOD CLA"]
$ions set chain i
#$ions writepdb ions.pdb
# Initially save all ions. That gave netCharge=+5.
[atomselect $tcr "name SOD"] num; # 542
[atomselect $tcr "name CLA"] num; # 320
# Remove 23 SOD to neutralize
set resSOD [[atomselect $tcr "name SOD"] get resid]
set saveSOD [lrange $resSOD 0 end-5]
# Now save fewer SOD
[atomselect $tcr "name CLA or (name SOD and resid $saveSOD)"] writepdb ions.pdb
###########################################################
# Use psfgen to build psf and pdb outputs.
package require psfgen
# Load the Charmm32m topology and parameters necessary for this system.
topology /home/jcannon/charm/toppar/top_all36_prot.rtf
topology /home/jcannon/charm/toppar/par_all36m_prot.prm
topology /home/jcannon/charm/toppar/top_all36_lipid.rtf
topology /home/jcannon/charm/toppar/par_all36_cgenff.prm
topology /home/jcannon/charm/toppar/par_all36_carb.prm
topology /home/jcannon/charm/toppar/par_all36_lipid.prm
topology /home/jcannon/charm/toppar/toppar_all36_lipid_cholesterol.str
topology /home/jcannon/charm/toppar/toppar_all36_carb_glycolipid.str
topology /home/jcannon/charm/toppar/toppar_all36_lipid_inositol.str
topology /home/jcannon/charm/toppar/toppar_water_ions.str
# The above reports many errors but they are OK (I hope so).
# Load segments with coordinates. 
# This defines the order in this model (same as 6JXR and tcr15 order)
pdbalias residue HIE HSD
# Load segments with coordinates.
# This defines the order in this model (same as 7FJD, tcr30, and tcr38 order)
segment a {pdb aChain.pdb}; coordpdb aChain.pdb a
segment b {pdb bChain.pdb}; coordpdb bChain.pdb b
segment d {pdb dChain.pdb}; coordpdb dChain.pdb d
segment e {pdb eChain.pdb}; coordpdb eChain.pdb e
segment f {pdb fChain.pdb}; coordpdb fChain.pdb f
segment g {pdb gChain.pdb}; coordpdb gChain.pdb g
segment m {pdb mChain.pdb}; coordpdb mChain.pdb m
segment n {pdb nChain.pdb}; coordpdb nChain.pdb n
# Structural cholesterol CHL1 201, DOPC 190, membrane, ions
segment s1 {pdb S1.pdb}; coordpdb S1.pdb s1
segment s2 {pdb S2.pdb}; coordpdb S2.pdb s2
segment j {pdb mem.pdb}; coordpdb mem.pdb j
segment k {pdb ions.pdb}; coordpdb ions.pdb k
# Load in water
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
# Water is in multiple pdb files
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
# Only one disulfide remains from 6jxr.pdb, between zeta chains.
patch DISU A:32 B:32
# Very important to add missing H atom coordinates for chain termini!
guesscoord
writepsf $outMod.psf
writepdb $outMod.pdb
# Check charge to see if ions need alteration.
set tcr3 [mol new $outMod.psf]
mol addfile $outMod.pdb waitfor all
set sel [atomselect $tcr3 all]
measure sumweights $sel weight charge
# ^The charge was originally +5, now zero.
################################ 
# Check for loop penetration
# From makeEdn2.tcl
proc checkThread {molA resN loopName} {
  # Check loop penetration: 
  # molA=molid; resN=residue numbers; loopName=loop atom names
  set OK 0
  foreach p $resN {
    set others [atomselect $molA "beta 1 and not residue $p"]
    set this [atomselect $molA "residue $p and sidechain"]
    set loop [atomselect $molA "residue $p and name $loopName"]
    set loopCnt [measure center $loop]
    set oops [measure contacts 2.1 $others $this]
    set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
    if {$nCon>2} {
      foreach i [lsort -unique [lindex $oops 0]] {
        set left [atomselect $molA "index $i"]
        # Get distance of left atom to loopCnt
        set leftPt [$left get {x y z}]
        set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
        if {$d<1.0} {set OK 1}
        $left delete
      }
    }
    $others delete; $this delete; $loop delete
  }  
  return $OK 
}
# Set the CD3 CTs beta 1 to screen for threading
set all3 [atomselect $tcr3 all]; $all3 set beta 0
set eCT3 [atomselect $tcr3 "chain E and resid 156 to 207"]
set dCT3 [atomselect $tcr3 "chain D and resid 130 to 171"]
set gCT3 [atomselect $tcr3 "chain G and resid 139 to 182"]
set fCT3 [atomselect $tcr3 "chain F and resid 156 to 207"]
set nCT3 [atomselect $tcr3 "chain N and resid 309 to 312"]
set aCT3 [atomselect $tcr3 "chain A and resid 58 to 164"]
set bCT3 [atomselect $tcr3 "chain B and resid 56 to 164"]
$eCT3 set beta 1; $dCT3 set beta 1;$gCT3 set beta 1; $fCT3 set beta 1
$aCT3 set beta 1; $bCT3 set beta 1
# Get residue numbers. These are global residue numbers (Amber numbering), 
# not like resid, which are chain-specific.
set ProRes [atomselect $tcr3 "name CA and resname PRO and beta 1"]
set pro [$ProRes get residue]; $ProRes delete
# No Phe in CD3 CTs or in Nck SH3.1
set TyrRes [atomselect $tcr3 "name CA and resname TYR and beta 1"] 
set tyr [$TyrRes get residue]; $TyrRes delete
set HisRes [atomselect $tcr3 "name CA and resname HSD and beta 1"]
set his [$HisRes get residue]; $HisRes delete
set TrpRes [atomselect $tcr3 "name CA and resname TRP and beta 1"]
set trp [$TrpRes get residue]; $TrpRes delete
# Start checking
if {[checkThread $tcr3 $pro "N CA CB CD CG"]} {puts "Pro Loop penetrated"}
if {[checkThread $tcr3 $tyr "CZ CE1 CD1 CG CD2 CE2"]} {puts "Tyr Loop penetrated"}
if {[checkThread $tcr3 $his "CG ND1 CE1 NE2 CD2"]} {puts "His Loop penetrated"}
if {[checkThread $tcr3 $trp "CG CD1 NE1 CE2 CD2"]} {puts "Trp1 Loop penetrated"}
if {[checkThread $tcr3 $trp "CE2 CD2 CE3 CZ3 CH2 CZ2"]} {puts "Trp2 Loop penetrated"}
puts "Passed the loop penetration test?"
###########################################################
# Clean up
file delete water.pdb; file delete ions.pdb; file delete mem.pdb;
file delete aChain.pdb; file delete bChain.pdb; file delete dChain.pdb;
file delete eChain.pdb; file delete fChain.pdb; file delete gChain.pdb;
file delete mChain.pdb; file delete nChain.pdb; file delete CHL.pdb
file delete S1.pdb; file delete S2.pdb
foreach f [glob water*.pdb] {file delete $f}
# Get zero-based protein residue and index numbers. 
# Change to one-based, Amber numbers for restraining during initial MD.
lindex [[atomselect $tcr3 "name CA and protein"] get residue] end; # max zero-based residue=681 
lindex [[atomselect $tcr3 protein] get index] end; # max zero-based index=10870
# These protein residues and atoms are the same as tcr15c and tcr21a models.
lindex [[atomselect $tcr3 "not water and not ion"] get residue] end; # residues=2283
# ^ Two more residues (structural CHL1) from tcr15c and tcr21a models.
quit

