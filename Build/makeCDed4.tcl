# makeCDed4.tcl: Make CD3 epsilon-delta dimers with random CTs fused to truncated subunits.
# The models made are called edt.*. They have protein, membrane, and water.
# This is derived from makeCDed.tcl, which saved only protein to make CD3edt3.
# NOTE: This fixes a water error in makeCDed3.tcl
# Use vmd -dispdev none 
# Environment variables for model range desired
set first $env(first)
set last $env(last)
cd /home/jcannon/tcr16/edt
set out [open log.dat w]
# Use psfgen to build edt.$mod.[psf,pdb]
# https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf
# http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node6.html
package require psfgen
# Use 6.875 ns equilibrated CD3edt4.eq7 for reference. 
# Autoimage using chain E anchor to bring the most water into cytoplasm.
set tcr2 [mol new /home/jcannon/tcr16/charm-gui3/amber/CD3edt4.psf]
mol addfile /home/jcannon/tcr16/charm-gui3/amber/CD3edt4.ref2.pdb waitfor all
# That psf marks chains by segname. chain D: segname PROA, chain E: PROB
# Save TM coordinates of reference, tcr2, so that chains can be returned to original position.
set tmE0 [atomselect $tcr2 "name CA and segname PROB and resid 125 to 154"]
set tmD0 [atomselect $tcr2 "name CA and segname PROA and resid 98 to 129"]
# CT residues
set eCT [atomselect $tcr2 "segname PROB and resid 155 to 207"]
set dCT [atomselect $tcr2 "segname PROA and resid 130 to 171"]
# Add back waters that current chain E and D CTs contact because those chains are removed.
# removeExt.tcl built a big water box to fill in removed protein.
# Use method in removeExt2.tcl
set wat [mol new /home/jcannon/tcr2/tcr10/bigWater3.pdb]
# Make sure eCT and dCT are within water box, otherwise move the box
set water [atomselect $wat all]
set limWat [measure minmax $water]
set limECT [measure minmax $eCT]
set limDCT [measure minmax $dCT]
# The water is fine
$water set beta 0;  # Mark all waters to not add, beta=0
# Add cytoplasmic water that collides with old chainD and E CT. 
set conD [measure contacts 1.0 $dCT $water]
set wats [lindex $conD 1]
set selD [atomselect $wat "same residue as index $wats"]
$selD set beta 2
set conE [measure contacts 1.0 $eCT $water]
set wats [lindex $conE 1]
set selE [atomselect $wat "same residue as index $wats"]
$selE set beta 3
set wat1 [atomselect $wat "beta 2 3"]
# Water diagnostics
set initWat [[atomselect $tcr2 "name O and water"] num]
set wat1Num [[atomselect $wat "beta 2 3 and name OH2"] num]
puts $out "Initial waters=$initWat, wat1 contact waters=$wat1Num"
# Later, add waters with beta 2 3 add to final output if they do not colide with new CTs.
# Save lipids and ions
set lipids [atomselect $tcr2 "segname MEMB"]
$lipids set chain J
$lipids writepdb lipids.pdb
set ions [atomselect $tcr2 "name SOD CLA"]
$ions writepdb ions.pdb
##########
# A second CD3edt3 copy to change. This is the "top" molecule because it is the last.
set tcr [mol new /home/jcannon/tcr16/charm-gui3/amber/CD3edt4.psf]
mol addfile /home/jcannon/tcr16/charm-gui3/amber/CD3edt4.ref2.pdb waitfor all
[atomselect $tcr all] move [trans z by 180]
set chainD [atomselect $tcr "segname PROA"] 
set chainE [atomselect $tcr "segname PROB"]
$chainD set chain D
$chainE set chain E
set water [atomselect $tcr water]
set mem [atomselect $tcr "segname MEMB"] 
# Transmembrane residues
set tmE1 [atomselect $tcr "name CA and segname PROB and resid 125 to 154"]
set tmD1 [atomselect $tcr "name CA and segname PROA and resid 98 to 129"]
# CT residues
set eCT [atomselect $tcr "segname PROB and resid 155 to 207"]
set dCT [atomselect $tcr "segname PROA and resid 130 to 171"]
set eCT2 [atomselect $tcr "segname PROB and resid 155 to 207 and mass>2"]
set dCT2 [atomselect $tcr "segname PROA and resid 130 to 171 and mass>2"]
set PRSphi [atomselect $tcr "name CA and segname PROB and resid 180 to 187"]
set PRSpsi [atomselect $tcr "name CA and segname PROB and resid 180 to 188"]
# The atomselections are moved out of loop so they are not created more than once.
#
proc setAng {chn ang res} {
# chn=segname; ang=HBGLE; res=resid.
# Use HGBLE angles defined in http://eire.umh.edu/mediawiki/index.php/MakeTorsions
# phi, psi H: -57, -47; G: -49, -26; B: -119, 113; L 57, 47; E: -150, 155
  if {$ang == "H"} {
    [atomselect top "name CA and segname $chn and resid $res"] set phi -57
    [atomselect top "name CA and segname $chn and resid $res"] set psi -47
  }
  if {$ang == "G"} {
    [atomselect top "name CA and segname $chn and resid $res"] set phi -49
    [atomselect top "name CA and segname $chn and resid $res"] set psi -26
  }
  if {$ang == "B"} {
    [atomselect top "name CA and segname $chn and resid $res"] set phi -119
    [atomselect top "name CA and segname $chn and resid $res"] set psi 133
  }
  if {$ang == "L"} {
    [atomselect top "name CA and segname $chn and resid $res"] set phi 57
    [atomselect top "name CA and segname $chn and resid $res"] set psi 47
  }
  if {$ang == "E"} {
    [atomselect top "name CA and segname $chn and resid $res"] set phi -150
    [atomselect top "name CA and segname $chn and resid $res"] set psi 155
  }
}
# Start of loop
for {set mod $first; set t 1} {$mod<=$last} {incr t} {
if {[expr $t % 10]==0} {puts "$t $mod" }
#
# The makeHGBLE program makes a random HGBLE string using the CT sequence. 
# and the table of frequencies I have used before.
# HGBLEstring=`makeHGBLE eps.seq HGBLEfreq2 ` was used in MakeEpsilon.sh
# Karplus PA. Experimentally observed conformation-dependent geometry and 
# hidden strain in proteins. Protein Sci. 1996 Jul;5(7):1406-20. 
# doi: 10.1002/pro.5560050719. PMID: 8819173; PMCID: PMC2143451.
set strE [exec makeHGBLE /home/jcannon/tcr2/cd3f/eps.seq /home/jcannon/tcr2/cd3f/HGBLEfreq2]
# That ^ CD3 epsilon started at Trp 151. I only need Lys156 to Ile207
# Changing phi and psi in the CT will move the whole subunit.
# Use strE to set new phi, psi in CD3 chain E.
for {set r 156} {$r<=207} {incr r} {
  set j [expr $r-151]
  set ang [string range $strE $j $j]
  if {$r<180 || $r>188} {setAng "PROB" $ang $r}
}
# Set PRS (Pro180-Tyr188) angles based on centroid of most populated 2jxb3 MD cluster.
# Same angles in eps3.RST.
$PRSphi set phi {-73 -65 -44 -87 -66 -53 -77 -82}
$PRSpsi set psi {155 155 159 145 169 124 5 92 -15}
# Check for chain E self-contacts
set conE [measure contacts 1.0 $eCT2 $eCT2]
set nConE [llength [lindex $conE 0]]
if {$nConE>0} {continue}; # No heavy atom contacts allowed!
# Move chain E back to original location using reference, tcr2, TM CA coordinates.
set M [measure fit $tmE1 $tmE0]
$chainE move $M 
# chainE cannot extend into cytoplasm more than 122 using CD3edt4 derived from cd3ed.37.pdb.
set limit [measure minmax $chainE]
set deltaZ [expr abs([lindex $limit 0 2] - [lindex $limit 1 2])]
if {$deltaZ>122} {continue}
# Check chain E collision with membrane
set con [measure contacts 1.0 $eCT $mem]
set flCon [llength [lindex $con 0]]
if {$flCon >0} {continue}
puts "found a chain E"
#####
# Make CD3 delta CT similarly. Try 100 times with a discovered chain E
set notfound 1
for {set k 0} {$k<100 && $notfound} {incr k} {
set strD [exec makeHGBLE /home/jcannon/tcr2/cd3d/delta.seq /home/jcannon/tcr2/cd3f/HGBLEfreq2]
# That ^ CT delta starts at Ala126. Last 6JXR delta residue Glu129.
# I need Thr130 to Lys171 
for {set r 130} {$r<=171} {incr r} {
  set j [expr $r-126]
  set ang [string range $strD $j $j]
  setAng "PROA" $ang $r
}
# Check for chain D self-contacts
set conD [measure contacts 1.0 $dCT2 $dCT2]
set nConD [llength [lindex $conD 0]]
if {$nConD>0} {incr t; continue}; # No heavy atom contacts allowed!
# Move chain D back to original location using reference, tcr2, TM CA coordinates.
set M [measure fit $tmD1 $tmD0]
$chainD move $M 
# chainD cannot extend into cytoplasm more than 122.
set limit [measure minmax $chainD]
set deltaZ [expr abs([lindex $limit 0 2] - [lindex $limit 1 2])]
if {$deltaZ>122} {incr t; continue}
# Check chain D collision with membrane
set con [measure contacts 1.0 $dCT $mem]
set dlCon [llength [lindex $con 0]]
if {$dlCon >0} {incr t; continue}
puts "found a chain D"
set notfound 0
}
# End of k loop for chain D
if {$notfound == 1} {puts "try another chain E"; continue}
#
# Finally check collison of eCT and dCT
set con [measure contacts 1.0 $eCT $dCT]
set edCon [llength [lindex $con 0]]
if {$edCon >0} {incr t; continue}
###
# In makeCDed.tcl, new chain D E protein coordinates were saved together.
# In this script, save separately for psfgen, remove water new CTs contact,
# and add back water old CTs contacted.
$chainD writepdb chainD.pdb
$chainE writepdb chainE.pdb
# Remove cytoplasmic water that collides with new chainD and E CT. 
# Set beta=1 for water to keep.
$water set beta 1
set conD [measure contacts 1.0 $dCT $water]
set wats [lindex $conD 1]
set selD [atomselect $tcr "same residue as index $wats"]
$selD set beta 8
set conE [measure contacts 1.0 $eCT $water]
set wats [lindex $conE 1]
set selE [atomselect $tcr "same residue as index $wats"]
$selE set beta 7
# Water diagnostics
set conNew [[atomselect $tcr "beta 7 8 and name O OH2"] num]
puts $out "Subtract waters contacting new chain ED = $conNew"
$selE delete; $selD delete; # Delete atom selections to conserve memory.
#
# Now save all waters in mol tcr with beta=1
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr "water and beta 1"]
$saveWat set chain z
set Res [lsort -unique -integer [$saveWat get residue]]
set lastWat [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$lastWat} {incr n} {
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
  $toSave delete
}
$saveWat delete
# water1.pdb to water[$n-1].pdb saved above
# Save waters added from removed old CTs, wat1, which are not in contact with new CTs.
set conD2 [measure contacts 1.0 $dCT $wat1]
set conE2 [measure contacts 1.0 $eCT $wat1]
# Indices of wat2 waters to not add
# FIXED an error in makeCDed3.tcl
unset -nocomplain notWat2
append notWat2 [lindex $conE2 1] [lindex $conD2 1]
set addWat [atomselect $wat "beta 2 3 and not same residue as index $notWat2"]
# Water diagnostics
set addWatNum [expr [$addWat num] / 3]
puts $out "Add waters from initial CTs = $addWatNum"
flush $out
# Now save addWat
$addWat writepdb waterExt.pdb 
$addWat delete
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
# The above reports many errors but they are OK.
# Load segments with coordinates.
segment d {pdb chainD.pdb}; coordpdb chainD.pdb d
segment e {pdb chainE.pdb}; coordpdb chainE.pdb e
segment j {pdb lipids.pdb}; coordpdb lipids.pdb j
segment k {pdb ions.pdb}; coordpdb ions.pdb k
# Load in water
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
# Water is in multiple pdb files
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
# Add extra water
incr i
segment w$i {pdb waterExt.pdb}; coordpdb waterExt.pdb w$i
# Only one disulfide remains from 6jxr.pdb, between zeta chains, not here.
writepsf edt.$mod.psf
writepdb edt.$mod.pdb
# Important to reset the psf for next model.
resetpsf
set numAS [llength [atomselect list]]
incr mod
}
# End of loop
###########################################################
# Clean up
file delete waterExt.pdb; file delete ions.pdb; file delete lipids.pdb;
file delete chainD.pdb; file delete chainE.pdb
foreach f [glob water*.pdb] {file delete $f}
# I used 1.0 angstrom as the distance to consider atoms colliding. This might
# be too short because Van der Waals radii of H 1.0-1.4, C 1.9, O 1.66, N 1.8
# are much larger and contacts would occur at the sum of VdW radii. Nevertheless,
# contacts are discovered in equilibrated systems with the 1.0 angstrom distance.
close $out
