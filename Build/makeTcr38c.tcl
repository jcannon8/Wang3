# makeTcr38c.tcl: Test adding CD3 zeta CTs to tcr38.* to make tcr38.*.**.
# This reports on collisions from models with various zeta3 frames.
# Constructions are made without regard to Nck binding.
# Derived from makeTcr37a.tcl
# vmd -dispdev none on Eire
cd /t/tcr38b
set mod $env(mod);   # input / output tcr38.[0-8] model number
set aRef $env(aRef); # chain a zeta reference (0-8), actually use only 0-2
set bRef $env(bRef); # chain b zeta reference (0-8), actually use only 0-2
set inMod tcr38.$mod; # input tcr38.* parent model to add zeta chains
set outMod tcr38.$mod.$aRef$bRef; # output model
set nlast 30; # Number of models to consider
# Use 211.2 ns, eq45a autoimaged frame for tcr38.0, tcr38.3
# Use 200 ns, eq49a autoimaged frame for tcr38.2
set eq eq45
set eqa eq45a
set m [file exists /t/tcr38/$inMod.$eqa.rst]
if {$m==0} {
  exec /t/tcr37/image1.sh /t/tcr38/$inMod $eq
}
set tcr [mol new /t/tcr38/$inMod.psf]
mol addfile /t/tcr38/$inMod.$eqa.rst type netcdf waitfor all
# Do not consider Nck binding 
# Atomselections for CTs
set CTfa [atomselect $tcr "chain F and resid 157 to 207"]
set CTda [atomselect $tcr "chain D and resid 130 to 171"]
set CTea [atomselect $tcr "chain E and resid 156 to 207"]
set CTga [atomselect $tcr "chain G and resid 139 to 182"]
# Testing membrane collision of heavy atoms.
set lipids [atomselect $tcr "chain J"]
llength [lsort -unique [$lipids get resid]]
# ^ There are 1600 lipid residues in membrane
set sCHL [atomselect $tcr "chain S"] 
llength [lsort -unique [$sCHL get resid]]
# ^ There are 2 structural cholesterols
# Get average membrane z-dimension for inner leaflet.
set memDim [measure minmax $lipids]
set min [lindex $memDim 0 2]
set max [lindex $memDim 1 2]
set aveMem [expr ($min+$max)/2]
# Get inner phosphate z-dimension
set phos [atomselect $tcr "name P"]
set memDim2 [measure minmax $phos]
set maxP31 [lindex $memDim2 1 2]
# cytoplasmic leaflet heavy atoms
set mem [atomselect $tcr "chain J and z>$aveMem and mass>2"]
set conX 1.0;  # Collisions with membrane heavy atoms are <= 1.0 angstrom
set maxConM 3;  # maximum number of membrane contacts acceptable
set maxCon 5;  # Maximum tolerable contacts between CTs
# Finding membrane collisions is tedious. Save non-collision frames after first pass
# in files called chain[A,B].$ref.
# Also check that CTs do not stick out of box, particularly chain a and b in z-dimension.
set all [atomselect $tcr all]
set box [measure minmax $all]
set maxZ [lindex $box 1 2]
####################
# Add zeta3 CT to chain a. From addTail14c.tcl
set m [file exists chainA.$mod.$aRef]
if {$m==0} {set outA [open chainA.$mod.$aRef w]}
set cd3a [mol new /home/jcannon/tcr15/zeta.top type parm7]
mol addfile /home/jcannon/tcr15/zeta3.c$aRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/zeta3.c$aRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $cd3a get numframes]
# Cytoplasmic tail resid 6-112 of zeta3. New addition to parent are residues 2-112.
set CTa [atomselect $cd3a "resid 2 to 112 and mass>2"]
set CTaa [atomselect $cd3a "resid 6 to 112"]
#[atomselect $tcr "chain A and name N"] get {resname resid}
# CD3 epsilon: 6JXR chain a resid 22 to 57 in parent
# No CT atoms should go deeper into the membrane than the last N atom in 6JXR.
set aLim [[atomselect $tcr "chain A and name N and resid 57"] get z]
# Fit between tcr and zeta3 uses LYS PHE SER ARG residues in each.
# Align resid 2-5 of cd3ex and resid 54 to 57 of parent
set End [atomselect $tcr "chain A and name CA and resid 54 to 57"]
set Start [atomselect $cd3a "resid 2 to 5 and name CA"]
set cd3CD [atomselect $cd3a all]
for {set f 0 } {$f < $lastframe} {incr f} {
  $Start frame $f; $cd3CD frame $f
  $CTa frame $f; $CTaa frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
	if {$m==0} {
    # Determine membrane collision	
    set limit [measure minmax $CTaa]
    set min [lindex $limit 0 2]
    set max [lindex $limit 1 2]
    # Save list of frames that did not collide with membrane and stay within the box.
    if {$min<$aLim} {continue}
    if {$max>$maxZ} {continue}
    set con [measure contacts $conX $CTa $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxConM} {
      puts [format "%d frame, min=%5.2f, aLim=%5.2f" $f $min $aLim]
      lappend aFL $f
    }
  } 			
}
if {$m==0} {puts $outA $aFL; close $outA}
###########################################################
# Add zeta3 CT to chain b
set m [file exists chainB.$mod.$bRef]
if {$m==0} {set outB [open chainB.$mod.$bRef w]}
set cd3b [mol new /home/jcannon/tcr15/zeta.top type parm7]
mol addfile /home/jcannon/tcr15/zeta3.c$bRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/zeta3.c$bRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $cd3b get numframes]
# Cytoplasmic tail resid 4-112 of zeta3. New addition to parent are residues 2-112.
set CTb [atomselect $cd3b "resid 2 to 112 and mass>2"]
set CTba [atomselect $cd3b "resid 4 to 112"]
#[atomselect $tcr "chain B and name N"] get {resname resid}
# CD3 epsilon: 6JXR chain b resid 24 to 55 in parent
# No CT atoms should go deeper into the membrane than the last N atom in 6JXR.
set bLim [[atomselect $tcr "chain B and name N and resid 55"] get z]
# Fit between tcr and zeta3 uses VAL LYS PHE residues in each.
# In retrospect, zeta3 should have started with one more N-terminal residue so
# that its start was ARG VAL LYS PHE, instead of VAL LYS PHE.
# Increase to "CA C N", the atoms aligned because it is short.
# Align resid 1-3 of cd3b and resid 53 to 55 of parent
set End [atomselect $tcr "chain B and name CA C N and resid 53 to 55"]
set Start [atomselect $cd3b "name CA C N and resid 1 to 3"]
set cd3CD [atomselect $cd3b all]
for {set f 0 } {$f < $lastframe} {incr f} {
  $Start frame $f; $cd3CD frame $f
  $CTb frame $f; $CTba frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
	if {$m==0} {
    # Determine membrane collision and distance from membrane.	
    set limit [measure minmax $CTba]
    set min [lindex $limit 0 2]
    set max [lindex $limit 1 2]
    # Save list of frames that did not collide with membrane and stay within the box.
    if {$min<$bLim} {continue}
    if {$max>$maxZ} {continue}
    set con [measure contacts $conX $CTb $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxConM} {
      puts [format "%d frame, min=%5.2f, bLim=%5.2f" $f $min $bLim]
      lappend bFL $f
    }
  } 			
}
if {$m==0} {puts $outB $bFL; close $outB}
###########################################################
# Now test collisions CTs that do not collide with membrane.
# Load lists of frames that did not collide with membrane.
set inA [open chainA.$mod.$aRef r]
set aFL [read $inA]; close $inA
set inB [open chainB.$mod.$bRef r]
set bFL [read $inB]; close $inB
set aLen [llen $aFL]
set bLen [llen $bFL]
if {$aLen==0||$bLen==0} {
  puts "Cannot make any models because of collisions or box escape."
  puts "Possible ab frames = $aLen $bLen"
  exit
}
set n 0;		  # number of successful framesets, found
# nlast=requested number of framesets at begining
set i 0;		  # number of ab pairs tested, tries
set memAll [atomselect $tcr "chain J and z>$aveMem"]
# Log of potential constructions
set out [open $outMod.dat w]
while {$n<$nlast } {
	incr i
  if {$i>150000 || ($i>1000 && $n==0)} {
   puts $out "Model $outMod is not feasible"
   close $out
   quit
  }
  # Output tries and found framesets every 100 tries.
  if {[expr $i % 100]==0} {puts "$i $n"}
  # set a b to random frame that does not collide with membrane.
  set ia [expr round(rand()*($aLen-1))]; set a [lindex $aFL $ia] 
	set ib [expr round(rand()*($bLen-1))]; set b [lindex $bFL $ib] 
  $CTaa frame $a; $CTba frame $b 
  # Test f and a collision 
	set con [measure contacts $conX $CTfa $CTaa]
	set fa [llength [lindex $con 0]]
  if {$fa>$maxCon} {continue}
  # Test f and b collision 
	set con [measure contacts $conX $CTfa $CTba]
	set fb [llength [lindex $con 0]]
  if {$fb>$maxCon} {continue}
  # Test d and a collision
	set con [measure contacts $conX $CTda $CTaa]
	set da [llength [lindex $con 0]]
  if {$da>$maxCon} {continue}
  # Test d and b collision 
	set con [measure contacts $conX $CTda $CTba]
	set db [llength [lindex $con 0]]
  if {$db>$maxCon} {continue}
  # Test g and a collision 
	set con [measure contacts $conX $CTga $CTaa]
	set ga [llength [lindex $con 0]]
  if {$ga>$maxCon} {continue}
	# Test g and b collision  
	set con [measure contacts $conX $CTga $CTba]
	set gb [llength [lindex $con 0]]
  if {$gb>$maxCon} {continue}
	# Test e and a collision  
	set con [measure contacts $conX $CTea $CTaa]
	set ea [llength [lindex $con 0]]
  if {$ea>$maxCon} {continue}
	# Test e and b collision  
	set con [measure contacts $conX $CTea $CTba]
	set eb [llength [lindex $con 0]]
  if {$eb>$maxCon} {continue} 
	# Test a and b collision 
	set con [measure contacts $conX $CTaa $CTba]
	set ab [llength [lindex $con 0]]
  if {$ab>$maxCon} {continue}
	# If it reaches here, there were few collisions between CTs
  # A frame pair was found that met the collision criteria.       
  #
  # Now check all atom contacts with membrane
  # Above and saved in the files are only heavy atom membrane contacts.
  # Additional max all-atom contacts.
  set strict 1 
  set con [measure contacts $conX $CTaa $memAll]
 	set amem [llength [lindex $con 0]]
  if {$strict==1 && $amem >$maxCon} {continue}    
  set con [measure contacts $conX $CTba $memAll]
 	set bmem [llength [lindex $con 0]]
  if {$strict==1 && $bmem >$maxCon} {continue}
  incr n
  set memSum [expr $amem+$bmem]
  set icSum [expr $fa+$fb+$da+$db+$ga+$gb+$ea+$eb+$ab]
  set sum [expr $memSum+$icSum]
#  # Output: <Total contacts><a frame><b frame> | <amem><bmem> |
# fa fb da db ga gb ea eb ab 
  puts -nonewline $out [format "%3d%6d%6d |%3d%3d |" \
    $sum $a $b $amem $bmem]  
  puts $out " $fa $fb $da $db $ga $gb $ea $eb $ab"
  flush $out 
}
puts $out "Total tries: $i"
close $out
quit
