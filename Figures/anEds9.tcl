# anEds9.tcl: Nck binding to eds (S7 fig panel A).
# Check Nck binding to eds 200 ns MD using NckPRS.dcd refs and 1.0 threshold.
# and combine with collision data to calculate Nck binding.
# Use the 7-27-24 criteria: P182 ψ <100°, P184 ψ <100°
# vmd -dispdev none on Eire.
cd /t/tcr16/eds
# Contacts are heavy atoms less than or equal to 1.0 angstroms. 
set dist 1.0;  # Nck collision distance for heavy atoms
# Output: edsPRS20.*.dat
##########################################
# Load the Nck-PRS references made by makeNckref.tcl
set nckRef [mol new /t/tcr39an/NckPRS.pdb]
mol addfile /t/tcr39an/NckPRS.dcd waitfor all
# Define reference atomselections
set PRSref [atomselect $nckRef "name CA and chain E and resid 180 to 186"]
set nckSH [atomselect $nckRef "chain K"]
set refAll [atomselect $nckRef all]
set numRef [expr [molinfo $nckRef get numframes]-1]
###########################################################
# eds models are edr models with Nck removed.
set edsModels "154 16 75 174 147 20 176 197 62 107 201"
# Process each eds model
foreach mod $edsModels {
# Use stripped topology with just 2448 protein atoms that were saved in trajectory,
# made by anEdt1.tcl
# Load trajectories
set tcr [mol new /t/tcr16/edt.prot.psf]
# The eds 106 ns MD eq2-24
for {set i 2} {$i<=24} {incr i} {
  set mdFile /t/tcr16/eds/eds.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
set out [open edsPRS20.$mod.dat w]
# PRS is Pro180-Pro186 in 6JXR chains E and F (see anEdt2.sh).
# Use resid here because edt.prot.psf supplies chain names and can use chain resid.
set PRS [atomselect $tcr "name CA and chain E and resid 180 to 186"]
# notPRS atoms are more than two residues away from PRS
set notPRS [atomselect $tcr "not hydrogen and not (chain E and resid 178 to 188)"] 
# CA of first chain E CT residue Ala157 used to judge membrane location. 
set Mem [atomselect $tcr "name CA and chain E and resid 157"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  $PRS frame $f
  $notPRS frame $f
  $Mem frame $f
  # Next Nck contacts for each nckRef
  for {set ref 0} {$ref<$numRef} {incr ref} {
    $PRSref frame $ref; $nckSH frame $ref; $refAll frame $ref;
    # defaults in case no Nck binding.
    set RMSD 0; set nckMem -200; set nRef -1
    # Move refAll to bind to ePRS
    set M [measure fit $PRSref $PRS]; # Fit the 7 PRS CA atoms
    $refAll move $M
    # Check contacts with Nck SH3.1
    set nckCon [measure contacts $dist $nckSH $notPRS] 
    set conNum [llen [lindex $nckCon 1]]
    if {$conNum==0} {
      set nRef $ref
      set RMSD [measure rmsd $PRSref $PRS]
      # Save Nck to membrane distance
      set loc [measure minmax $nckSH]
      set minZ [lindex $loc 0 2]
      set memZ [$Mem get z]
      set eNckMem [expr $minZ- $memZ]
      break; # Stop searching once the first Nck binding without contact found.
    }
  }
  puts [format "%6d%4d %7.2f %7.2f" $f $nRef $RMSD $eNckMem]
  puts $out [format "%6d%4d %7.2f %7.2f" $f $nRef $RMSD $eNckMem]
  # output <frame><eNck binding><eNck membrane>
  # Nck collision with membrane if MemNck < -6.0 (check on this).
}
close $out
mol delete $tcr
}
###########################################################
# Analyze frequency of Nck binding
# This uses the NckPRS.dcd Nck-PRS references above.
tclsh<<"eof"
# Disregard the PRS RMSD
set models {154 16 75 174 147 20 176 197 62 107 201} 
foreach mod $models {
puts $mod
set in [open edsPRS20.$mod.dat r]
set inData [read $in]
close $in
set out [open edsPRS21.$mod.dat w]
# Example output of edtPRS20.$mod.dat.
#   0   0    2.33   21.55
#   1   0    2.60   21.52
set nCon 0; set mCon 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  incr frame
  set tm [lindex $line 0]
  incr tm
  # Get Nck-PRS ref index
  set nRef [lindex $line 1]
  incr r($nRef); # track how many times each Nck-PRS ref used.
  set RMSD [lindex $line 2]
  set dMem [lindex $line 3]
  # Nck can bind if no contacts with other proteins and 
  # no Nck contact with membrane.
  if {$nRef>-1 && $dMem>-6.0 } {incr nCon}
  # If no protein contact 
  if {$nRef>-1 } {incr mCon}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # output <time (ns)><percent Nck bound><percent no protein contact>
    puts $out [format "%4d %4.3f %4.3f" [expr $tm/100] [expr $nCon/10.0] [expr $mCon/10.0]]
    set nCon 0; set mCon 0
  }
}
close $out
}
# Output how many times each Nck reference bound.
parray r
puts "Total frames=$frame"
eof
###########################################################
# Plot Nck binding frequency for eleven eds models in multiplot
# edt models incluster prevalence order:
models="154 16 75 174 147 20 176 197 62 107 201"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'eds.nck8.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold Potential Nck binding to eds CD3εδ}"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 6
set xtics out nomirror 0,25,110
set ytic out nomirror offset 1,0
set xrange [0:110]
set yrange [0:100]
set key bottom center textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("eds.%s",n) offset 0,-1
  if (n==147) {set ylabel "{/:Bold Potential Nck binding frequency}" offset 3.4,0}
  else {set ylabel " "}
#  if (n>56) {set xlabel "{/:Bold Time (ns)}"}
#  else {unset xlabel}
  name2 = sprintf("edsPRS21.%s.dat",n)
  plot name2 using 1 : 2 with linespoints ls 2 t "pro + mem",\
  '' using 1 : 3 with linespoints ls 1 t "pro only"
}
quit
eof
###########################################################
# Check CD3 PRS dihedral angles for Nck binding criteria
# vmd -dispdev none
# Using wrapped angles, these limits prevent Nck binding
# P182 ψ <100°, P184 ψ <100°
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
set edsModels "154 16 75 174 147 20 176 197 62 107 201"
# Process each eds model
foreach mod $edsModels {
# Use stripped topology with just 2448 protein atoms that were saved in trajectory,
# made by anEdt1.tcl
# Load trajectories
set tcr [mol new /t/tcr16/edt.prot.psf]
# The eds 106 ns MD eq2-24
for {set i 2} {$i<=24} {incr i} {
  set mdFile /t/tcr16/eds/eds.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
set out [open edsPRS24.$mod.dat w]
# Use edt resid from edt.prot.psf, which is different from edt.*.top
set P182 [atomselect $tcr "chain E and name CA and resid 182"]
set P184 [atomselect $tcr "chain E and name CA and resid 184"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  set badA 0; set bad182 0; set bad184 0
  $P182 frame $f; $P184 frame $f
  # Use the 7-27-24 criteria: P182 ψ <100°, P184 ψ <100°
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 100} {incr badA; incr bad182}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 100} {incr badA; incr bad184}
  # Output: <frame><total bad angles><bad 182><bad184>
  puts $out [format "%6d%4d%4d%4d" $f $badA $bad182 $bad184]
}
close $out
mol delete $tcr
}
###########################################################
# Combine collision in edsPRS20.*.dat with bad angles in edsPRS24.*.dat
for mod in 154 16 75 174 147 20 176 197 62 107 201; do
paste edsPRS20.$mod.dat edsPRS24.$mod.dat >temp
# Save <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad 183><bad184>
awk '{printf("%6d%4d%7.2f%7.2f%4d%4d%4d\n",$1,$2,$3,$4,$6,$7,$8)}' \
temp>edsPRS24.$mod.dat 
done
rm -f temp
# Now edsPRS24.$mod.dat has
# <frame><nRef><RMSD fit><membrane dist><bad angles><bad 182><bad 183><bad184>
###########################################################
# Analysis for frequency of Nck binding in 10 ns blocks.
# This reports individual criteria as well as cumulative criteria. 
tclsh<<"eof"
set models {154 16 75 174 147 20 176 197 62 107 201} 
foreach mod $models {
puts $mod
set in [open edsPRS24.$mod.dat r]
set inData [read $in]
close $in
set out [open edsPRS25.$mod.dat w]
# Example output of edsPRS24.$mod.dat. 
#      0   0   2.33  21.55   1   0   1   0
set nBind 0; set pCon 0; set mCon 0; set aCon 0
foreach line [split $inData "\n"] {
  if {$line == ""} break
  incr frame
  set tm [lindex $line 0]
  incr tm
  # Get Nck-PRS ref index
  set nRef [lindex $line 1]; #  Nck-PRS ref index 
  incr r($nRef); # track how many times each Nck-PRS ref used.
  set RMSD [lindex $line 2]; # Disregard the PRS RMSD
  set dMem [lindex $line 3]; # Distance to membrane
  set badA [lindex $line 4]; # number of bad PRS angles
  set bad182 [lindex $line 5]; # bad 182
  set bad184 [lindex $line 6]; # bad 184
  # No protein contact exlcuding fgr.201
  if {$nRef>-1 && $nRef<24} {incr pCon}
  # No membrane contact
  if {$dMem>-6.0} {incr mCon}
  # Just look at P182 and P184
  if {$bad182==0 && $bad184==0} {incr aCon}
  # Nck can bind if no contacts with other proteins, 
  # no Nck contact with membrane, and no bad angles.
  if {$nRef>-1 && $nRef<24 && $dMem>-6.0 && $bad182==0 && $bad184==0} {incr nBind}
  # In 10 ns blocks, each frame is 10 ps.
  if {[expr $tm % 1000]==0} {
    # Output in percents 
    # <time (ns)><no protein contact><no mem contacts><angles OK><binding possible>
    puts $out [format "%4d %4.3f %4.3f %4.3f %4.3f" \
        [expr $tm/100] [expr $pCon/10.0] [expr $mCon/10.0] \
        [expr $aCon/10.0] [expr $nBind/10.0]]
    set nBind 0; set pCon 0; set mCon 0; set aCon 0
  }
}
close $out
}
# Output how many times each Nck reference bound.
# This was useful in the ordering of references in anEdt7.tcl.
parray r
puts "Total frames=$frame"
eof
###########################################################
cd /t/tcr16/eds
# Plot Nck binding frequency for eleven eds models in multiplot
# Figure S7A
# edt models in cluster prevalence order:
models="154 16 75 174 147 20 176 197 62 107 201"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'eds.nck9.jpg'
# layout: rows, columns
#set multiplot layout 3,4 title "{/:Bold Potential Nck binding to CD3εδ}"
set multiplot layout 3,4 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set xtics out nomirror 0,50,220
set ytic out nomirror offset 1,0
set xrange [0:110]
set yrange [0:100]
# Legend placed in lower right 
set key Left reverse textcolor variable samplen -1 font "arial.ttf,16"
# These have a vairety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
#
L=0.75
do for [n in "$models"] {
  set title sprintf("eds.%s",n) offset 0,-1
  fname = sprintf("edsPRS25.%s.dat",n)
  # Titles in screen coordinates. Last points shifted to expose overlaps.   
  plot fname using 1:2 with linespoints ls 1 t "No protein collision" at L,0.19,\
  '' using 1:3 with linespoints ls 2 t "No membrane collision" at L, 0.16,\
  '' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at L, 0.13,\
  '' using (column(1)+2):5 with linespoints ls 4 t "Potential Nck binding" at L, 0.10
}
quit
eof
###########################################################
# Use R to consolidate Nck binding for all models of eds.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# edsPRS25.*.dat has
# <time (ns)><no protein contact><no mem contacts><angles OK><binding possible>
# I want the above in edtSum.dat with the averages over all 11 models. 
models=c("154","16","75","174","147","20","176","197","62","107","201")
# Matrices of times in rows, models in columns
P<-matrix(,nrow=10,ncol=11)
M<-matrix(,nrow=10,ncol=11)
A<-matrix(,nrow=10,ncol=11)
T<-matrix(,nrow=10,ncol=11)
# Fill above matrices with frequencies for each criteria.
for(i in seq(1,11,1)) {
  read.table(sprintf("edsPRS25.%s.dat",models[i]))->a
  P[,i] <- a[,2]
  M[,i] <- a[,3]
  A[,i] <- a[,4]
  T[,i] <- a[,5]
}
# Output has average over all models for binding criteria for each row.
out<-matrix(,nrow=10,ncol=5) 
for(t in seq(1,10,1)) {
  out[t,1]= t*10
  out[t,2] <- round(mean(P[t,]),2)
  out[t,3] <- round(mean(M[t,]),2)
  out[t,4] <- round(mean(A[t,]),2)
  out[t,5] <- round(mean(T[t,]),2)
}
write.table(out,file="edsSum.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
# Plot eds Nck binding averages,
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'eds.nck9b.jpg'
set title "{/:Bold Potential Nck binding to CD3εδ in eds}"
set border 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set xrange [0:110]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These have a vairety of point types and high contrast colors.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 2 pt 13 pi -1 ps 1.0
plot '/t/tcr16/eds/edsSum.dat' u 1:2 with linespoints ls 1 t "No protein collision" at 0.25,0.29,\
'' using 1:3 with linespoints ls 2 t "No membrane collision" at 0.25,0.26,\
'' using 1:4 with linespoints ls 3 t "Favorable PRS angles" at 0.25,0.23,\
'' using 1:5 with linespoints ls 4 t "Potential Nck binding" at 0.25,0.20 
quit
eof
