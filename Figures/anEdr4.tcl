# anEdr4.tcl: Analyze Nck-PRS interactions in edr (S3 Fig, Fig 3A).
# This makes the four plots for S3 Fig.
cd /t/tcr16/eds
# vmd -dispdev none on Eire.
# Get the one-based residue numbers for edn/edr/edt
set tcr [mol new /t/tcr16/edn.prot.pdb waitfor all]
[atomselect $tcr "chain E and name CA"] get {resname resid}
# PRS is Pro180-Tyr188
[atomselect $tcr "chain E and name CA and resid 180 to 188"] get {resname resid residue}
# ^That shows PRS are residues 129-137 (zero-based), 130-138 (one-based)
[atomselect $tcr "chain E and name CA"] get {resname resid residue}
# ^That shows chain E are residues 74-156 (zero-based), 75-147 (one-based)
[atomselect $tcr "chain D and name CA"] get {resname resid residue}
# ^That shows chain D are residues 0-73 (zero-based), 1-74 (one-based)
###########################################################
# Analyze Nck and PRS RMSD, LIE, and distances for edr models.
# edr reference models: 154 16 75 174 147 20 176 116 197 62 107 201
for mod in 154 16 75 174 147 20 176 116 197 62 107 201 ; do
name=edr.$mod
echo $name
rm -f trajSum
# Get 425 ns trajectories.
for ((i=13; i<=96;i++)); do
if [ -e /t/tcr16/edr/edr.$mod.eq$i.cdf ]; then
echo "trajin /t/tcr16/edr/edr.$mod.eq$i.cdf" >>trajSum
fi
done
# Use cpptraj to get three outputs: $name.[rms,lie,dist].dat
cpptraj<<eof
# Use the edm topology
parm /t/tcr16/edm/edm.$mod.top
# Only 3405 protein atoms in trajectory
parmstrip !@1-3405
readinput trajSum
# Get RMS
# See https://ambermd.org/tutorials/analysis/tutorial1/index.php
# for rmsd tutorial to see how to use references.
# Nck reference is frame 9803 of 2jxb3 MD.
parm /t/tcr16/eds/Nck9803.pdb [ref0]
# Use ref0 topology and name reference "ref2".
reference /t/tcr16/eds/Nck9803.pdb parm [ref0] [ref2]
# PRS reference is also frame 9803
parm /t/tcr16/eds/NckPRS9803.pdb [ref3]
reference /t/tcr16/eds/NckPRS9803.pdb parm [ref3] [ref4]
# Note order of masks: first mask is for trajectory, second for reference.
rmsd Nck ref [ref2] ^3@CA @CA out $name.rms.dat
# PRS Amber residue numbers (one-based) determined via VMD above
rmsd PRS ref [ref4] :130-138@CA :6-14@CA out $name.rms.dat
# Linear interaction energy between Nck and chain E (PRS or all) and chain D
# Nck (chain K) residues 158-213 (one-based)
# CD3e PRS residues 130-138 (one-based), CD3e residues 75-157 (one-based)
# CD3d (chain D) residues 1-74 (one-based)
lie PRS :130-138 :158-213 out $name.lie.dat
lie CD3e :75-157 :158-213 out $name.lie.dat
lie CD3d :1-74 :158-213 out $name.lie.dat
# Restraint distances see anNck3.sh
distance PN1 :131@O :209@HH out $name.dist.dat
distance PN2 :135@O :192@HE1 out $name.dist.dat
distance PN3 :130@CD,CA,CB,CG :192@CE2,CD2,CE3,CZ3,CZ2,CH2 out $name.dist.dat
distance PN4 :138@HH :173@OE1 out $name.dist.dat
# This next one has atoms in two residues
distance PN5 :138@CG,CD1,CE1,CZ,CD2,CE2 \
:192@CE2,CD2,CE3,CZ3,CZ2,CH2,:204,CG,CD1,CE1,CZ,CD2,CE2 out $name.dist.dat
# PRS Val183 carbonyl to Asn185 amide
distance HB12 :133@O :135@HN out $name.dist.dat
# PRS Asp187 to Nck Lys36
distance salt :137@CG :164@NZ out $name.dist.dat
eof
done
rm -f trajSum
###########################################################
# Figure S3 Fig panel A
# Plot all Nck LIE for edr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr.lie.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^R Nck-CD3εδ interaction energy}"
set border 3
set rmargin 1.5
set lmargin 4.3
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out -250,50,0 nomirror offset 1,0
set yrange [-250:10]
set grid ytics
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
g(e,v) = 0.18 * v + 0.33 * e
#
L=0.45
do for [n in "154 16 75 174 147 20 176 116 197 62 107 201" ] {
  set title sprintf("edr.%s",n) offset 0,-1
  fname = sprintf("edr.%s.lie.dat",n)
  # Titles in screen coordinates.
  plot fname u (\$1/100) : (g(\$2,\$3)) with lines ls 2 t "Nck-PRS" at L,0.19,\
  '' u (\$1/100) : (g(\$4,\$5)) with lines ls 1 t "Nck-CD3ε" at L, 0.16,\
  '' u (\$1/100) : (g(\$6,\$7)) with lines ls 3 t "Nck-CD3εδ" at L, 0.13
}
quit
eof
###########################################################
# S3 Fig panel B
# Plot edr Nck-PRS distances
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr4.dis.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^R Nck-PRS distances}"
set border 3
set rmargin 1.5
set lmargin 3
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out nomirror 0,10,40 offset 1,0
set yrange [0:40]
set key Left reverse textcolor variable samplen -1 font "arial.ttf,18"
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
L=0.50
do for [n in "154 16 75 174 147 20 176 116 197 62 107 201" ] {
  set title sprintf("edr.%s",n) offset 0,-1
  filename = sprintf("edr.%s.dist.dat",n)
  plot filename using (\$1/100) : 2 with lines ls 1 title \
  "PRS Pro181 O – Nck Tyr82 HH" at L,0.25,\
  '' using (\$1/100) : 3 with lines ls 2 title \
  "PRS Asn185 O – Nck Trp65 HE1" at L, 0.22,\
  '' using (\$1/100) : 4 with lines ls 3 title \
  "PRS Pro180 - Nck Trp65" at L, 0.19,\
  '' using (\$1/100) : 6 with lines ls 4 title \
  "PRS Tyr188 - Nck Trp65,Tyr77" at L, 0.16,\
  '' using (\$1/100) : 8 with lines ls 5 title \
  "PRS Asp187 - Nck Lys37" at L, 0.13
}
quit
eof
###########################################################
# S3 Fig panel C
# Plot edr RMSD
models="154 16 75 174 147 20 176 116 197 62 107 201"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr3.rms.jpg'
# layout: rows, columns
set multiplot layout 3,5 title "{/:Bold CD3εδ^R Nck and PRS RMSD}"
set border 3
# margin units are character heights or widths
set rmargin 1.5
set lmargin 3
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out nomirror 0,1,4 offset 1,0
set yrange [0:4.5]
set key left horizontal textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
L=0.50
do for [n in "$models"] {
  set title sprintf("edr.%s",n) offset 0,-1
  fname = sprintf("edr.%s.rms.dat",n)
  plot fname using (\$1/100):2 with lines ls 1 t "Nck" at L,0.16,\
  '' using (\$1/100):3 with lines ls 2 t "PRS" at L,0.19
}
quit
eof
###########################################################
# Check CD3 PRS dihedral angles for Nck binding criteria in edr
# Using wrapped angles, these limits prevent Nck binding
# P182 ψ <100°, P184 ψ <100° 
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
set edrModels "154 16 75 174 147 20 176 116 197 62 107 201"
foreach mod $edrModels {
set out [open edrPRS16.$mod.dat w]
# Use stripped topology with just 3405 protein atoms that were saved in trajectory.
set tcr [mol new /t/tcr16/edn.prot.psf]
#[atomselect $tcr "name CA"] get {chain resname resid}
# ^To check resid
# The NPT MD started with edr.*.eq13, 425 ns total.
for {set i 13} {$i<=100} {incr i} {
  set mdFile /t/tcr16/edr/edr.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
# Use edr resid from /t/tcr16/edn.prot.psf
set P182 [atomselect $tcr "chain E and name CA and resid 182"]
set V183 [atomselect $tcr "chain E and name CA and resid 183"]
set P184 [atomselect $tcr "chain E and name CA and resid 184"]
set lastframe [molinfo $tcr get numframes]
for {set f 0} {$f<$lastframe} {incr f} {
  set badA 0; set bad182 0; set bad183 0; set bad184 0
  $P182 frame $f; $V183 frame $f; $P184 frame $f
  # Use the 7-27-24 criteria
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 100} {incr badA; incr bad182}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 100} {incr badA; incr bad184}
  # Output: <frame><total bad angles><bad 182><bad184>
  puts $out [format "%6d%4d%4d%4d" $f $badA $bad182 $bad184]
  # Could output raw or wrapped angles as well.
}
close $out
mol delete $tcr
}
# 7-18-24 angle criteria
if {0} {
  set P182psi [$P182 get psi] 
  if {[wrapPsi $P182psi] < 50} {incr badA; incr bad182}
  set P183phi [$V183 get phi]
  if {[wrapPhi $P183phi] < 240} {incr badA; incr bad183}
  set P184psi [$P184 get psi]
  if {[wrapPsi $P184psi] < 50} {incr badA; incr bad184}
}
###########################################################
# S3 Fig panel C
# Plot PRS angle violations for edr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr.ang.jpg'
set multiplot layout 3,5 title "{/:Bold CD3εδ^R PRS angle violations}"
set border 3
set rmargin 1.5
set lmargin 4.3
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out 0,1,2 nomirror offset 1,0
set yrange [-0.5:2.5]
set grid ytics
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
#
do for [n in "154 16 75 174 147 20 176 116 197 62 107 201" ] {
  set title sprintf("edr.%s",n) offset 0,-1
  fname = sprintf("edrPRS16.%s.dat",n)
  plot fname u (\$1/100):2 with points ls 2 t ""
}
quit
eof
###########################################################
# Plot PRS angle which violation angles for edr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr2.ang.jpg'
set multiplot layout 3,5 title "{/:Bold Edr which PRS angle violations}"
set border 3
set rmargin 1.5
set lmargin 4.3
set xtic 0,200,440 out nomirror
# Minor tics
set mxtics 2
set xrange [0:440]
set ytic out 0,1,3 nomirror offset 1,0
set yrange [-0.5:3.5]
set grid ytics
# These all have thin lines.
set style line 1 lc rgb "red" lt 1 lw 1 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 1 pt 1 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 1 pt 7 pi -1 ps 1.0
#
do for [n in "154 16 75 174 147 20 176 116 197 62 107 201" ] {
  set title sprintf("edr.%s",n) offset 0,-1
  fname = sprintf("edrPRS16.%s.dat",n)
  plot fname u (\$1/100):(\$3==1?1:1/0) with points ls 2 t "", \
  '' u (\$1/100):(\$4==1?2:1/0) with points ls 2 t "", \
  '' u (\$1/100):(\$5==1?3:1/0) with points ls 2 t ""
}
quit
eof

###########################################################
cat edrPRS16.*.dat > edrPRS16.sum.dat
# use R to get summary of violations
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
read.table("edrPRS16.sum.dat")->a
length(a[,1])
[1] 483000
> sum(a[,3])
[1] 4260 0.88%
> sum(a[,4])
[1] 10370 2.1%
> sum(a[,5])
[1] 7390 1.5%
###########################################################
# LIE distribution for edr ensemble
# Concatenate data from all models in ensemble.
# Concatenate only data from models with for CD3d[EELEC] <-10.
cat edr.*.lie.dat> edrEns.lie.dat
# Concatenate only data from models with CD3δ energies <-10.
awk '{if($6<-10 && $7<-10) print $0}' edr.*.lie.dat>edrEns2.lie.dat
# Fig 3A
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 500,375
set out 'edr2c.lie.jpg'
#set title "{/:Bold CD3εδ^R Nck-CD3εδ interaction energy}\n{/:Bold distribution}"
set border 3
set xtics out nomirror 
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
# Use linestyle here
set linestyle 1 lc rgb "red" lw 2 
set linestyle 2 lc rgb "blue" lw 2
set linestyle 3 lc rgb "dark-green" lw 2
set style fill transparent pattern 0 noborder
set xlabel "{/:Bold Interaction energy (kcal/mol)}"
set ylabel "{/:Bold Relative frequency}" offset 2,0
# Function for binning
bin(x,s) = s*int(x/s)
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
g(e,v) = 0.18 * v + 0.33 * e
# Plots with smooth, frequency, with lines
# y-axis is relative frequency 
# with r=0.01, max y-axis was 250
# with r=0.001, ,ax y-axis is 25
r=0.001
plot "edrEns.lie.dat" u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 2 t "Nck-PRS",\
'' u (bin((g(\$4,\$5)),1.0)):(r) s f w lines ls 1 t "Nck-CD3ε",\
"edrEns2.lie.dat" u (bin((g(\$6,\$7)),1.0)):(r) s f w lines ls 3 t "Nck-CD3δ"
quit
eof
# The y-axis needs some help because the units of frequency are ambiguous.
#


