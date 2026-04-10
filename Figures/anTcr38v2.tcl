# anTcr38v2.tcl: Plot CD3 CT anchor configuration and membrane position (Fig 16).
# in tcr38b and tcr39b models.
# Single 405 ns (eq85) frame analyzed.
# vmd -dispdev none
cd /t/tcr39an
###########################################################
set m tcr38; # =tcr38 for tcr38b; =tcr39 for tcr39b
# Get shortest anchor distance to membrane within 20
proc getSmem {anch mem20} {
# parameters: atomselections anch for anchor atom and mem20 within 20 of anchor
global tcr
set cMem [measure contacts 20 $anch $mem20]
set cNum [llength [lindex $cMem 1]]
if {$cNum==0} {return 99}
unset -nocomplain dl     
foreach a2 [lindex $cMem 0] b2 [lindex $cMem 1] {
  set ax [atomselect $tcr "index $a2"]
  set bx [atomselect $tcr "index $b2"]
  set aa [$ax get {x y z}]
  set bb [$bx get {x y z}]
  set d [vecdist [lindex $aa 0] [lindex $bb 0]]; # distance between atom pair
  lappend dl $d
  set minMem [lindex [lsort $dl] 0]; # Shortest distance to membrane           
  $ax delete; $bx delete
}
return $minMem
}
#
set out1 [open ${m}v85a.dat w]; # anchor configuration
set out2 [open ${m}v85b.dat w]; # anchor membrane position using local Z dimension
set out3 [open ${m}v85c.dat w]; # anchor membrane position using membrane distance 
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fName $m.$j.$k 
# Skip missing  models. These are 405 ns frames.
if {[file exist /t/${m}b/$fName.psf]==0} continue
if {[file exist /t/${m}b/$fName.eq85.rst ]==0} continue
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
# Anchor points for six CD3 subunits, second CT residue CA.
set e157 [atomselect $tcr "name CA and chain E and resid 157"]
set f157 [atomselect $tcr "name CA and chain F and resid 157"]
set d131 [atomselect $tcr "name CA and chain D and resid 131"]
set g140 [atomselect $tcr "name CA and chain G and resid 140"]
set a59 [atomselect $tcr "name CA and chain A and resid 59"]
set b59 [atomselect $tcr "name CA and chain B and resid 59"]
# tcr38b has S2 DOPC
set s2P [atomselect $tcr "name P and chain S and resname DOPC"]
# Get anchor XY coordinates
# https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node163.html
set e157c [list [$e157 get x] [$e157 get y]]
set f157c [list [$f157 get x] [$f157 get y]]
set d131c [list [$d131 get x] [$d131 get y]]
set g140c [list [$g140 get x] [$g140 get y]]
set a59c [list [$a59 get x] [$a59 get y]] 
set b59c [list [$b59 get x] [$b59 get y]]
set s2Pc [list [$s2P get x] [$s2P get y]]
# Calculate anchor XY center
set cen [vecscale [expr 1.0/6] [vecadd $e157c $f157c $d131c $g140c $a59c $b59c]]
# anchor distances from center
set eDis [vecdist $e157c $cen]
set fDis [vecdist $f157c $cen]
set dDis [vecdist $d131c $cen]
set gDis [vecdist $g140c $cen]
set aDis [vecdist $a59c $cen]
set bDis [vecdist $b59c $cen]
# Get maximum distance for anchor circle radius, max.
set max 0; set min 100
# Get anchor nearest to center, near.
foreach d "$eDis $fDis $dDis $gDis $aDis $bDis" chn {E F D G A B} {
puts "$d $chn"
if {$d > $max} {set max $d}
if {$d < $min} {set min $d; set near $chn}
}
# Output anchor configuration: <model><distances><max radius><center anchor>
puts $out1 [format "%6s%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f | %7.2f %4s" \
 $fName $eDis $fDis $dDis $gDis $aDis $bDis $max $near]
#
# Get anchor membrane postions like dimer.mem.tcl
# Using atomselections initially to get membrane heavy atoms within 20 angstroms.
# Then get the minimum Z-dimension distance.
set mem [atomselect $tcr "mass>2 and chain J"]; # Same for all anchors
set memZ [lindex [measure minmax $mem] 1 2]; # This is the maximum Z-dimension of membrane.
# Chain E
set deltaE [expr $memZ - [$e157 get z]];  # E anchor distance to total membrane.
# ^ deltaE>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near Ala157 in chain E
set memE2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain E and resid 157)"]
set memEZ2 [lindex [measure minmax $memE2] 1 2] 
# ^ Maximum membrane Z-dimension near E chain anchor
set deltaE2 [expr $memEZ2 - [$e157 get z]]; # E anchor distance to local membrane
set deltaE3 [getSmem $e157 $memE2]; # Shortest E anchor to membrane.
# Chain F
set deltaF [expr $memZ - [$f157 get z]];  # F anchor distance to total membrane.
# ^ deltaF>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near Ala157 in chain F
set memF2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain F and resid 157)"]
set memFZ2 [lindex [measure minmax $memF2] 1 2] 
# ^ Maximum membrane Z-dimension near F chain anchor
set deltaF2 [expr $memFZ2 - [$f157 get z]]; # E anchor distance to local membrane
set deltaF3 [getSmem $f157 $memF2]; # Shortest F anchor to membrane.
# Chain D
set deltaD [expr $memZ - [$d131 get z]];  # D anchor distance to total membrane.
# ^ deltaD>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near Gly131 in chain D
set memD2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain D and resid 131)"]
set memDZ2 [lindex [measure minmax $memD2] 1 2] 
# ^ Maximum membrane Z-dimension near D chain anchor
set deltaD2 [expr $memDZ2 - [$d131 get z]]; # D anchor distance to local membrane
set deltaD3 [getSmem $d131 $memD2]; # Shortest D anchor to membrane.
# Chain G
set deltaG [expr $memZ - [$g140 get z]];  # G anchor distance to total membrane.
# ^ deltaG>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near D140 in chain G
set memG2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain G and resid 140)"]
set memGZ2 [lindex [measure minmax $memG2] 1 2] 
# ^ Maximum membrane Z-dimension near G chain anchor
set deltaG2 [expr $memGZ2 - [$g140 get z]]; # G anchor distance to local membrane
set deltaG3 [getSmem $g140 $memG2]; # Shortest G anchor to membrane.
# Chain A
set deltaA [expr $memZ - [$a59 get z]];  # A anchor distance to total membrane.
# ^ deltaA>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near Ala59 in chain A
set memA2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain A and resid 59)"]
set memAZ2 [lindex [measure minmax $memA2] 1 2] 
# ^ Maximum membrane Z-dimension near A chain anchor
set deltaA2 [expr $memAZ2 - [$a59 get z]]; # A anchor distance to local membrane
set deltaA3 [getSmem $a59 $memA2]; # Shortest A anchor to membrane.
# Chain B
set deltaB [expr $memZ - [$b59 get z]];  # B anchor distance to total membrane.
# ^ deltaB>0 if inserted into membrane cytoplasmic side leaflet
# Get membrane near Ala59 in chain B
set memB2 [atomselect $tcr "mass>2 and chain J and within 20 of (name CA and chain B and resid 59)"]
set memBZ2 [lindex [measure minmax $memB2] 1 2] 
# ^ Maximum membrane Z-dimension near B chain anchor
set deltaB2 [expr $memBZ2 - [$b59 get z]]; # B anchor distance to local membrane
set deltaB3 [getSmem $b59 $memB2]; # Shortest B anchor to membrane.
# S2 DOPC phosphate
set deltaS 0; set deltaS2 0; # defaults for tcr39b
if {$m == "tcr38"} {
  set deltaS [expr $memZ - [$s2P get z]];  # DOPC S2 distance to total membrane.
  # ^ deltaS>0 if inserted into membrane cytoplasmic side leaflet
  # Get membrane near DOPC S2
  set memS2 [atomselect $tcr "mass>2 and chain J and within 20 of (name P and chain S and resname DOPC)"]
  set memSZ2 [lindex [measure minmax $memS2] 1 2] 
  # ^ Maximum membrane Z-dimension near DOPC S2
  set deltaS2 [expr $memSZ2 - [$s2P get z]]; # DOPC S2 distance to local membrane
  set deltaS3 [getSmem $s2P $memS2]; # Shortest DOPC S2 distance to membrane.
}
# Output anchor distance to total membrane and distance to local membrane for each anchor.
# Output: <model><E total><F total> ... <E local ><F local> ...
puts $out2 [format "%6s%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f   %7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f" \
 $fName $deltaE $deltaF $deltaD $deltaG $deltaA $deltaB $deltaS \
 $deltaE2 $deltaF2 $deltaD2 $deltaG2 $deltaA2 $deltaB2 $deltaS2]
#
puts $out3 [format "%6s%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f   %7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f" \
 $fName $deltaE3 $deltaF3 $deltaD3 $deltaG3 $deltaA3 $deltaB3 $deltaS \
 $deltaE2 $deltaF2 $deltaD2 $deltaG2 $deltaA2 $deltaB2 $deltaS2]
# Delete atomselections
$e157 delete;$f157 delete;$d131 delete;$g140 delete;$a59 delete;$b59 delete
$memE2 delete;$memF2 delete;$memD2 delete;$memG2 delete;$memA2 delete;$memB2 delete
if {$m == "tcr38"} {$s2P delete;$memS2 delete}
mol delete $tcr
}
}
close $out1;close $out2;close $out3
#
###########################################################
# boxplot graphs CD3 CT anchor data
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'tcr38v1.jpg'
set title "{/:Bold Anchor circle radii }"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 8
set xtic nomirror 
set ytic out nomirror
set xrange [0.5:12.5]
set grid ytics
set ylabel "{/:Bold Distance (Å)}" offset 1,0
set style boxplot nooutliers pointtype 7
set style data boxplot
set boxwidth 0.5
#set key textcolor variable samplen -1
set key horizontal textcolor variable samplen -1
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 7 pi -1 ps 1.0
# 
plot 'tcr39v85a.dat' u (1):9 with boxplot ls 1 t "Tcr39b resting",\
'tcr38v85a.dat' u (2):9 with boxplot ls 2 t "Tcr38b GOF"
eof
# The range of radii in Tc38b GOF is greater than Tcr39b.
###########################################################
# boxplot graphs CD3 CT2 membrane positions and distances
# Fig 17A,B
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'tcr38v2.jpg'
# layout: rows, columns
set multiplot layout 1,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 5
set rmargin 3
set xtic out nomirror 
set ytic out nomirror offset 1,0
set yrange [-15:18]
set grid ytics
set ylabel "{/:Bold Distance (Å)}" offset 2.5,0
set style boxplot nooutliers pointtype 7
set style data boxplot
set boxwidth 0.5
# Reduce key font size so that horizontal will work.
set key left horizontal textcolor variable samplen -1 font "arial.ttf,16" 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set xtics ("εE" 1.5, "εF" 4.5, "δD" 7.5, "γG" 10.5, "ζA" 13.5, "ζB" 16.5)
#
set title "{/:Bold CT2 cytoplasm positions ❷ }"
plot 'tcr39v85b.dat' u (1):9 with boxplot ls 1 t "",\
'tcr38v85b.dat' u (2):9 with boxplot ls 2 t "",\
'tcr39v85b.dat' u (4):10 with boxplot ls 1 t "",\
'tcr38v85b.dat' u (5):10 with boxplot ls 2 t "",\
'tcr39v85b.dat' u (7):11 with boxplot ls 1 t "",\
'tcr38v85b.dat' u (8):11 with boxplot ls 2 t "",\
'tcr39v85b.dat' u (10):12 with boxplot ls 1 t "",\
'tcr38v85b.dat' u (11):12 with boxplot ls 2 t "",\
'tcr39v85b.dat' u (13):13 with boxplot ls 1 t "",\
'tcr38v85b.dat' u (14):13 with boxplot ls 2 t "",\
'tcr39v85b.dat' u (16):14 with boxplot ls 1 t "TCR",\
'tcr38v85b.dat' u (17):14 with boxplot ls 2 t "TCR-GOF"
#
set title "{/:Bold CT2 Membrane distances ❸}"
# Showing outliers is better here.
set style boxplot outliers pointtype 7
set ytic out nomirror add ("" 9.5,"" 10.5,"" 11.5,"" 12.5) offset 1,0 
set yrange [9:13]
set ylabel "{/:Bold Distance (Å)}" offset 5,0
plot 'tcr39v85c.dat' u (1):2 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (2):2 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (4):3 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (5):3 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (7):4 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (8):4 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (10):5 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (11):5 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (13):6 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (14):6 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (16):7 with boxplot ls 1 t "TCR",\
'tcr38v85c.dat' u (17):7 with boxplot ls 2 t "TCR-GOF"
eof
###########################################################
# boxplot graphs CD3 CT2 membrane distances
gnuplot<<eof
set term jpeg font "arial.ttf,18" 
set out 'tcr38v3.jpg'
set title "{/:Bold CT2 Membrane distances }"
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 8
set xtic nomirror 
set ytic out nomirror
set yrange [9:13]
set grid ytics
set ylabel "{/:Bold Distance (Å)}" offset 1,0
#set style boxplot nooutliers pointtype 7
# Showing outliers is better here.
set style boxplot pointtype 7
set style data boxplot
set boxwidth 0.5
# Reduce key font size so that horizontal will work.
set key left horizontal textcolor variable samplen -1 font "arial.ttf,16" 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 7 pi -1 ps 1.0
set xtics ("εE" 1.5, "εF" 4.5, "δD" 7.5, "γG" 10.5, "ζA" 13.5, "ζB" 16.5)
# 
plot 'tcr39v85c.dat' u (1):2 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (2):2 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (4):3 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (5):3 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (7):4 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (8):4 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (10):5 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (11):5 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (13):6 with boxplot ls 1 t "",\
'tcr38v85c.dat' u (14):6 with boxplot ls 2 t "",\
'tcr39v85c.dat' u (16):7 with boxplot ls 1 t "Tcr39b resting",\
'tcr38v85c.dat' u (17):7 with boxplot ls 2 t "Tcr38b GOF"
eof

###########################################################
# Use R to get t-tests to Compare circle radii
read.table("tcr38v85a.dat")->t38
read.table("tcr39v85a.dat")->t39
t.test(t39[,9],t38[,9])
# ^ p-value = 0.2113. 
# Anchor (CT2) Circle radii difference is insignificant.
#########################
# Compare CT2 membrane positions
read.table("tcr38v85b.dat")->t38
read.table("tcr39v85b.dat")->t39
t.test(t39[,9],t38[,9])
p <- numeric()
m38 <- numeric()
m39 <- numeric()
dis <- numeric()
for(i in seq(9,14,1)) {
  # two-sided t-test of columns 9-14
  t.test(t39[,i],t38[,i])->q  
  p[i-8] <- format(q$p.value, digits =3)
  m39[i-8] <- round(mean(t39[, i]),digits=3) 
  m38[i-8] <- round(mean(t38[, i]),digits=3)
  dis[i-8] <- round((mean(t38[, i]) - mean(t39[, i])),digits=3) 
}
n<-c("E","F","D","G","A","B")
d<-data.frame(n,m39,m38,p,dis)
write.table(d,"tcr38v85d.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
read.table("tcr38v85c.dat")->t38
read.table("tcr39v85c.dat")->t39
t.test(t39[,9],t38[,9])
p <- numeric()
m38 <- numeric()
m39 <- numeric()
dis <- numeric()
for(i in seq(2,7,1)) {
  # two-sided t-test of columns 2-7
  t.test(t39[,i],t38[,i])->q  
  p[i-1] <- format(q$p.value, digits =3)
  m39[i-1] <- round(mean(t39[, i]),digits=3) 
  m38[i-1] <- round(mean(t38[, i]),digits=3)
  dis[i-1] <- round((mean(t38[, i]) - mean(t39[, i])),digits=3) 
}
n<-c("E","F","D","G","A","B")
d<-data.frame(n,m39,m38,p,dis)
write.table(d,"tcr38v85e.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#####################
# CT2 cytoplasm positions
cat tcr38v85d.dat
E 0.515 0.83 0.59 0.315
F 1.481 2.66 0.0368 1.179 *
D 1.36 3.16 0.00237 1.8 *
G 3.465 3.955 0.391 0.489
A -0.839 0.236 0.038 1.074 *
B -2.603 -1.827 0.319 0.777
###########################################################
#
