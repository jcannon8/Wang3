# anPRS3.tcl: Analyze PRS dihedral angles in edr and fgr (Fig 4 C,D,G,H).
cd /t/tcr16/eds
# Get backbone dihedral angles with VMD.
# vmd -dispdev none
# Wrap angles as decribed by Hollingsworth and Karplus 2010
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3061398/
# Wrapped Ramachandran plots reset phi<0 to phi+360 and psi< -120 to psi+360.
# Then plot -120 to 240 psi on y-axis, 0 to 360 phi on x-axis.
proc phiOut {phi outFile last} {
# Output one line of phi angles with wrapping.
  for {set j 0} {$j<=8} {incr j} {
    set a [lindex $phi $j]
    set b [expr $a<0?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
proc psiOut {psi outFile last} {
# Output one line of psi angles with wrapping.
  for {set j 0} {$j<=8} { incr j} {
    set a [lindex $psi $j]
    set b [expr $a<-120?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
# anEdr6.tcl used cpptraj to get the dihedral angles but did not wrap them.
set outPhi [open "edrfgr.phi.dat" w]
set outPsi [open "edrfgr.psi.dat" w]
# 11 edr references. Does not include 116 but does include 201.
set Edrmod {147 16 176 154 107 201 20 62 174 75 197}
foreach mod $Edrmod {
  set edr [mol new /t/tcr16/edm/edm.$mod.top type parm7 waitfor all]
  # Use the 215 ns frame
  mol addfile /t/tcr16/edr/edr.$mod.eq54.rst type netcdf waitfor all
  set phi [[atomselect $edr "name CA and resid 130 to 138"] get phi]
  phiOut $phi $outPhi 1
  set psi [[atomselect $edr "name CA and resid 130 to 138"] get psi]
  psiOut $psi $outPsi 1
  mol delete $edr
}
# 14 fgr references. Does not include 2 but does include 201.
set Fgrmod {105 177 199 93 10 83 168 142 30 55 154 117 86 201}
foreach mod $Fgrmod {
  set fgr [mol new /t/tcr16/fgm/fgm.$mod.top type parm7 waitfor all]
  # Use the 210 ns frame
  mol addfile /t/tcr16/fgr/fgr.$mod.eq53.rst type netcdf waitfor all
  set phi [[atomselect $fgr "name CA and resid 58 to 66"] get phi]
  phiOut $phi $outPhi 2
  set psi [[atomselect $fgr "name CA and resid 58 to 66"] get psi]
  psiOut $psi $outPsi 2
  mol delete $fgr
}
# Original NckPRS9803 reference
set 2jxb3 [mol new /home/jcannon/tcr2/cd3e/2jxb3.c0.pdb waitfor all]
set phi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get phi]
phiOut $phi $outPhi 3
set psi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get psi]
psiOut $psi $outPsi 3
#
close $outPhi
close $outPsi
###########################################################
# Plot angles for edr,fgr,edt,fgt (four plots).
# Fig 4 C,D,G,H 
gnuplot<<"eof"
# size is xscale, yscale
#set term jpeg font "arial.ttf,18" size 1280,960
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edr3.dih.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
# linetype controls the color of points rather than line style
set linetype 1 lc rgb "blue" 
set linetype 2 lc rgb "red" 
set linetype 3 lc rgb "green" 
set xrange [0:9.5]
set grid ytics
# offset points other than model 1 (edr)
jit(x,m)=m==1?x:x+0.1
# Do not label offset angle name
lab(s,m)=m==1?s:""
# plot-specific parameters for edr/fgr psi in uper left
root="edrfgr"
set title "{/:Bold Nck-bound CD3εδ^R and CD3εγ^R}\n{/:Bold MD}" offset -1,-1
set ylabel "{/:Bold ψ angle (degrees)}" offset 1,0
ext="psi"
set yrange [-120:240]
fname=sprintf("%s.%s.dat",root,ext)
# Each column has a different angle, each row has a different model.
# Tenth colum has 1-3 for point color
plot fname using (jit(1,$10)):1:10:xtic(lab("P180",$10)) with points lc var pt 7 t "",\
'' u (jit(2,$10)):2:10:xtic(lab("P181",$10)) with points lc var pt 7 t "",\
'' u (jit(3,$10)):3:10:xtic(lab("P182",$10)) with points lc var pt 7 t "",\
'' u (jit(4,$10)):4:10:xtic(lab("V183",$10)) with points lc var pt 7 t "",\
'' u (jit(5,$10)):5:10:xtic(lab("P184",$10)) with points lc var pt 7 t "",\
'' u (jit(6,$10)):6:10:xtic(lab("N185",$10)) with points lc var pt 7 t "",\
'' u (jit(7,$10)):7:10:xtic(lab("P186",$10)) with points lc var pt 7 t "",\
'' u (jit(8,$10)):8:10:xtic(lab("D187",$10)) with points lc var pt 7 t "",\
'' u (jit(9,$10)):9:10:xtic(lab("Y188",$10)) with points lc var pt 7 t ""
# Second plot upper right: edr/fgr phi
set ylabel "{/:Bold φ angle (degrees)}" offset 1,0
ext="phi"
set yrange [0:360]
fname=sprintf("%s.%s.dat",root,ext)
# repeat the last plot command with new arguments
replot
# Third plot lower left: edt/fgt psi
root="edtfgt"
set title "{/:Bold Nck-free CD3εδ and CD3εγ}\n{/:Bold MD}" offset 0,-1
set ylabel "{/:Bold ψ angle (degrees)}" offset 1,0
ext="psi"
set yrange [-120:240]
fname=sprintf("%s.%s.dat",root,ext)
replot
# Fourth plot lower right: edt/fgt phi
set ylabel "{/:Bold φ angle (degrees)}" offset 1,0
ext="phi"
set yrange [0:360]
fname=sprintf("%s.%s.dat",root,ext)
replot
#
quit
eof

