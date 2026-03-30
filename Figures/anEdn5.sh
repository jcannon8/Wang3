# anEdn5.sh: PRS dihedral cluster of 200-member edn and fgn ensembles (S5 Fig).
# Script anEdn1.tcl made protein-only edn.prot.pdb and ednEns.dcd
# RunTest.sh has dihedral clustering examples. 
# /usr/share/amber18/AmberTools/src/cpptraj/test/Test_Cluster_DataSet/RunTest.sh 
###########################################################
# Two ways of dihedral clustering here.
# Edn ensemble clustering without intermediate dihedral files.
cpptraj<<eof
parm /t/tcr16/edn.prot.pdb
trajin /t/tcr16/ednEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :132@N :132@CA :132@C :133@N
dihedral P184psi :134@N :134@CA :134@C :135@N
dihedral P186psi :136@N :136@CA :136@C :137@N
# Make COORDS set crd1
createcrd crd1
# It is unclear why a COORDS set is needed for the clustering because it is the dihedral data
# that is clustered, which is irrelevant to coordinates at this point.
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out anEdn5.out \
summary anEdn5.summary.dat info anEdn5.info.dat 
eof
###########################################################
# Edn ensemble clustering using dihedrals via readdata
# First write dihedrals
cpptraj<<eof
parm /t/tcr16/edn.prot.pdb
trajin /t/tcr16/ednEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :132@N :132@CA :132@C :133@N
dihedral P184psi :134@N :134@CA :134@C :135@N
dihedral P186psi :136@N :136@CA :136@C :137@N
# This "go" is necessary before the write command.
go
write eP182psi.dat P182psi
write eP184psi.dat P184psi
write eP186psi.dat P186psi
eof
# Second read dihedrals
cpptraj<<eof
parm /t/tcr16/edn.prot.pdb
# Read previously calculated dihedrals
readdata eP182psi.dat name P182psi
readdata eP184psi.dat name P184psi 
readdata eP186psi.dat name P186psi
# So it knows above data is periodic, very important!
# http://archive.ambermd.org/201902/0251.html
dataset mode torsion P182psi
dataset mode torsion P184psi
dataset mode torsion P186psi
# Make COORDS set crd1
createcrd crd1
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out anEdn5b.out \
summary anEdn5b.summary.dat info anEdn5b.info.dat 
eof
# anEdn5.* and anEdn5b.* output is the same, so I know I did the two-step correctly.
###########################################################
# Combine cluster assignment and Nck binding energy into temp2 for edn.
# Energies per model are in col 29, 30 of edn.AvsE4.dat from anEdn4.tcl.
echo "ELECT    VDW" > temp
awk '{printf("%7.3f %7.3f\n", $29,$30)}' edn.AvsE4.dat >> temp
# Cluster# per model are in col 2 of anEdn5.out
paste anEdn5b.out temp > temp2
###########################################################
# Fgn ensemble clustering without intermediate dihedral files.
cpptraj<<eof
parm /t/tcr16/fgn.prot.pdb
trajin /t/tcr16/fgnEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :60@N :60@CA :60@C :61@N
dihedral P184psi :62@N :62@CA :62@C :63@N
dihedral P186psi :64@N :64@CA :64@C :65@N
# Make COORDS set crd1
createcrd crd1
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out fgn5.out \
summary fgn5.summary.dat info fgn5.info.dat 
eof
# Combine cluster assignment and Nck binding energy into temp2 for fgn.
echo "ELECT    VDW" > temp
awk '{printf("%7.3f %7.3f\n", $29,$30)}' fgn.AvsE4.dat >> temp
# Cluster# per model are in col 2 of anEdn5.out
paste fgn5.out temp > temp2
###########################################################
# Save fgn dihedrals
cpptraj<<eof
parm /t/tcr16/fgn.prot.pdb
trajin /t/tcr16/fgnEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :60@N :60@CA :60@C :61@N
dihedral P184psi :62@N :62@CA :62@C :63@N
dihedral P186psi :64@N :64@CA :64@C :65@N
# This "go" is necessary before the write command.
go
write fP182psi.dat P182psi
write fP184psi.dat P184psi
write fP186psi.dat P186psi
eof
###########################################################
# Combine edn and fgn dihedals
cat eP182psi.dat > cP182psi.dat
sed '1d' fP182psi.dat >> cP182psi.dat
cat eP184psi.dat > cP184psi.dat
sed '1d' fP184psi.dat >> cP184psi.dat
cat eP186psi.dat > cP186psi.dat
sed '1d' fP186psi.dat >> cP186psi.dat
###########################################################
# Cluster edn and fgn together
cpptraj<<eof
parm /t/tcr16/edn.prot.pdb
# Read previously calculated dihedrals
readdata cP182psi.dat name P182psi
readdata cP184psi.dat name P184psi 
readdata cP186psi.dat name P186psi
# So it knows above data is periodic.
# http://archive.ambermd.org/201902/0251.html
dataset mode torsion P182psi
dataset mode torsion P184psi
dataset mode torsion P186psi
# Make COORDS set crd1
createcrd crd1
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out anEdn5c.out \
summary anEdn5c.summary.dat info anEdn5c.info.dat 
eof
# Combine cluster assignment and Nck binding energy into temp2 for combined edn and fgn.
echo "ELECT    VDW" > temp
awk '{printf("%7.3f %7.3f\n", $29,$30)}' edn.AvsE4.dat >> temp
awk '{printf("%7.3f %7.3f\n", $29,$30)}' fgn.AvsE4.dat >> temp
# Cluster# per model are in col 2 of anEdn5c.out
paste anEdn5c.out temp > temp2
awk '{printf("%7.3f %7.3f\n", $29,$30)}' edn.AvsE4.dat > temp0
paste anEdn5c.out temp0 > edn.AvsE5.dat
awk '{printf("%7.3f %7.3f\n", $29,$30)}' fgn.AvsE4.dat > temp0
paste anEdn5c.out temp0 > fgn.AvsE5.dat
###########################################################
# Plot cluster results
gnuplot<<"eof"
set term jpeg font "arial.ttf,18"
# Function for linear interaction energy (alpha=0.18, beta=0.5)
g(e,v) = 0.18 * v + 0.5 * e
#set out 'edn.AvsE3.jpg'
#set out 'fgn.AvsE3.jpg'
set out 'edn.AvsE4.jpg'
set title "{/:Bold CD3εδ^N and CD3εγ^N 5 ns PRS dihedral clustering}"
set border 3
set xtic in nomirror 0,1,4
set xrange [-0.5:5]
set ytic in nomirror
set xzeroaxis
set xlabel "{/:Bold Cluster}"
set ylabel "{/:Bold Nck binding energy}"
set linetype 1 lc rgb "blue" lw 2 pt 7
set linetype 2 lc rgb "red" lw 2 pt 7
set linetype 3 lc rgb "green" lw 2 pt 7
# The combined temp2 below depends on which clustering to be plotted.
plot 'temp2' using 2:(g($3,$4)) with points lt 1 t "" 
quit
eof
###########################################################
# Plot cluster results: Energy vs cluster
# S5 Fig panel C
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1500, 560
# Function for linear interaction energy (alpha=0.18, beta=0.5)
g(e,v) = 0.18 * v + 0.5 * e
set out 'ang.clust.jpg'
# layout: rows, columns
set multiplot layout 1,2
#set title "{/:Bold Edn and Fgn 5 ns PRS dihedral clustering}"
set border 3
set xtic in nomirror 0,1,4
set xrange [-0.5:5]
set ytic in nomirror
set xzeroaxis
set xlabel "{/:Bold Cluster}"
set ylabel "{/:Bold Nck interaction energy}\n{/:Bold (kcal/mole)} "
set linetype 1 lc rgb "blue" lw 2 pt 7
set linetype 2 lc rgb "red" lw 2 pt 7
set linetype 3 lc rgb "green" lw 2 pt 7
fname="edn.AvsE5.dat"
set title "{/:Bold CD3εδ^N PRS dihedral clustering}"
# The combined temp2 below depends on which clustering to be plotted.
plot fname using 2:(g($3,$4)) with points lt 1 t ""
#
fname="fgn.AvsE5.dat"
set title "{/:Bold CD3εγ^N PRS dihedral clustering}"
replot 
quit
eof
###########################################################
# Get angles for centroids
# anEdn5.summary.dat has the centroid frame in column 5
# eP182psi, eP184psi.dat have the angles in column 2 for frames listed in column 1
paste eP182psi.dat eP184psi.dat anEdn5.out > eAng2.dat
paste fP182psi.dat fP184psi.dat fgn5.out > fAng2.dat
###########################################################
# Plot angles in color-coded by cluster 
gnuplot<<"eof"
set term jpeg font "arial.ttf,18"
set out 'edn.ang2.jpg'
#set out 'fgn.ang2.jpg'
set title "{/:Bold Edn angles}"
#set title "{/:Bold Fgn angles}"
set border 3
set xtic in nomirror 
set xrange [-100:200]
set yrange [-100:200]
set ytic in nomirror
set xzeroaxis
set yzeroaxis
set xlabel "{/:Bold P182ψ}"
set ylabel "{/:Bold P184ψ}"
# For the "lc var" to work in the plot command, use "linetypes", not "style line".
# Colors: https://i.stack.imgur.com/x6yLm.png
set linetype 1 lc rgb "blue" lw 3 pt 7 pi -1 ps 1.0
set linetype 2 lc rgb "red" lw 3 pt 7 pi -1 ps 1.0
set linetype 3 lc rgb "forest-green" lw 3 pt 7 pi -1 ps 1.0
set linetype 4 lc rgb "purple" lw 3 pt 7 pi -1 ps 1.0
set linetype 5 lc rgb "dark-cyan" lw 3 pt 7 pi -1 ps 1.0
# P182psi in column 2, P184psi in column 4, cluster# (0-4) in column 6
plot 'eAng2.dat' using 2:4:($6+1) with points ls 1 lc var t ""
#plot 'fAng2.dat' using 2:4:($6+1) with points ls 1 lc var t ""
quit
eof
###########################################################
# Use vmd to render cluster centroid examples.
# Fig 5 panel E
set tcr [mol new /t/tcr16/edn.prot.pdb]
mol addfile /t/tcr16/ednEns.dcd waitfor all
# PRS: chain E and resid 180 to 188
display projection Orthographic
display depthcue off
color Display Background white
axes location Off
# Set up representation
mol delrep 0 $tcr 
mol representation Licorice 0.20 12.0  12.0
mol color Name
mol selection mass>2 and chain E F and resid 180 to 188
mol addrep $tcr
# Set the display to the desired size. 
# Ten peptides across 2250 pixels wide -> 225 pixel width each.
display resize 225 400
# Change frame to cluster centroid of chain E
# Edn centroids: 173 137 62 28 54
animate goto 173
animate goto 137
animate goto 62
animate goto 28
animate goto 54
# Manually, move PRS to center of display window with Pro180 at top.
# Copy to clipboard with Win+Shift+S, paste in Illustrator.
#
set tcr [mol new /t/tcr16/fgn.prot.pdb]
mol addfile /t/tcr16/fgnEns.dcd waitfor all
display projection Orthographic
display depthcue off
color Display Background white
axes location Off
# Set up representation
mol delrep 0 $tcr 
mol representation Licorice 0.20 12.0  12.0
mol color Name
mol selection mass>2 and chain E F and resid 180 to 188
mol addrep $tcr
# Set the display to the desired size. 
# Ten peptides across 2250 pixels wide -> 225 pixel width each.
display resize 225 400
# Change frame to cluster centroid of chain F
# Fgn centroids: 63 171 93 55 94
animate goto 63
animate goto 171
animate goto 93
animate goto 55
animate goto 94
# Manually, move PRS to center of display window with Pro180 at top.
# Copy to clipboard with Win+Shift+S, paste in Illustrator.
