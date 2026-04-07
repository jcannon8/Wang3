# anTcr38ss.sh: Analyze secondary structure of CT residues in tcr38.*.** MD (S12 fig).
cd /t/tcr38an
# Use VMD to get 1-based Amber residue numbers for CD3 CT residues. 
vmd -dispdev none
set tcr [mol new /t/tcr38b/tcr38.0.00.psf]
mol addfile /t/tcr38b/tcr38.0.00.eq45.rst type netcdf waitfor all
set chnA2 [atomselect $tcr "chain A and resid 111 to 126  and name CA"]
$chnA2 get {resname resid residue}
# ITAM2 helix 90-105, 1-based residues
set chnB2 [atomselect $tcr "chain B and resid 111 to 126  and name CA"]
$chnB2 get {resname resid residue}
# ITAM2 helix 231-246, 1-based residues
set chnA3 [atomselect $tcr "chain A and resid 142 to 156 and name CA"]
$chnA3 get {resname resid residue}
# ITAM3 helix 121-135, 1-based residues
$chnA3 get helix  
# 1-based residues 1-143
set chnB3 [atomselect $tcr "chain B and resid 142 to 156 and name CA"]
$chnB3 get {resname resid residue}
# ITAM3 helix 262-276, 1-based residues
# ^^^ The 1-based residue numbers for tcr38b are the same as those for tcr39b.
###########################################################
# Use cpptraj to get secondary structure every 5 ns.
mem=":685-2284"
ncParm="noimage mindist"
res2="Y N E L Q K D K M A E A Y S E I"
res3="Y Q G L S T A T K D T Y D A L"
# Yes, res2 has 15 residues and res3 has 14 residues.
rm -f ssA2.h.dat ssB2.h.dat ssA3.h.dat ssB3.h.dat
rm -f ssA2.dat ssB2.dat ssA3.dat ssB3.dat
# Print residue name in first column
for r in $res2; do
echo $r>>ssA2.dat
echo $r>>ssB2.dat
done
for r in $res3; do
echo $r>>ssA3.dat
echo $r>>ssB3.dat
done
# Process all tcr38b models.
for ((i=0;i<=11;i++)); do
for j in 00 01 02 10 11 12 20 21 22; do
name=tcr38.$i.$j
p=/t/tcr38b/$name.top
if [ -e $p ]; then
# Valid tcr38b model, collect data from 206 ns (eq45), 405 ns (eq85)
# I will call this the 400-425 ns interval.
#
cpptraj<<eof
parm $p
# last 50 ns (375-425 ns)
trajin /t/tcr38b/$name.eq76.rst
trajin /t/tcr38b/$name.eq77.rst
trajin /t/tcr38b/$name.eq78.rst
trajin /t/tcr38b/$name.eq79.rst
trajin /t/tcr38b/$name.eq80.rst
# Last 25 ns (385-405 ns)
trajin /t/tcr38b/$name.eq81.rst
trajin /t/tcr38b/$name.eq82.rst
trajin /t/tcr38b/$name.eq83.rst
trajin /t/tcr38b/$name.eq84.rst
trajin /t/tcr38b/$name.eq85.rst
# Backbone amide hydrogen is HN, complete CD3ζ
secstruct :1-143   nameh HN out ssA 
secstruct :144-284 nameh HN out ssB
# Get distance to membrane
nativecontacts name chnA2 :90-105 $mem $ncParm
nativecontacts name chnB2 :231-246 $mem $ncParm
nativecontacts name chnA3 :121-135 $mem $ncParm
nativecontacts name chnB3 :262-276 $mem $ncParm
go
# Report membrane distance statistics.
avg chnA2[mindist] out chnA2.avg
avg chnB2[mindist] out chnB2.avg
avg chnA3[mindist] out chnA3.avg
avg chnB3[mindist] out chnB3.avg
go
eof
############
# Chain A ITAM2 alpha helix
# Tally number of helical residues.
n=`awk 'BEGIN {n=0};{if($1>=90 && $1<=105) n=$5+n};END {print n}' ssA.sum`
# Get shortest to membrane distance
d=`awk '{print $4}' chnA2.avg|tail -n1`
echo $name $d $n >>ssA2.h.dat
awk '{if($1>=90 && $1<=105) print $5}' ssA.sum> temp
# Paste to ssA2.dat
paste ssA2.dat temp>t
mv t ssA2.dat
##########
# Chain B ITAM2 alpha helix
# Tally number of helical residues.
n=`awk 'BEGIN {n=0};{if($1>=231 && $1<=246) n=$5+n};END {print n}' ssB.sum`
# Get shortest to membrane distance
d=`awk '{print $4}' chnB2.avg|tail -n1`
echo $name $d $n >>ssB2.h.dat
awk '{if($1>=231 && $1<=246) print $5}' ssB.sum> temp
# Paste to ssB2.dat
paste ssB2.dat temp>t
mv t ssB2.dat
############
# Chain A ITAM3 alpha helix
# Tally number of helical residues.
n=`awk 'BEGIN {n=0};{if($1>=121 && $1<=135) n=$5+n};END {print n}' ssA.sum`
# Get shortest to membrane distance
d=`awk '{print $4}' chnA3.avg|tail -n1`
echo $name $d $n >>ssA3.h.dat
awk '{if($1>=121 && $1<=135) print $5}' ssA.sum> temp
# Paste to ssA3.dat
paste ssA3.dat temp>t
mv t ssA3.dat
##########
# Chain B ITAM3 alpha helix
# Tally number of helical residues.
n=`awk 'BEGIN {n=0};{if($1>=262 && $1<=276) n=$5+n};END {print n}' ssB.sum`
# Get shortest to membrane distance
d=`awk '{print $4}' chnB3.avg|tail -n1`
echo $name $d $n >>ssB3.h.dat
awk '{if($1>=262 && $1<=276) print $5}' ssB.sum> temp
# Paste to ssB3.dat
paste ssB3.dat temp>t
mv t ssB3.dat
# Done with this model
fi
# Next model
done
done
# Delete files
rm -f temp ssA ssB ssA.sum ssB.sum
rm -f chnA2.avg chnB2.avg
rm -f chnA3.avg chnB3.avg
mv ssA2.dat tcr38.ssA2.dat
mv ssB2.dat tcr38.ssB2.dat
mv ssA3.dat tcr38.ssA3.dat
mv ssB3.dat tcr38.ssB3.dat
###########################################################
# Output from above for ITAM3 during 380-405 ns: 
# tcr38.ssA3.dat, tcr38.ssB3.dat: 
# percent helix for each model (column) and residue (row)
# ssA3.h.dat, ssB3.h.dat:
# <model><ITAM to membrane distance><# helical residues>
########
# Count how many ITAMs have any helices.
cat ssA2.h.dat ssB2.h.dat ssA3.h.dat ssB3.h.dat> temp
awk 'BEGIN {n=0};{if($3>0) n=n+1};END {print n}' temp
# 33/248 = 13% same as tcr39b
# Look at helix membrane distances 
awk '{if($3>0) print $0}' temp|sort -nk3
# Calculate mean to membrane distance for 4 resisue or greater helices.
awk 'BEGIN {n=0;d=0};{if($3>=4){n=n+1;d=d+$2}};END{print d/n}' temp
########
# If I have to transpose rows and columns:
# https://www.thelinuxrain.org/articles/transposing-rows-and-columns-3-methods
###########################################################
# Resuts for one ITAM
gnuplot<<eof
set term jpeg font "arial.ttf,18"
set out 'tcr38.helixB2.jpg'
unset key
set border 3
set xtics 0,1,44 nomirror
set ytics rotate 90 offset 1.5,0 nomirror
# Do not label xtics (model #)
set format x "" 
set format cb "%0.1f"
# Grid on top of plot
set grid x y front
set cblabel "relative α-helix frequency" 
# 
set palette model CMY rgbformulae 7,5,15
#
set view map
plot 'tcr38.ssA3.dat' matrix rowheaders with image
eof
###########################################################
# Multiplot with four panels with four ITAM results.
# S12 fig panel A
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr38.helix2.jpg'
set multiplot layout 2,2
unset key
set border 3
set xtics 0,1,61 nomirror
set ytics rotate 90 offset 1.5,0 nomirror
# Do not label xtics (model #), set colorbox format
set format x "" 
set format cb "%0.1f"
#set cbtics 0,0.2,1.0 rotate 90 
# Grid on top of plot
set grid x y front
set cblabel "relative α-helix frequency" 
# Sequential color map with zero white.
set palette model CMY rgbformulae 7,5,15
#
do for [n in "ssA2 ssB2 ssA3 ssB3"] {
  file = sprintf("tcr38.%s.dat",n)
  plot file matrix rowheaders with image
}
eof
