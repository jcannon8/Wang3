# anTcr33.sh: Analyze CHL affinity in tcr33 (Fig 8).
cd /t/tcr33
# First look at CHL1 S1 to S2 distance.
# Use VMD to get atom masks for C25-C25 distance
vmd -dispdev none
set name tcr33
set tcr [mol new $name.psf]
mol addfile $name.eq6.rst type netcdf waitfor all
[atomselect $tcr "chain S"] get {resname segid}; # Structural CHL1 are segid S1 and S2
[atomselect $tcr "segid S1 and name C25"] get index; # zero-based index 16793
[atomselect $tcr "segid S2 and name C25"] get index; # zero-based index 16867
###########################################################
# Use cpptraj to get distances
name=tcr33b
cd /t/tcr33
# Build collection of trajectory *.rst inputs
rm -f trajx
echo time>time.dat; echo eq>eq.dat
# Up to 402 ns, eq86 for tcr33, eq6-90 for tcr33b.
for ((j=6; j<=90;j++)); do
if [ -e /t/tcr33/$name.eq$j.rst ]; then
echo trajin /t/tcr33/$name.eq$j.rst >> trajx
# Get restart file time (ps) from last frame.
grep TIME /t/tcr33/$name.eq$j.out|tail -n1|awk '{printf("%7d\n",$6)}' >>time.dat
echo "eq$j" >> eq.dat
fi
done
# Run ccptraj to get distance to membrane
cpptraj<<eof
# Here the Charmm format topology is being read without problem.
parm /t/tcr33/$name.top
# Load trajectories like "source"
readinput trajx
# Use one-based C25 indices from from above VMD code
distance CHLdist @16794 @16868 out $name.c25.dat
#autoimage anchor ^1
# Output autoimaged trajectory for tcr33.eq[6-86], 81 frames.
#trajout tcr33.eq46a.cdf cdf
eof
paste time.dat eq.dat $name.c25.dat > $name.c25a.dat
###########################################################
# Plot C25-C25 CHL1 distance 
# Fig 8B
gnuplot<<eof
set term jpeg font "arial.ttf,15" size 500,375
set out 'tcr33.c25b.jpg'
set xlabel "{/:Bold Time (ns)}"
set ylabel "{/:Bold C25-C25 CHL1 distance (\305)}" offset 1,0
set border 3
set xtics out nomirror 
set ytics out nomirror
set xtics 0,100,400 
set key left textcolor variable samplen -1
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
plot 'tcr33.c25a.dat' using (\$1/1000) : 4 with linespoints ls 1 t "TCR^{MA}", \
'tcr33b.c25a.dat' using (\$1/1000) : 4 with linespoints ls 2 t "TCR^{MB}"
eof
###########################################################
cpptraj<<eof
parm /t/$name/$name.top
# Only protein plus structural CHL saved.
parmstrip !@1-16877
# Load eq44-45 trajectory
trajin /t/$name/$name.eq44.cdf
trajin /t/$name/$name.eq45.cdf
# Use one-based C25 indices from above VMD code
distance CHLdist @16794 @16868 out tcr33.c25b.dat
eof
# from less tcr33.c25a.dat
191875 eq44          39       6.2718
196875 eq45          40      10.6632
201875 eq46          41      16.6547
206875 eq47          42      18.7601
###########################################################
# Autoimage three frames
cpptraj<<eof
parm /t/$name/$name.top
trajin /t/$name/$name.eq44.rst
trajin /t/$name/$name.eq45.rst
trajin /t/$name/$name.eq46.rst
autoimage anchor ^1
trajout tcr33.eq44a.cdf cdf
eof
###########################################################
# View with vmd 
set name tcr33
set tcr [mol new $name.psf]
#mol addfile $name.eq44a.cdf type netcdf waitfor all
mol addfile $name.pdb
display projection Orthographic
display depthcue off
color Display Background white
rotate x by 90
# Find residues near CHL1 S1
set nearCHL [atomselect $tcr "same residue as within 4 of segname S1"]
lsort -unique [$nearCHL get {chain resname resid}]
{B ASP 36} {B GLY 37} {B LEU 34} {B LYS 30} {B PHE 40} {B TYR 33} 
{G ALA 113} {G ALA 121} {G GLY 117} {G PHE 120} {G THR 114} {G VAL 124} 
{J POPC 1161} {M LEU 256} {M LEU 257} {N ALA 278} {N GLY 274} {N LEU 281} 
{N LEU 285} {N SER 277} {N TYR 282} {S CHL1 201} {W TIP3 6047}
# 20 total: 6 chain B, 6 chain G, 2 chain M, 6 chain N
###########################################################
# Get table of S1 distances from tcr33.eq[1-46]
# vmd -dispdev none 
set name tcr33
set tcr [mol new $name.psf]
if {0} {
  # Try un-autoimaged files first
  set out [open tcrS1.dat w]
  for {set i 6} {$i<=46} {incr i} {
  mol addfile $name.eq$i.rst type netcdf waitfor all
  }
}
if {1} {
  # Use autoimaged files second
  set out [open tcrS1a.dat w]
  mol addfile $name.eq46a.cdf type netcdf waitfor all
}
for {set i 1} {$i<=46} {incr i} {
# NZO3: K30b NZ to S1 O3
set K30bNZ [[atomselect $tcr "name NZ and chain B and resid 30" frame $i] get {x y z}] 
set S1O3 [[atomselect $tcr "name O3 and segid S1" frame $i] get {x y z}] 
set NZO3 [vecdist [lindex $K30bNZ 0] [lindex $S1O3 0]]
# Y33bS1: Y33b sidechain to S1
set Y33b [measure center [atomselect $tcr "chain B and sidechain and resid 33" frame $i]] 
set S1 [measure center [atomselect $tcr "segid S1" frame $i]] 
set Y33bS1 [vecdist $Y33b $S1]
# F40bS1: F40b sidechain to S1
set F40b [measure center [atomselect $tcr "chain B and sidechain and resid 40" frame $i]] 
set F40bS1 [vecdist $F40b $S1]
# Y282nS1: Y282n sidechain to S1
set Y282n [measure center [atomselect $tcr "chain N and sidechain and resid 282" frame $i]] 
set Y282nS1 [vecdist $Y282n $S1]
# F120gS1: F120g sidechain to S1
set F120g [measure center [atomselect $tcr "chain G and sidechain and resid 120" frame $i]] 
set F120gS1 [vecdist $F120g $S1]
# L281nS1: L281n sidechain to S1
set L281n [measure center [atomselect $tcr "chain N and sidechain and resid 281" frame $i]] 
set L281nS1 [vecdist $L281n $S1]
# L285nS1: L285n sidechain to S1
set L285n [measure center [atomselect $tcr "chain N and sidechain and resid 285" frame $i]] 
set L285nS1 [vecdist $L285n $S1]
# POPCs1: POPC 1161 to S1
set POPC [measure center [atomselect $tcr "resname POPC and resid 1161" frame $i]] 
set POPCs1 [vecdist $POPC $S1]
# Output <eq><NZO3><Y33bS1><F40bS1><Y282nS1><F120gS1><L281nS1><L285nS1><POPCs1>
puts $out [format "%3d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" \
 [expr $i+5] $NZO3 $Y33bS1 $F40bS1 $Y282nS1 $F120gS1 $L281nS1 $L285nS1 $POPCs1]
}
close $out
quit
# tcr33.eq44 shows a larger deviation of distances to the preceeding eq43 frame.
###########################################################
# Compare original tcr33 to tcr33.eq44 with vmd
set name tcr33
set tcr [mol new $name.psf]
mol addfile $name.pdb
if {[file exists $name.eq44.rst]==0} {
  exec /t/tcr37/image1.sh $name eq44
}
# Alter display
mol addfile $name.eq44a.rst type netcdf waitfor all
display projection Orthographic
display depthcue off
color Display Background white
rotate x by 90
if {0} {
# Move frame 1 protein atoms so that CHL1 S1 fits frame 0 CHL1 S1
set protF1 [atomselect $tcr "protein or segid S1" frame 1]
set S1F0 [atomselect $tcr "segid S1 and not hydrogen" frame 0]
set S1F1 [atomselect $tcr "segid S1 and not hydrogen" frame 1]
set M [measure fit  $S1F1 $S1F0]
$protF1 move $M
measure rmsd $S1F0 $S1F1; # rmsd=1.0282
# That shows that S1 made a big change in orientation and protein did not.
}
# Move frame 1 protein and CHL1 S1 so that they fit protein near CHL in frame 0.
# This gets CA of residues near CHL1 S1.
set nearCHL0 [atomselect $tcr "name CA and same residue as within 4 of segid S1" frame 0]
# Use the indices of above to get coordinates of frame 1 atoms.
set nearCHL1 [atomselect $tcr "index [$nearCHL0 get index]" frame 1]
set protF1 [atomselect $tcr "protein or segid S1" frame 1]
set M [measure fit $nearCHL1 $nearCHL0]
$protF1 move $M
measure rmsd $nearCHL0 $nearCHL1; # rmsd=1.143
# Get 492 ns frame, eq86
mol addfile $name.eq86a.rst type netcdf waitfor all
set nearCHL2 [atomselect $tcr "index [$nearCHL0 get index]" frame 2]
# Move frame 2 protein and CHL1 S1 so that they fit protein near CHL in frame 0.
set protF2 [atomselect $tcr "protein or segid S1" frame 2]
set M [measure fit $nearCHL2 $nearCHL0]
$protF2 move $M
measure rmsd $nearCHL0 $nearCHL2; # rmsd=1.33
# Get autoimaged eq6 frame
mol addfile $name.eq46a.cdf type netcdf waitfor all
set nearCHL3 [atomselect $tcr "index [$nearCHL0 get index]" frame 3]
# Move frame 2 protein and CHL1 S1 so that they fit protein near CHL in frame 0.
set protF3 [atomselect $tcr "protein or segid S1" frame 3]
set M [measure fit $nearCHL3 $nearCHL0]
$protF3 move $M
measure rmsd $nearCHL0 $nearCHL3; # rmsd=0.468
###########################################################
# Bind the extracellular-side fgt cholesterol to tcr33.
# From chlTest3.tcl
set fgt [mol new /t/tcr16/fgt/fgt.185.psf]
# Use fgt185.eq54 (215 ns). 
mol addfile /t/tcr16/fgt/fgt.185.eq54.rst type netcdf waitfor all
# Get tcr33 model derived from 7jfd.
set tcr [mol new tcr33.psf]
mol addfile tcr33.pdb
display projection Orthographic
display depthcue off
color Display Background white
rotate x by 90
# Fit only chain G resid 108 to 138 
set CAfgtG [atomselect $fgt "name CA and chain G and resid 108 to 138"]
set ca33G [atomselect $tcr "name CA and chain G and resid 108 to 138"]
# Move fgt to fit tcr33
set allFgt [atomselect $fgt all]
set M [measure fit $CAfgtG $ca33G]
$allFgt move $M
# Now fgt CHL1 519 could replace tcr33 CHL1 S1
set S1 [atomselect $fgt "resname CHL1 and resid 519"] 
$S1 set chain S
$S1 writepdb CHL519.pdb
#
set tcr [mol new tcr33b.psf]
mol addfile tcr33b.min.rst type netcdf waitfor all
###########################################################
# Examine cholesterol orientation in tcr34
cd /t/tcr34
set tcr [mol new tcr34.psf]
mol addfile tcr34.pdb
set nearCHL [atomselect $tcr "same residue as within 4 of segname S1"]
lsort -unique [$nearCHL get {chain resname resid}]
{B GLY 37} {B ILE 41} {B LEU 34} {B LYS 30} {B PHE 40} {B TYR 33} {G ALA 113} 
{G ALA 121} {G GLY 117} {G PHE 118} {G SER 125} {G THR 114} {G VAL 124} 
{J DOPC 1164} {J DOPE 1167} {J POPS 364} {M LEU 256} {N ALA 278} {N ALA 281} 
{N ALA 285} {N GLY 274} {N SER 277} {N TYR 282} {S CHL1 201}
mol addfile tcr34.eq50a.rst type netcdf waitfor all
