# anFgt8.sh: Analyze fgt membrane association of KR CT residues (Fig 1C,D).
# Remember VMD residue and index are zero-based. 
# Amber residue and index are one-based.
cd /t/tcr16/fgtAn
# Analyze PRS and ITAM domain distance to membrane in chains F and G.
# fgt models have F and G residues in that order.
# Use VMD to convert chain resid into Amber residue using fgt.84.pdb. 
vmd -dispdev none
#set fgt [mol new /t/tcr16/fgt/fgt.75.pdb]
# ^ Either of these examples give the same residue numbers below.
set fgt [mol new /t/tcr16/fgt/fgt.75.psf]
# Look for Arg and Lys landing on membrane.
# Chain F CT resid 157-207, chain G resid 139-182
set fRK [atomselect $fgt "chain F and resid 157 to 207 and resname ARG LYS and name CA"]
$fRK get {resname resid residue}; # 12 residues
set gRK [atomselect $fgt "chain G and resid 139 to 182 and resname ARG LYS and name CA"]
$gRK get {resname resid residue}; # 7 residues
# Get the first and last membrane residues
set mem [atomselect $fgt "chain J"] 
set memI [$mem get residue]
lindex $memI 0; lindex $memI end; # Residues 161-1000
quit
###########################################################
# Models with 215 ns or 415 ns
models="185 121 191 180 137 75 14 153 93 28 109 122 42 45 "
for n in $models;do 
name=fgt.$n
rm -f trajx time.dat
echo time>time.dat
# The heating ended with fgt.*.eq1.
#for ((i=2; i<=55;i++)); do
for ((i=2; i<=94;i++)); do
# only if file exists
if [ -e /t/tcr16/fgt/fgt.$n.eq$i.rst ]; then
echo trajin /t/tcr16/fgt/fgt.$n.eq$i.rst >> trajx
# Get time from restart file, just integer
ncdump -v time /t/tcr16/fgt/fgt.$n.eq$i.rst| \
grep "time ="|awk '{printf "%7.0f\n", $3}' >>time.dat
fi
done
#
out=$name.mem3.dat
cpptraj<<eof
# cpptraj cannot handle CHARMM31 format topology files, use pdb instead
parm /t/tcr16/fgt/fgt.$n.pdb noconect
# Load trajectories like "source"
readinput trajx
# Minimum distance between Chain F RK and membrane (:217-1056).
# Name is 6JXR residue
nativecontacts name R179 :57 :161-1000&!@H= mindist
nativecontacts name R192 :70 :161-1000&!@H= mindist
nativecontacts name K193 :71 :161-1000&!@H= mindist
nativecontacts name R196 :74 :161-1000&!@H= mindist
nativecontacts name R205 :83 :161-1000&!@H= mindist
nativecontacts name R206 :84 :161-1000&!@H= mindist
# Minimum distance between Chain G RK and membrane (:217-1056).
nativecontacts name R143 :121 :161-1000&!@H= mindist
nativecontacts name R146 :124 :161-1000&!@H= mindist
nativecontacts name K150 :128 :161-1000&!@H= mindist
nativecontacts name K164 :142 :161-1000&!@H= mindist
nativecontacts name R166 :144 :161-1000&!@H= mindist
nativecontacts name R180 :158 :161-1000&!@H= mindist
nativecontacts name R181 :159 :161-1000&!@H= mindist
# This "go" is necessary before the write command.
go
# These data files cannot be combined into one output in ccptraj
# Chain F  
write R179.dat R179[mindist]
write R192.dat R192[mindist]
write K193.dat K193[mindist]
write R196.dat R196[mindist]
write R205.dat R205[mindist]
write R206.dat R206[mindist]
# Below are for Chain G
write R143.dat R143[mindist] 
write R146.dat R146[mindist]
write K150.dat K150[mindist] 
write K164.dat K164[mindist]
write R166.dat R166[mindist]
write R180.dat R180[mindist]
write R181.dat R181[mindist]
go
eof
paste time.dat R179.dat R192.dat K193.dat R196.dat R205.dat R206.dat \
R143.dat R146.dat K150.dat K164.dat R166.dat R180.dat R181.dat > $out
rm R179.dat R192.dat K193.dat R196.dat R205.dat R206.dat
rm R143.dat R146.dat K150.dat K164.dat R166.dat R180.dat R181.dat
done
rm -f trajx time.dat
###########################################################
# Plot set of 14 fgt membrane contact by residue
models="185 121 191 180 137 75 14 153 93 28 109 122 42 45 "
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
# Adjust the size so that the height uses 320 pixels per graph
set term jpeg font "arial.ttf,18" size 1280,1280
set out 'fgt.mem04.jpg'
# layout: rows, columns 
set multiplot layout 4,4 title "{/:Bold CD3ε membrane association}"
set border 3
set rmargin 0
set lmargin 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:230]
set yrange [0:60]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 7 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("fgt.%s",n) offset 0,-1
  if (n==147) {set ylabel "{/:Bold Minimum distance to membrane (\305)}" offset 2.5,0}
  else {set ylabel " "}
#  if (n>49) {set xlabel "{/:Bold Time (ns)}"}
#  else {unset xlabel}
  filename = sprintf("fgt.%s.mem3.dat",n)
  # Plot only last six RK, BRS ends at R180, label with 6JXR resid
  plot filename using (\$1/1000) : 3 with lines ls 1 t "R179",\
  '' using (\$1/1000) : 5 with lines ls 2 t "R192",\
  '' using (\$1/1000) : 7 with lines ls 3 t "K193",\
  '' using (\$1/1000) : 9 with lines ls 4 t "R196",\
  '' using (\$1/1000) : 11 with lines ls 5 t "R205",\
  '' using (\$1/1000) : 13 with lines ls 7 t "R206" 
}
quit
eof
###########################################################
# Plot Chain G membrane contact by residue
gnuplot<<eof
#set term jpeg font "arial.ttf,18" size 1280,960
# Adjust the size so that the height uses 320 pixels per graph
set term jpeg font "arial.ttf,18" size 1280,1280
set out 'fgt.mem02.jpg'
# layout: rows, columns
set multiplot layout 4,4 title "{/:Bold CD3γ membrane association}"
set border 3
set rmargin 0
set lmargin 3
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:230]
set yrange [0:60]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 7 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("fgt.%s",n) offset 0,-1
  if (n==147) {set ylabel "{/:Bold Minimum distance to membrane (\305)}" offset 2.5,0}
  else {set ylabel " "}
  filename = sprintf("fgt.%s.mem3.dat",n)
  # Plot all seven RK, label with 6JXR resid
  plot filename using (\$1/1000) : 15 with lines ls 1 t "R143",\
  '' using (\$1/1000) : 17 with lines ls 2 t "R146",\
  '' using (\$1/1000) : 19 with lines ls 3 t "K150",\
  '' using (\$1/1000) : 21 with lines ls 4 t "K164",\
  '' using (\$1/1000) : 23 with lines ls 5 t "R166",\
  '' using (\$1/1000) : 27 with lines ls 7 t "R181"
}
quit
eof
# Eliminate the R180 so that plot is not so crowded.
#    '' using (\$1/1000) : 25 with lines ls 6 t "R180" 
###########################################################
# Use R to get average distances for all models at each time
# Files edt.%s.mem3.dat have 12 RK to membrane distances in columns 3,5,7,...
# For each row (frame), average data in each columns across all files.
models=c("185","121","191","180","137","75","14","153","93","28","109","122","42","45")
fResidues=c("R179","R192","K193","R196","R205","R206")
gResidues=c("R143","R146","K150","K164","R166","R180")
# D array has dimensions model, residue, frame
#nFrame=41; # For 0-215 ns
nFrame=83; # For 0-415 ns
D<-array(,dim=c(14,12,nFrame))
for(i in seq(1,14,1)) {
  read.table(sprintf("fgt.%s.mem3.dat",models[i]),skip=1)->a
  k=1
  # residues are in every other column starting with column 3.
  for(j in seq(3,25,2)) {
    D[i,k,]<-a[,j]
    k=k+1
}
}
# Output mean and SD of residue (r) distance for each frame (c) over all models 
M<-matrix(,nrow=nFrame,ncol=12)
S<-matrix(,nrow=nFrame,ncol=12)
for(c in seq(1,12,1)) {
  for(r in seq(1,nFrame,1)) {
    # r=residue index, c=time index
    mean(D[,c,r])->M[r,c]
    sd(D[,c,r])/2->S[r,c]
  }
}
#write.table(M,file="fgt.memAve4.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
# ^ For 0-215 ns 
write.table(M,file="fgt.memAve5.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(S,file="fgt.memSD5.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
# combine mean and SD
paste fgt.memAve5.dat fgt.memSD5.dat >fgt.memAS.dat
# Plot average distances for CD3ε Figure 1C
residues="R179 R192 K193 R196 R205 R206"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgt.mem07.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3ε membrane association in CD3εγ dimer}"
set border 3
set rmargin 1
set lmargin 5
set xtics out nomirror 0,100,400
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:430]
set yrange [0:30]
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set errorbars small
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
# Time is row number * 5 ns, CD3ε residues in columns 1-6, SD in columns 13-18.
meanCol = 1
SDcol =13
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2
  # This puts one ylabel per set of six plots.
  if (meanCol==1) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-6}
  else {set ylabel " "}
  plot 'fgt.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'fgt.memAve5.dat' u (\$0*5+5):meanCol with lines ls 2 t "",\
  meanCol = meanCol+1
  SDcol = SDcol+1
}
quit
eof
###########################################################
# Plot average distances for CD3γ Figure 1D
residues="R143 R146 K150 K164 R166 R180"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'fgt.mem08.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3γ  membrane association in CD3εγ dimer}"
set border 3
set rmargin 1
set lmargin 5
set xtics out nomirror 0,100,400
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:430]
set yrange [0:30]
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set errorbars small
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
# Time is row number * 5 ns, CD3δ residues in columns 7-12, SD in coluns 19-24.
meanCol = 7
SDcol = 19 
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2 
  # This puts one ylabel per set of six plots.
  if (meanCol==7) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-6}
  else {set ylabel " "}
  plot 'fgt.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'fgt.memAve5.dat' u (\$0*5+5):meanCol with lines ls 2 t "",\
  meanCol = meanCol+1
  SDcol = SDcol+1
}
quit
eof
# Above EPS output via Adobe Illustrator: make 1000x750 pixel frames



