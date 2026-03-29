# anEdt8.sh: Analyze edt membrane association of KR CT residues (Fig 1A,B).
# derived from tcr15An3.sh.sh
cd /t/tcr16/edtAn
# Analyze PRS and ITAM domain distance to membrane in chains E and D.
# edt models have D and E residues in that order.
# Use VMD to convert chain resid into Amber residue using edt.84.pdb. 
vmd -dispdev none
#set edt [mol new /t/tcr16/edt/edt.84.pdb]
# ^ Either of these examples give the same residue numbers below.
set edt [mol new /t/tcr16/edt/edt.84.psf]
# Look for Arg and Lys landing on membrane.
# Chain E CT resid 156-207, chain D resid 130-171
set eRK [atomselect $edt "chain E and resid 156 to 207 and resname ARG LYS and name CA"]
$eRK get {resname resid residue}; # 13 residues
set dRK [atomselect $edt "chain D and resid 130 to 171 and resname ARG LYS and name CA"]
$dRK get {resname resid residue}; # 6 residues
# Get the first and last membrane residues
set mem [atomselect $edt "chain J"] 
set memI [$mem get residue]
lindex $memI 0; lindex $memI end; # Residues 158-997 in Amber 1-based numbers
# These membrane residue numbers are the same for all edt models.
quit
###########################################################
# Reference edt models with 205 or 405 ns
models="84 171 30 115 35 90 13 41 32 153 198"
for n in $models;do 
name=edt.$n
rm -f trajx time.dat
echo time>time.dat
# The heating ended with edt.*.eq1.
#for ((i=2; i<=53;i++)); do
for ((i=2; i<=92;i++)); do
# only if file exists
if [ -e /t/tcr16/edt/edt.$n.eq$i.rst ]; then
echo trajin /t/tcr16/edt/edt.$n.eq$i.rst >> trajx
# Get time from restart file, just integer
ncdump -v time /t/tcr16/edt/edt.$n.eq$i.rst| \
grep "time ="|awk '{printf "%7.0f\n", $3}' >>time.dat
fi
done
#
# Run ccptraj to get distance to membrane
out=$name.mem3.dat
cpptraj<<eof
# cpptraj cannot handle CHARMM31 format topology files, use pdb instead
parm /t/tcr16/edt/edt.$n.pdb noconect
# Load trajectories like "source"
readinput trajx
# Minimum distance between Chain E RK and membrane (:158-997).
# name is 6JXK residue name
nativecontacts name R179 :129 :158-997&!@H= mindist
nativecontacts name R192 :142 :158-997&!@H= mindist
nativecontacts name K193 :143 :158-997&!@H= mindist
nativecontacts name R196 :146 :158-997&!@H= mindist
nativecontacts name R205 :155 :158-997&!@H= mindist
nativecontacts name R206 :156 :158-997&!@H= mindist
# Minimum distance between Chain D RK and membrane (:158-997).
nativecontacts name R132 :35 :158-997&!@H= mindist
nativecontacts name R144 :47 :158-997&!@H= mindist
nativecontacts name R153 :56 :158-997&!@H= mindist
nativecontacts name K155 :58 :158-997&!@H= mindist
nativecontacts name R169 :72 :158-997&!@H= mindist
nativecontacts name K171 :74 :158-997&!@H= mindist
# This "go" is necessary before the write command.
go
# These data files cannot be combined into one output in ccptraj
# Chain E  
write R179.dat R179[mindist]
write R192.dat R192[mindist]
write K193.dat K193[mindist]
write R196.dat R196[mindist]
write R205.dat R205[mindist]
write R206.dat R206[mindist]
# Below are for Chain D
write R132.dat R132[mindist] 
write R144.dat R144[mindist]
write R153.dat R153[mindist] 
write K155.dat K155[mindist]
write R169.dat R169[mindist]
write K171.dat K171[mindist]
go
eof
paste time.dat R179.dat R192.dat K193.dat R196.dat R205.dat R206.dat \
R132.dat R144.dat R153.dat K155.dat R169.dat K171.dat > $out
rm time.dat R179.dat R192.dat K193.dat R196.dat R205.dat R206.dat
rm R132.dat R144.dat R153.dat K155.dat R169.dat K171.dat
done
rm -f trajx time.dat
###########################################################
# Plot Chain E membrane contact by model
models="84 171 30 115 35 90 13 41 32 153 198"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.mem04.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold CD3ε membrane association}"
set border 3
set rmargin 0
set lmargin 5
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:220]
set yrange [0:30]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 7 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  if (n==147) {set ylabel "{/:Bold Minimum distance to membrane (\305)}" offset 2.5,0}
  else {set ylabel " "}
#  if (n>49) {set xlabel "{/:Bold Time (ns)}"}
#  else {unset xlabel}
  filename = sprintf("edt.%s.mem3.dat",n)
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
# Plot Chain D membrane contact by model
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.mem02.jpg'
# layout: rows, columns
set multiplot layout 3,4 title "{/:Bold CD3δ membrane association}"
set border 3
set rmargin 0
set lmargin 5
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:220]
set yrange [0:30]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 7 lc rgb "black" lt 1 lw 3 pt 7 pi -1 ps 1.0
do for [n in "$models"] {
  set title sprintf("edt.%s",n) offset 0,-1
  if (n==147) {set ylabel "{/:Bold Minimum distance to membrane (\305)}" offset 2.5,0}
  else {set ylabel " "}
  filename = sprintf("edt.%s.mem3.dat",n)
  # Plot all six RK, label with 6JXR resid
  plot filename using (\$1/1000) : 15 with lines ls 1 t "R132",\
  '' using (\$1/1000) : 17 with lines ls 2 t "R144",\
  '' using (\$1/1000) : 19 with lines ls 3 t "R153",\
  '' using (\$1/1000) : 21 with lines ls 4 t "R155",\
  '' using (\$1/1000) : 23 with lines ls 5 t "R169",\
  '' using (\$1/1000) : 25 with lines ls 7 t "R169", 
}
quit
eof
###########################################################
# Plot Chain E membrane contact by residue
models="84 171 30 115 35 90 13 41 32 153 198"
residues="R179 R192 K193 R196 R205 R206"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.mem05.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3e membrane association}"
set border 3
set rmargin 0
set lmargin 5
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:220]
set yrange [0:30]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}"
#set ylabel "{/:Bold Minimum distance to membrane (\305)}"
dataCol = 3
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-1
  # no line style, ls, use defaults
  plot for [m in "$models"] sprintf("edt.%s.mem3.dat",m) u (\$1/1000):dataCol with lines lw 3 t ""
  dataCol=dataCol+2 
}
quit
eof
###########################################################
# Plot Chain D membrane contact by residue
models="84 171 30 115 35 90 13 41 32 153 198"
residues="R132 R144 R153 R155 R169 K171"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1280,960
set out 'edt.mem06.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3d membrane association}"
set border 3
set rmargin 0
set lmargin 5
set xtics out nomirror 0,50,220
set ytics out nomirror offset 1,0
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:220]
set yrange [0:30]
set style line 1 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "green" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "purple" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 5 lc rgb "orange" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 6 lc rgb "brown" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}"
dataCol = 15
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-1
  # no line style, ls, use defaults
  plot for [m in "$models"] sprintf("edt.%s.mem3.dat",m) u (\$1/1000):dataCol with lines lw 3 t ""
  dataCol=dataCol+2 
}
quit
eof
###########################################################
# Use R to get average distances for all models at each time
# Files edt.%s.mem3.dat have 12 RK to membrane distances in columns 3,5,7,...
# For each row (frame), average data in each columns across all files.
models=c("84","171","30","115","35","90","13","41","32","153","198")
eResidues=c("R179","R192","K193","R196","R205","R206")
dResidues=c("R132","R144","R153","R155","R169","R169")
# D array has dimensions model, residue, frame
#nFrame=41; # For 0-205 ns
nFrame=81; # For 0-405 ns
#D<-array(,dim=c(11,12,41))
D<-array(,dim=c(11,12,nFrame))
for(i in seq(1,11,1)) {
  read.table(sprintf("edt.%s.mem3.dat",models[i]),skip=1)->a
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
#write.table(M,file="edt.memAve4.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(M,file="edt.memAve5.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
# ^edt.memAve4.dat has 205 ns, edt.memAve5.dat has 405 ns.
write.table(S,file="edt.memSD5.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
q("no")
###########################################################
# combine mean and SD
paste edt.memAve5.dat edt.memSD5.dat >edt.memAS.dat
# Plot average distances for CD3ε 
# Figure 1A
residues="R179 R192 K193 R196 R205 R206"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.mem07.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3ε membrane association in CD3εδ dimer}"
set border 3
set rmargin 1
set lmargin 5
set xtics out nomirror 0,100,400
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:420]
set yrange [0:30]
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set errorbars small
# Time is row number * 5 ns, CD3ε residues in columns 1-6, SD in columns 13-18.
meanCol = 1
SDcol =13
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2
  # This puts one ylabel per set of six plots.
  if (meanCol==1) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-6}
  else {set ylabel " "}
  plot 'edt.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'edt.memAve5.dat' u (\$0*5+5):meanCol with lines ls 2 t ""
  meanCol = meanCol+1
  SDcol = SDcol+1
}
quit
eof
###########################################################
# Plot average distances for CD3δ 
# Figure 1B
residues="R132 R144 R153 R155 R169 K171"
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.mem08.jpg'
# layout: rows, columns
set multiplot layout 2,3 title "{/:Bold CD3δ membrane association in CD3εδ dimer}"
set border 3
set rmargin 1
set lmargin 5
set xtics out nomirror 0,100,400
set ytics out nomirror offset 1,0
set grid ytics
set key left top horizontal textcolor variable samplen -1 font "arial.ttf,12"
set xrange [0:420]
set yrange [0:30]
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
set xlabel "{/:Bold Time (ns)}" offset 0,0.5
set errorbars small
# Time is row number * 5 ns, CD3δ residues in columns 7-12, SD in coluns 19-24.
meanCol = 7
SDcol = 19 
do for [n in "$residues"] {
  set title sprintf("%s",n) offset 0,-2 
  # This puts one ylabel per set of six plots.
  if (meanCol==7) {set ylabel "{/:Bold Average distance to membrane (\305)}" offset 3,-6}
  else {set ylabel " "}
  plot 'edt.memAS.dat' u (\$0*5+5):meanCol:SDcol with yerrorbars ls 1 t "",\
  'edt.memAve5.dat' u (\$0*5+5):meanCol with lines ls 2 t ""
  meanCol=meanCol+1
  SDcol = SDcol+1
}
quit
eof
# Above EPS output via Adobe Illustrator: make 1000x750 pixel frames
