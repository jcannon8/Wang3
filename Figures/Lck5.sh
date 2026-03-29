# Lck5.sh: Analysis of E188, E199 trajectories (Fig 3C,D).
# derived from Lckan.sh
cd /home/jcannon/lck
# Make no water topologies
cpptraj<<eof
parm E188.top
# Last protein atom is 4641.
parmstrip !@1-4641
parmwrite out E188.nowat.top
eof
cpptraj<<eof
parm E199.top
# Last protein atom is 4666.
parmstrip !@1-4666
parmwrite out E199.nowat.top
eof
# 
rm -f trajx
for ((i=2; i<=11;i++)); do
# trajectories are 100 ns at 10 ps/frame
echo trajin E188.eq$i.cdf >> trajx
done
# Calculate peptide to Lck LIE and RMSD of peptide CA RMSD
cpptraj<<eof
parm E188.nowat.top 
# Load trajectories like "source"
readinput trajx
# Lck residues 1-273, peptide 274-288, pTyr281
lie bind :1-273 :274-288 out E188.lie.dat
rmsd pep9 :277-285@CA out E188.rms.dat first
go
eof
#
rm -f trajx
for ((i=2; i<=11;i++)); do
# trajectories are 100 ns at 10 ps/frame
echo trajin E199.eq$i.cdf >> trajx
done
# Calculate peptide to Lck LIE and RMSD of peptide CA RMSD
cpptraj<<eof
parm E199.nowat.top 
# Load trajectories like "source"
readinput trajx
# Lck residues 1-273, peptide 274-288, pTyr281
lie bind :1-273 :274-288 out E199.lie.dat
rmsd pep9 :277-285@CA out E199.rms.dat first
go
eof
rm -f trajx
###########################################################
# Fig 3C,D
# Plot RMSD and LIE time courses and distributions.
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'lck5.lie.jpg'
# layout: rows, columns
set multiplot layout 1,2 
# Function for binning
bin(x,s) = s*int(x/s)
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
g(e,v) = 0.18 * v + 0.33 * e
set border 3
set key right textcolor variable samplen -1
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "green" lt 1 lw 2 pt 11 pi -1 ps 1.0
# LIE distribution
set xrange [-150:0]
set yrange [0:100]
#set xtic -140,20,0 out nomirror
set xtic -150,30,0 out nomirror
set ytic 0,20,100 out nomirror offset 1,0
set xlabel "{/:Bold Interaction energy (kcal/mol)}"
set ylabel "{/:Bold Relative frequency}" offset 2,0
r=0.4
# ^ r modifies the y-axis scale, but not the distribution shape.
plot "lckPep.lie.dat" u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 4 t "ζ111",\
'E188.lie.dat' u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 1 t "ε188",\
'E199.lie.dat' u (bin((g(\$2,\$3)),1.0)):(r) s f w lines ls 2 t "ε199"
#
set xrange [0:4]
set yrange [-140:0]
set ytic -140,20,0 out nomirror offset 1,0
set grid ytics
unset xlabel
set ylabel "{/:Bold Interaction energy}\n{/:Bold (kcal/mol)}" offset 1,0
set xtics ("ζ111" 1, "ε188" 2, "ε199" 3)
set style boxplot nooutliers pointtype 7
set style data boxplot
set boxwidth 0.5
set key off
# 0-100 ns data
plot "lckPep.lie.dat" u (1):(g(\$2,\$3))  w boxplot ls 4 t "",\
'E188.lie.dat' u (2):(g(\$2,\$3)) w boxplot ls 1 t "",\
'E199.lie.dat' u (3):(g(\$2,\$3)) w boxplot ls 2 t ""
# 50-100 data
#plot "lckPep.lie.dat" u (1):(\$1>5000?g(\$2,\$3):1/0) w boxplot ls 4 t "",\
#'E188.lie.dat' u (2):(\$1>5000?g(\$2,\$3):1/0) w boxplot ls 1 t "",\
#'E199.lie.dat' u (3):(\$1>5000?g(\$2,\$3):1/0) w boxplot ls 2 t ""
# 90-100 data
#plot "lckPep.lie.dat" u (1):(\$1>9000?g(\$2,\$3):1/0) w boxplot ls 4 t "",\
#'E188.lie.dat' u (2):(\$1>9000?g(\$2,\$3):1/0) w boxplot ls 1 t "",\
#'E199.lie.dat' u (3):(\$1>9000?g(\$2,\$3):1/0) w boxplot ls 2 t ""
quit
eof
###########################################################
# Use R to get the mean and standard deviation of LIE distribution.
# https://cran.r-project.org/doc/contrib/Short-refcard.pdf
read.table("lckPep.lie.dat")->a
# Function for linear interaction energy (alpha=0.18, beta=0.33)
# Aqvist et al. 1994; Carlsson et al. 2006
a[,2]*0.33+a[,3]*0.18->b
summary(b)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -137.36  -86.60  -73.15  -72.61  -58.69  -11.72
sd(b)
# 19.15845
# Get statistics of RMSD distributions
read.table("lckPep.rms.dat")->c
# 15-residue peptide statistics
summary(c[,2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   1.918   2.207   2.311   2.622   4.963
sd(c[,2])
# 0.5462873
# 9-residue statistics
summary(c[,3])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   1.115   1.255   1.273   1.398   3.474
sd(c[,3])
# 0.262934
######################################
# Compare the Lck-ITAM LIE for 0-100 ns
read.table("lckPep.lie.dat")->a
read.table("E188.lie.dat")->b
read.table("E199.lie.dat")->g
# g(e,v) = 0.18 * v + 0.33 * e
0.18*a[,3]+0.33*a[,2]->c
0.18*b[,3]+0.33*b[,2]->d
0.18*g[,3]+0.33*g[,2]->j
t.test(c,d)
# t = 33.064, df = 19993, p-value < 2.2e-16
# 90-100 ns
0.18*a[9000:10000,2]+0.33*a[9000:10000,3]->c
0.18*b[9000:10000,2]+0.33*b[9000:10000,3]->d
t.test(c,d)
# t = 33.064, df = 19993, p-value < 2.2e-16

