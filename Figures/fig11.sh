# fig11.sh: Cluster analysis of CD3 CT ensembles (Fig 10).
# This has output for fig11A-D
# Figure 10A: CD3ε from clusterEps.sh
cd /home/jcannon/tcr2/cd3f
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11a.jpg'
# layout: rows, columns
set multiplot layout 1,2
set border 3
set title "{/:Bold CD3ε CT}" offset 0,-1 
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/home/jcannon/tcr2/cd3f/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/home/jcannon/tcr2/cd3f/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
# Modified from 
# https://stackoverflow.com/questions/31896718/generation-of-pie-chart-using-gnuplot
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 14 k-means clusters}"
plot '/home/jcannon/tcr2/cd3f/whole.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10B: CD3δ from clusterDelta3.sh
#cd /home/jcannon/tcr2/cd3d
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11b.jpg'
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3δ CT" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror 
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/home/jcannon/tcr2/cd3d/hier.dat' using 1 : 3 with linespoints ls 1 title "hierarchical", \
'/home/jcannon/tcr2/cd3d/means.dat' using 1 : 3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 12 k-means clusters}"
plot '/home/jcannon/tcr2/cd3d/means12.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10C: CD3γ from clusterGamma3.sh
#cd /home/jcannon/tcr2/cd3g
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11c.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3γ CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/home/jcannon/tcr2/cd3g/hier.dat' using 1 : 3 with linespoints ls 1 title "hierarchical", \
'/home/jcannon/tcr2/cd3g/means.dat' using 1 : 3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 10 k-means clusters}"
plot '/home/jcannon/tcr2/cd3g/means10.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10D: CD3ζ from clusterZeta3.sh 
#cd /home/jcannon/tcr2/cd3z
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11d.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3ζ CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/home/jcannon/tcr2/cd3z/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/home/jcannon/tcr2/cd3z/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 19 k-means clusters}"
plot '/home/jcannon/tcr2/cd3z/means19.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10E: CD3εδ from anEdt1.tcl
#cd /home/jcannon/tcr16/edt
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11e.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3εδ CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/t/tcr16/edt/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/t/tcr16/edt/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 14 k-means clusters}"
plot '/t/tcr16/edt/edt14.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10F: CD3εγ from anFgt1.tcl
#cd /home/jcannon/tcr16/fgt
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11f.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3εγ CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/t/tcr16/fgt/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/t/tcr16/fgt/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 19 k-means clusters}"
plot '/t/tcr16/fgt/fgt19.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10G: CD3εδR from anEdr1.tcl
#cd /home/jcannon/tcr16/edr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11g.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3εδ^R CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/t/tcr16/edr/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/t/tcr16/edr/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 16 k-means clusters}"
plot '/t/tcr16/edr/edr16.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10H: CD3εγR from anFgr1.tcl
#cd /home/jcannon/tcr16/fgr
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11h.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3εγ^R CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/t/tcr16/fgr/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/t/tcr16/fgr/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 20 k-means clusters}"
plot '/t/tcr16/fgr/fgr20.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
###########################################################
# Figure 10I: CD3εδN from anEdr1.tcl
#cd /home/jcannon/tcr16/edn
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'fig11h.jpg'  
# layout: rows, columns
set multiplot layout 1,2
set border 3 
set title "{/:Bold CD3εδ^N CT}" offset 0,-1   
set xlabel "{/:Bold Cluster number}" offset 0,0.5
set ylabel "{/:Bold DBI}" offset 1,0
set format y "%.2f"
set xtics out nomirror
set ytics out nomirror
set grid y
set key right top textcolor variable samplen -1
set mxtics 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
plot '/t/tcr16/edn/hier.dat' using 1:3 with linespoints ls 1 title "hierarchical", \
'/t/tcr16/edn/means.dat' using 1:3 with linespoints ls 2  title "means"
##########################
# Make pie chart of cluster distribution.
unset xlabel
unset ylabel
set lmargin at screen 0.50
# ^ This gives the pie chart a smaller share of panel width.
# function for angle
ang(x)=x*360
# dimensions for pie-chart
cenX=0
cenY=0.2
# ^ This leaves room below pie chart to descibe final ensemble.
radius=0.8
set style fill solid 1     # filled pie-chart
unset key                  # no automatic labels
unset tics                 # remove tics
unset border
set size ratio -1          # equal scale length
set xrange [-1:1] 
set yrange [-1:1]
pos = 0            # init angle
color = 0          # init color
set title "{/:Bold 16 k-means clusters}"
plot '/t/tcr16/edn/edn16.dat' u (cenX):(cenY):(radius):(pos):(pos=pos+ang(\$3)):(color=color+1)\
 w circle lc var
quit
eof
