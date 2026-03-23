# anEdt22.tcl: Get CD3ε χ1 and χ2 ITAM angles for edt 200 ns MD.
# Output edtLck3.$mod.dat used by runE188.sh.
# Derived from anEdt20.tcl 
# vmd -dispdev none on Eire.
cd /t/tcr16/edtAn
###########################################################
# Load edt 5.4-205 ns MD trajectories for 11 edt models.
set edtModels "84 171 30 115 35 90 13 41 32 153 198"
foreach mod $edtModels {
set out [open edtLck3.$mod.dat w]
# Use stripped topology with just 2448 protein atoms that were saved in trajectory.
# made by anEdt1.tcl
set tcr [mol new /t/tcr16/edt.prot.psf]
# The NPT MD started with edt.*.eq2.
for {set i 2} {$i<=53} {incr i} {
  set mdFile /t/tcr16/edt/edt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
set lastframe [molinfo $tcr get numframes]
# Tyr chi1 angle defined by N-CA-CB-CG, chi2 is CA-CB-CG-CD1
set 188N [[atomselect $tcr "chain E and resid 188 and name N"] get index]
set 188CA [[atomselect $tcr "chain E and resid 188 and name CA"] get index]
set 188CB [[atomselect $tcr "chain E and resid 188 and name CB"] get index]
set 188CG [[atomselect $tcr "chain E and resid 188 and name CG"] get index]
set 188CD1 [[atomselect $tcr "chain E and resid 188 and name CD1"] get index]
set 199N [[atomselect $tcr "chain E and resid 199 and name N"] get index]
set 199CA [[atomselect $tcr "chain E and resid 199 and name CA"] get index]
set 199CB [[atomselect $tcr "chain E and resid 199 and name CB"] get index]
set 199CG [[atomselect $tcr "chain E and resid 199 and name CG"] get index]
set 199CD1 [[atomselect $tcr "chain E and resid 199 and name CD1"] get index]
#
for {set f 0} {$f<$lastframe} {incr f} {
# Tyr chi1 angle defined by N-CA-CB-CG, chi2 is CA-CB-CG-CD1
set 188chi1 [measure dihed [list $188N $188CA $188CB $188CG] frame $f]
set 188chi2 [measure dihed [list $188CA $188CB $188CG $188CD1] frame $f]
set 199chi1 [measure dihed [list $199N $199CA $199CB $199CG] frame $f]
set 199chi2 [measure dihed [list $199CA $199CB $199CG $199CD1] frame $f]
puts $out [format "%6d%7.1f%7.1f%7.1f%7.1f" $f $188chi1 $188chi2 $199chi1 $199chi2]
}
close $out
foreach n [atomselect list] {$n delete}
mol delete $tcr
}
###########################################################
# Combine ITAM RMSD with χ1 and χ2 angles to find a low RMSD model with 
# χ1 and χ2 close to -60 (gauche-) 
tclsh<<"eof"
# A single output with all models
set out [open edtLck4.dat w]
set models {84 171 30 115 35 90 13 41 32 153 198} 
foreach mod $models {
puts $mod
# paste two files together
exec paste edtLck1.$mod.dat edtLck3.$mod.dat > temp
set in [open temp r]
set inData [read $in]
close $in
#
foreach line [split $inData "\n"] {
  if {$line == ""} break
  set tm [lindex $line 0]
  # chain E Tyr188
  set rms188 [lindex $line 13]
  # chain E Tyr199
  set rms199 [lindex $line 19]
  set 188chi1 [lindex $line 26]
  set 188chi2 [lindex $line 27]
  set 199chi1 [lindex $line 28]
  set 199chi2 [lindex $line 29]
  # Output:
  # <mod><zero-based frame><rms188><188chi1><188chi2><rms199><199chi1><199chi2>
  puts $out [format "%4d %5d %5.2f %7.1f %7.1f %5.2f %7.1f %7.1f" \
   $mod $tm $rms188 $188chi1 $188chi2 $rms199 $199chi1 $199chi2]
  }
}
close $out
exec rm temp
eof
# ^ The above terminated early before going through all edt models.
# Nevertheless, sort -nk3 edtLck4.dat|head -n20 gave a Tyr188 model with two gauche- angles.
 35  1941  1.61   -64.4   -65.7  2.29   179.4    90.1
# sort -nk6 edtLck4.dat|head -n20 gave a Ty199 model
 84 13851  3.55   -77.5    99.3  1.47   -50.5   -75.5
# The early termination was because the anEdt20.tcl output was in error in a line that
# had a larger number of proximal contacts, it left no white space between fields on 
# a line so that data on that line was not parsed correctly. 
# Edited anEdt20.tcl and added code to check the number of fields on the line.

