#!/bin/bash
# makeTcr38aa.sh: Make several tcr38.*.** models
export mod=11
export aRef=0
export bRef=0
vmd -dispdev none -e makeTcr38c.tcl
export aRef=0
export bRef=1
vmd -dispdev none -e makeTcr38c.tcl
export aRef=0
export bRef=2
vmd -dispdev none -e makeTcr38c.tcl
export aRef=1
export bRef=0
vmd -dispdev none -e makeTcr38c.tcl
export aRef=1
export bRef=1
vmd -dispdev none -e makeTcr38c.tcl
export aRef=1
export bRef=2
vmd -dispdev none -e makeTcr38c.tcl
export aRef=2
export bRef=0
vmd -dispdev none -e makeTcr38c.tcl
export aRef=2
export bRef=1
vmd -dispdev none -e makeTcr38c.tcl
export aRef=2
export bRef=2
vmd -dispdev none -e makeTcr38c.tcl
#
for m in 00 01 02 10 11 12 20 21 22; do
echo model $m
sort -nrk1 tcr38.$mod.$m.dat|tail -n2
echo
done
##################
# For tcr38.2.??
#model 00
#  0  4082  8727 |  0  0 | 0 0 0 0 0 0 0 0 0 for old/tcr38.2.00
#  1  2148   307 |  0  0 | 0 0 0 0 0 1 0 0 0 for tcr38.2.00
#model 01
#  0  8208  8163 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0  6349  4243 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 10
#  0  5719  7071 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 11
#  0  9108  7869 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 12
#  0  6269  6934 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 20
#  0  8321  5222 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 21
#  0  5588  9095 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0  9078 10989 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# for tcr38.3.??
#model 00
#  0  5739 10090 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 01
#  0  8044  2068 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0  2262  2707 |  0  0 | 0 0 0 0 0 0 0 0 0 Early termination (eq3), periodic box crossing?
#  0 10087  9013 |  0  0 | 0 0 0 0 0 0 0 0 0 
#model 10
#  1  2265  6690 |  0  1 | 0 0 0 0 0 0 0 0 0
#model 11
#  0  4771  2068 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 12
#  0  4780  2707 |  0  0 | 0 0 0 0 0 0 0 0 0 Pro loop penetrated
#  1  4323  7673 |  0  1 | 0 0 0 0 0 0 0 0 0 Early terminations (eq4-11), periodic box crossing?
#  1  5718  7673 |  0  1 | 0 0 0 0 0 0 0 0 0
#model 20
#  0  4468  4184 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 21
#  0 10146 10719 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0 10806  8583 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# For tcr38.5.??
#model 00
#  0  4566  4070 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 01
#  0  6349  6381 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0  1196  8131 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 10
#  0  8300 10165 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 11
#  0  4261  8480 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 12
#  0  5396  1336 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 20
#  0  4879 10165 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 21
#  0  6940  2479 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0  3895  4734 |  0  0 | 0 0 0 0 0 0 0 0 0
# For tcr38.7.??
###########################################################
# For tcr38.7.??
#model 00
#Model tcr38.7.00 is not feasible
#model 01
#  4  3103 10826 |  0  0 | 0 0 0 3 0 1 0 0 0
# model 02
#  5  7737  6712 |  0  0 | 0 0 0 3 0 0 0 2 0
#model 10
#Model tcr38.7.10 is not feasible
#model 11
#Model tcr38.7.11 is not feasible
#model 12
#  5  4217  6712 |  0  0 | 0 0 0 3 0 0 0 2 0
#model 20
#Model tcr38.7.20 is not feasible
#model 21
#Model tcr38.7.21 is not feasible
#model 22
#  5  9061  6712 |  0  0 | 0 0 0 3 0 0 0 2 0
###########################################################
# For tcr38.8.??
#model 00
#  0  1925  7266 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 10
#  2  9644  7864 |  1  0 | 0 0 0 0 0 0 0 1 0 Use this one
#model 20
#  0  6176  7266 |  0  0 | 0 0 0 0 0 0 0 0 0 Use this one
#model 01
#  2  8611   150 |  0  2 | 0 0 0 0 0 0 0 0 0 Use this one
#model 11
#  2  9224  3011 |  0  0 | 0 0 0 0 0 0 0 2 0 Use this one
#model 21
#  2  5750  7645 |  0  0 | 0 0 0 0 0 0 0 2 0
#model 02
#  1  7799  7496 |  0  0 | 0 1 0 0 0 0 0 0 0 Use this one
#model 12
#  1  7120  9012 |  0  0 | 0 0 0 0 0 0 0 1 0 Use this one
#model 22
#  0 10146  8696 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# For tcr38.10.??
# model 00
#  1  2148  8105 |  0  0 | 0 0 0 0 0 1 0 0 0
#model 01
#  0  8217  9523 |  0  0 | 0 0 0 0 0 0 0 0 0 
#model 02
#  0  5172  8372 |  0  0 | 0 0 0 0 0 0 0 0 0 **
#model 10
#Model tcr38.10.10 is not feasible
#model 11
#  2  1498  3305 |  0  0 | 0 0 0 0 0 2 0 0 0 Use this one
#model 12
#  0  2155  8372 |  0  0 | 0 0 0 0 0 0 0 0 0 **
#model 20
#  1  3416  8105 |  0  0 | 0 0 0 0 0 1 0 0 0
#model 21
#  0  5964  9523 |  0  0 | 0 0 0 0 0 0 0 0 0 
#model 22
#  0  3875  8372 |  0  0 | 0 0 0 0 0 0 0 0 0 **
###########################################################
# For tcr38.11.??
model 00
  0  6030  6987 |  0  0 | 0 0 0 0 0 0 0 0 0
model 01
  0  7097 10680 |  0  0 | 0 0 0 0 0 0 0 0 0
model 02
  0 10733  6110 |  0  0 | 0 0 0 0 0 0 0 0 0
model 10
  0  6810 10481 |  0  0 | 0 0 0 0 0 0 0 0 0
model 11
  0  6439  5329 |  0  0 | 0 0 0 0 0 0 0 0 0
model 12
  0  8649  7057 |  0  0 | 0 0 0 0 0 0 0 0 0
model 20
  0  7907  6987 |  0  0 | 0 0 0 0 0 0 0 0 0
model 21
  0  8146  1664 |  0  0 | 0 0 0 0 0 0 0 0 0
model 22
  0  3212  6680 |  0  0 | 0 0 0 0 0 0 0 0 0

