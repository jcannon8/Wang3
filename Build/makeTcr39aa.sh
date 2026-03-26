#!/bin/bash
# makeTcr39aa.sh: Make several tcr39.*.** models
export mod=8
export aRef=0
export bRef=0
vmd -dispdev none -e makeTcr39c.tcl
export aRef=0
export bRef=1
vmd -dispdev none -e makeTcr39c.tcl
export aRef=0
export bRef=2
vmd -dispdev none -e makeTcr39c.tcl
export aRef=1
export bRef=0
vmd -dispdev none -e makeTcr39c.tcl
export aRef=1
export bRef=1
vmd -dispdev none -e makeTcr39c.tcl
export aRef=1
export bRef=2
vmd -dispdev none -e makeTcr39c.tcl
export aRef=2
export bRef=0
vmd -dispdev none -e makeTcr39c.tcl
export aRef=2
export bRef=1
vmd -dispdev none -e makeTcr39c.tcl
export aRef=2
export bRef=2
vmd -dispdev none -e makeTcr39c.tcl
#
for m in 00 01 02 10 11 12 20 21 22; do
echo model $m
sort -nrk1 tcr39.$mod.$m.dat|tail -n2
echo
done
###########################################################
# tcr39.0.??
#model 00
#  0 10759  5787 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 01
#  0 10414  9867 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0 10477  7966 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 10
#  0  1575  9159 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 11
#  0  2640  1327 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 12
#  0  1183 10754 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 20
#  0  5687  8636 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 21
#  0 10974  5987 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0  5718  3787 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# tcr39.1.?? Those did not fuse
###########################################################
# tcr39.2.??
#model 00
#  2  6014  3119 |  0  0 | 0 0 0 0 0 0 0 2 0 Use this one
#model 10
#  1  5452  4657 |  0  1 | 0 0 0 0 0 0 0 0 0 Use this one
#model 20
#  2 10622 10088 |  0  1 | 0 0 0 1 0 0 0 0 0
#model 01
#  1  8000  1520 |  0  0 | 0 0 0 1 0 0 0 0 0 Use this one
#model 11
#  1  3290  6510 |  1  0 | 0 0 0 0 0 0 0 0 0 Use this one
#model 21
#  0 10929  6510 |  0  0 | 0 0 0 0 0 0 0 0 0 Use this one
#model 02
#  1  8000  2815 |  0  0 | 0 0 0 0 0 1 0 0 0 Use this one
#model 12 
#  1  5452  2815 |  0  0 | 0 0 0 0 0 1 0 0 0 Use this one
#model 22
#  2  6979  7466 |  0  1 | 0 1 0 0 0 0 0 0 0
###########################################################
# tcr39.3.??
#model 00
#Model tcr39.3.00 is not feasible, >1000 attempts
#model 01
#Model tcr39.3.01 is not feasible, >1000 attempts
#model 02
#  2  6008  3899 |  0  1 | 0 0 0 0 0 0 0 1 0
#model 10
#Model tcr39.3.10 is not feasible, >1000 attempts
#model 11
#  5  1132  5444 |  0  1 | 0 0 0 2 0 0 0 2 0
#model 12
#  2  1132  3899 |  0  1 | 0 0 0 0 0 0 0 1 0
#model 20
#Model tcr39.3.20 is not feasible, >1000 attempts
#model 21
#  5   210  5444 |  0  1 | 0 0 0 2 0 0 0 2 0
#model 22
  2   210  3899 |  0  1 | 0 0 0 0 0 0 0 1 0
###########################################################
# tcr39.4.??
#model 00
#  0  5968 10058 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 01
#  0  5953  3974 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0  5968  7328 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 10
#  0  1168  5357 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 11
#  0  6068  6075 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 12
#  0  3329  7663 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 20
#  0  8981 10057 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 21
#  0  8311  3970 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0  8981  7663 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# tcr39.5.??
#model 00
#  0  7762  4463 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 01
#  0  3641  6683 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 02
#  0  5963  7033 |  0  0 | 0 0 0 0 0 0 0 0 0 use
#model 10
#  0  1964  8557 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 11
#  0 10386  4658 |  0  0 | 0 0 0 0 0 0 0 0 0 use
#model 12
#  1   913  7584 |  0  0 | 0 0 0 0 0 1 0 0 0 use
#model 20
#  0 10358  5635 |  0  0 | 0 0 0 0 0 0 0 0 0 Use
#model 21
#  0  8955  2982 |  0  0 | 0 0 0 0 0 0 0 0 0
#model 22
#  0  6676  7415 |  0  0 | 0 0 0 0 0 0 0 0 0
###########################################################
# tcr39.7.??
#model 00
#  5  3014  8670 |  0  0 | 0 0 0 5 0 0 0 0 0
#model 01
#Model tcr39.7.01 is not feasible, >1000 attempts
#model 02
#  8 10814  7057 |  0  4 | 0 4 0 0 0 0 0 0 0
#model 10
#  5 10477  8670 |  0  0 | 0 0 0 5 0 0 0 0 0
#model 11
#  2  1182 10680 |  0  1 | 0 1 0 0 0 0 0 0 0
#model 12
#Model tcr39.7.12 is not feasible, >1000 attempts
#model 20
#  5  3937  8670 |  0  0 | 0 0 0 5 0 0 0 0 0
#model 21
#  2 10097 10680 |  0  1 | 0 1 0 0 0 0 0 0 0
#model 22
#  8 10097  7057 |  0  4 | 0 4 0 0 0 0 0 0 0
