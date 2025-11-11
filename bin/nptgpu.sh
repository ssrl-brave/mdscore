#!/bin/bash

TLEAPDIR=$1
HEATDIR=$2
TEMP=$3

cat > npt1.in << EOF
1ns of equilibration
 &cntrl
   imin=0, irest=1, ntx=5,
   ntt=3, tempi=$TEMP, temp0=$TEMP, gamma_ln=1.0,
   ntp=1, pres0=1.0, taup=1.0, 
   ntb=2, ntc=2, ntf=2, cut=8.0,
   ntr=1, restraint_wt=5.0, restraintmask='!(@H=)&!(:WAT|:ETH|@Na+=)',
   nstlim=50000, dt=0.001, t=0.0, nscm=5000,
   ntpr=2500, ntwx=2500, ntwe=0, ntwr=25000,
   ioutfm=1,
 /
EOF

pmemd.cuda -O -i npt1.in -o npt1.out -p $TLEAPDIR/complex_solvated.prmtop -c $HEATDIR/heat.rst -r npt.rst -inf npt.mdinfo -ref $HEATDIR/heat.rst

