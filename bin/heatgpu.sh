#!/bin/bash

TLEAPDIR=$1
MINDIR=$2
TEMP=$3

cat > mdt.in << EOF
1ns of thermalization
 &cntrl
   imin=0, irest=0, ntx=1,
   ntt=3, tempi=0.0, temp0=${TEMP}, gamma_ln=1.0,
   ntb=1, ntc=2, ntf=2, cut=8.0,
   ntr=1, restraint_wt=5.0, restraintmask='!(@H=)&!(:WAT|:ETH|@Na+=)',
   nstlim=1000000, dt=0.001, t=0.0, nscm=5000,
   ntpr=500000, ntwx=500000, ntwe=0, ntwr=50000,
   ioutfm=1,
 /
EOF

pmemd.cuda -O -i mdt.in -o heat.out -p ${TLEAPDIR}/complex_solvated.prmtop -c $MINDIR/min2.rst -r heat.rst -x heat.nc -inf heat.mdinfo -ref $MINDIR/min2.rst
