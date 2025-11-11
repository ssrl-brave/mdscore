#!/bin/bash

TLEAPDIR=$1

cat > min1.in << EOF
 polyA-polyT 10-mer: initial minimization solvent + ions
 &cntrl
  imin   = 1,
  maxcyc = 4000,
  ncyc   = 3000,
  ntb    = 1,
  ntr    = 1,
  cut    = 10.0
  restraintmask= '!(@H=) &!(:WAT| :ETH| @Na+=)',
  restraint_wt=20,
/
EOF

pmemd.cuda -O -i min1.in -o min1.out -p ${TLEAPDIR}/complex_solvated.prmtop -c ${TLEAPDIR}/complex_solvated.inpcrd -r min1.rst -inf min1.mdinfo -ref ${TLEAPDIR}/complex_solvated.inpcrd

cat > min2.in << EOF
polyA-polyT 10-mer: initial minimization solvent + ions
 &cntrl
  imin   = 1,
  maxcyc = 2000,
  ncyc   = 1000,
  ntb    = 1,
  ntr    = 0,
  cut    = 10.0
 /
EOF

pmemd.cuda -O -i min2.in -o min2.out -p ${TLEAPDIR}/complex_solvated.prmtop -c min1.rst -r min2.rst -inf min2.mdinfo -ref min1.rst

#cp complex_solvated.prmtop min2.rst system_hmass.top trial_1_boltz
#cp complex_solvated.prmtop min2.rst system_hmass.top trial_2_boltz
#cp complex_solvated.prmtop min2.rst system_hmass.top trial_3_boltz
