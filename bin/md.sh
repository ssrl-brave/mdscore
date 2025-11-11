#!/bin/bash
TLEAPDIR=$1
NPT_CPU_DIR=$2
PARMDIR=$3
TEMP=298.15
DUR_NS=100
DT_PS=0.004
NSTEPS=$(echo "scale=0; (${DUR_NS} * 1000) / ${DT_PS}" | bc)

cat > npt1.in  << EOF
1ns of equilibration
 &cntrl
   imin=0, irest=1, ntx=5,
   ntt=3, tempi=$TEMP, temp0=$TEMP, gamma_ln=1.0,
   ntp=1, pres0=1.0, taup=1.0, 
   ntb=2, ntc=2, ntf=2, cut=8.0,
   
   ntr=1, restraint_wt=5.0, restraintmask='!(@H=)&!(:WAT|:ETH|@Na+=)',
   nstlim=3000000, dt=0.001, t=0.0, nscm=5000,
   ntpr=25000, ntwx=25000, ntwe=0, ntwr=500000,
   ioutfm=1,
 /
EOF

#pmemd.cuda -O -i npt1.in -o npt1.out -p $TLEAPDIR/complex_solvated.prmtop -c $NPT_CPU_DIR/npt_cpu.rst -r npt1.rst -x npt1.nc -inf npt1.mdinfo -ref $NPT_CPU_DIR/npt_cpu.rst

cat > npt2.in << EOF
1ns of equilibration
 &cntrl
   imin=0, irest=1, ntx=5,
   ntt=3, tempi=$TEMP, temp0=$TEMP, gamma_ln=1.0,
   ntp=1, pres0=1.0, taup=1.0, 
   ntb=2, ntc=2, ntf=2, cut=8.0,

   ntr=0,
   nstlim=5000000, dt=0.001, t=0.0, nscm=5000,
   ntpr=25000, ntwx=25000, ntwe=0, ntwr=500000,
   ioutfm=1,
 /
EOF

#pmemd.cuda -O -i npt2.in -o npt2.out -p $TLEAPDIR/complex_solvated.prmtop -c npt1.rst -r npt2.rst -x npt2.nc -inf npt2.mdinfo -ref npt1.rst

cat > md.in << EOF
500ns of production - NTP - Langevin thermostat - gamma 0.01 - Monte Carlo barostat - BOOST
 &cntrl
   imin=0, irest=1, ntx=5,ig=-1,
   ntt=3, tempi=$TEMP, temp0=$TEMP, gamma_ln=1.0,
   ntp=1, pres0=1.0, taup=1.0,
   ntb=2, ntc=2, ntf=2, cut=8.0,
   ntr=0, tol=0.0000001,
   nstlim=$NSTEPS, dt=${DT_PS}, nscm=5000,barostat=2,
   ntpr=25000, ntwx=25000, ntwe=0, ntwr=-5000000, iwrap=1,
   ioutfm=1,
 /
&ewald
  netfrc=0,
  skin_permit=0.75,
&end
EOF

pmemd.cuda -O -i md.in -o md.out -p $PARMDIR/system_hmass.top -c npt2.rst -r md.rst -x md.nc -inf md.mdinfo -ref npt2.rst

