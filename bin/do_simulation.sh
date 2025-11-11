#!/bin/bash

######## INPUTS ########
# input Ligand and receptor
LIG=$(realpath $1)
PROT=$(realpath $2)
# resired output  rootdir
trialdir=$3
# temp of the simulation Kelvin
TEMP_K=298.5
# duration of the simulation nanoseconds
DUR_NS=100
# subdirs will start with this prefix
pref=mdout
###### END  INPUTS #######

# make a trial dir and move there
mkdir -p $trialdir
trialdir=$(realpath $trialdir)
cd $trialdir

# copy ligand and receptor for bookkeeping
cp $LIG input_ligand.pdb
cp $PROT input_receptor.pdb

# use this pattern to run a command sequence where we make a folder, and run a command within the folder
# the first arg $1 is the dirname
# the rest of the args $@ should be the command
run_step () {
    local outdir=$1
    shift # Remove the first argument (directory name) so $@ now holds only the command
    echo "--- Running step in: ${outdir} ---"
    mkdir -p ${outdir}
    outdir=$(realpath $outdir)
    cd $outdir
    # after the shift, the command sequence remains
    "$@"
    cd -
    # echo the output dir so we can grab it for later use
    echo $outdir
}

# Antechamber
antdir=$(run_step $pref.antechamber antechamber -i $LIG -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2 -at gaff2 -nc 0 | tail -n 1)

# parmchk2
parmchkdir=$(run_step $pref.parmchk2 parmchk2 -i $antdir/ligand.mol2 -f mol2 -o ligand.frcmod | tail  -n 1)

# tleap
tleapdir=$(run_step $pref.tleap tleap.in.sh "$antdir/ligand.mol2" "${parmchkdir}/ligand.frcmod" "$PROT" | tail -n 1)

#parmed
parmdir=$(run_step $pref.parmed parmed.sh "$tleapdir/complex_solvated.prmtop" | tail -n 1 )

# initial minimization
mindir=$(run_step $pref.min min.sh $tleapdir | tail -n 1)

# heat system to desired temp
heatdir=$(run_step $pref.heat heatgpu.sh $tleapdir $mindir $TEMP_K | tail -n 1)

# equilibration
nptdir=$(run_step $pref.npt nptgpu.sh $tleapdir $heatdir $TEMP_K | tail -n 1)

# dynamics simulation
mddir=$(run_step $pref.md md.sh $tleapdir $nptdir $parmdir $TEMP_K $DUR_NS | tail -n 1)

echo Final results written to $mddir
echo Done.

