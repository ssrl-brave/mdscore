# MD simulation workflow

### Installing Amber (for `pmemd.cuda`)  and AmbterTools

This was tested on debian 12. First, one needs a valid cuda installation. You will need `nvcc` in path, e.g.

```bash
nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2024 NVIDIA Corporation
Built on Tue_Oct_29_23:50:19_PDT_2024
Cuda compilation tools, release 12.6, V12.6.85
Build cuda_12.6.r12.6/compiler.35059454_0
```

Then, create a conda environment with amberTools and the necessary components to later build pmemd:

```bash
# Download conda with mamba backend
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash ./Miniforge3-Linux-x86_64.sh -b -u -p $PWD/amberforge
source amberforge/etc/profile.d/conda.sh

# make the environment
mamba create -n amber -c conda-forge flex bison plumed patch bzip2 libzip ambertools openmpi gfortran=12.2.0 cxx_compiler conda -y

# tweak the numpy installation
conda activate amber
conda install numpy=1.24.4 -y
```

Note the specific versioning of cxx_compiler, gfortran (and dep gcc) should all be compatible with the cuda version. In my case, I knew that nvcc version 12.6 is incompatible with gcc versions >= 13. Also, the downgrade of numpy in my case was needed for `parmed` to work properly, but the [devs are aware](https://github.com/ParmEd/ParmEd/issues/1386). 

Next, download the Amber24 (pmemd) package by going [here](https://ambermd.org/GetAmber.php) and filling out the form under "Getting Amber24 for non-commercial use". To install: 

```bash
# untar the archive
tar -xjvf pmemd24.tar.bz2
cd pmemd24_src
./update_pmemd  --update
cd build
# open the file run_cmake in text editor
```

At this point, for Debian 12, it was necessary to tell `cmake` about `kmmd`, which should be in the amber CONDA env. For that `-DKMMD_DIR=$CONDA_PREFIX` was added to the `cmake` command in `run_cmake`, e.g.

```bash
#  Assume this is Linux:

  cmake $AMBER_PREFIX/pmemd24_src \
    -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/pmemd24 \
    -DCOMPILER=GNU  \
    -DKMMD_DIR=$CONDA_PREFIX \
    -DMPI=TRUE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
    -DDOWNLOAD_MINICONDA=FALSE -DBUILD_PYTHON=FALSE \
    -DBUILD_PERL=FALSE -DBUILD_GUI=FALSE \
    -DPMEMD_ONLY=TRUE -DCHECK_UPDATES=FALSE \
    2>&1 | tee  cmake.log

```

With that change, configuration and building should commence properly

```
run_cmake
make -j 4 install
```

Good luck! If successful, this will create an install folder `pmemd24` alongside `pmemd24_src` and there will be a startup file ```pmemd24/amber.sh ``` that is meant to be sourced to set the amber environment.


### Running the MD simulation pipeline

To run the pipeline, simply source your new amber build (assuming you havent already) and prepend the `bin` folder odf this repository to your `PATH` variable:

```
# load the amber tools and conda env:
source /path/to/miniforge/etc/profile.d/conda.sh
conda activate amber
# load the pmemd24 env:
source /path/to/pmemd24/amber.sh
# load the md pipeline env from mdscore:
export PATH=/path/to/mdscore/bin:$PATH
```

The main pipeline script is `do_simuilation.sh` and it simply takes as input a `ligand.pdb` and `receptor.pdb` file. Prepared are two such files for testing purposes:

```
wget https://smb.slac.stanford.edu/~dermen/receptor_boltz.pdb
wget https://smb.slac.stanford.edu/~dermen/boltz_lig_mod.pdb
```

(TODO:add required preprocessing steps for preparing the ligand / receptor PDB files)

Now, the simulation pipeline can be run using:

```
do_simulation.sh ligand.pdb recetor.pdb sysAlpha.trial1
```

The resulting `sysAlpha.trial1` folder contains the md simulation trajectory `sysAlpha.trial1/mdout.md/md.nc`, as well as all of the intermediate input/output files created along the way. Note each major pipeline step produces its own unique output folder: 

```bash
sysAlpha.trial1
├── ligand.pdb
├── mdout.antechamber
│   ├── ANTECHAMBER_AC.AC
│   ├── ANTECHAMBER_AC.AC0
│   ├── ANTECHAMBER_AM1BCC.AC
│   ├── ANTECHAMBER_AM1BCC_PRE.AC
│   ├── ANTECHAMBER_BOND_TYPE.AC
│   ├── ANTECHAMBER_BOND_TYPE.AC0
│   ├── ATOMTYPE.INF
│   ├── ligand.mol2
│   ├── sqm.in
│   ├── sqm.out
│   └── sqm.pdb
├── mdout.heat
│   ├── heat.mdinfo
│   ├── heat.nc
│   ├── heat.out
│   ├── heat.rst
│   └── mdt.in
├── mdout.md
│   ├── md.in
│   ├── md.mdinfo
│   ├── md.nc
│   ├── md.out
│   ├── md.rst
│   ├── md.rst_125000
│   ├── npt1.in
│   ├── npt1.mdinfo
│   ├── npt1.nc
│   ├── npt1.out
│   ├── npt1.rst
│   ├── npt2.in
│   ├── npt2.mdinfo
│   ├── npt2.nc
│   ├── npt2.out
│   └── npt2.rst
├── mdout.min
│   ├── min1.in
│   ├── min1.mdinfo
│   ├── min1.out
│   ├── min1.rst
│   ├── min2.in
│   ├── min2.mdinfo
│   ├── min2.out
│   └── min2.rst
├── mdout.npt
│   ├── mdcrd
│   ├── npt1.in
│   ├── npt1.out
│   ├── npt.mdinfo
│   └── npt.rst
├── mdout.parmchk2
│   └── ligand.frcmod
├── mdout.parmed
│   ├── parmed.in
│   ├── parmed.log
│   └── system_hmass.top
├── mdout.tleap
│   ├── complex_solvated.inpcrd
│   ├── complex_solvated.pdb
│   ├── complex_solvated.prmtop
│   ├── leap.log
│   ├── ligand.lib
│   ├── ligand.prmtop
│   ├── ligand.rst7
│   ├── tleap_0.in
│   ├── tleap_0.log
│   ├── tleap.in
│   └── tleap.log
└── receptor.pdb

9 directories, 62 files
```



