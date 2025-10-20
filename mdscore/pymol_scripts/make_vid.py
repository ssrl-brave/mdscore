


def parse_arguments():
    from argparse import ArgumentParser
    """Sets up and parses command line arguments."""
    parser = ArgumentParser(
        description="Processes a molecular dynamics trajectory file for PyMOL visualization.",
    )

    parser.add_argument(
        '--topfile',
        type=str,
        required=True,
        help="Path to the topology file (e.g., .prmtop or .pdb)."
    )

    parser.add_argument(
        '--ncfile',
        type=str,
        required=True,
        help="Path to the trajectory file (e.g., .nc, .dcd, or .xtc)."
    )

    parser.add_argument(
        '--prefix',
        type=str,
        required=True,
        help="Basename for output files (e.g., 'system' creates system.dcd and system.pdb)."
    )

    return parser.parse_args()


def main(ncfile, topfile, prefix, ligand_color="magenta"):
    import mdtraj
    import os
    from pymol import cmd
    import tempfile
    import shutil

    cmd.reinitialize()

    temp_dir = tempfile.mkdtemp()
    temp_filename = os.path.basename(topfile) + ".prmtop"
    topfile_path = os.path.join(temp_dir, temp_filename)
    shutil.copy(topfile, topfile_path)
    print(f"Temporary .prmtop created at: {topfile_path}")
    md = mdtraj.load(ncfile, top=topfile_path)

    dcdf = prefix+".dcd"
    pdbf = prefix + ".pdb"
    md.save_dcd(dcdf)
    md[0].save_pdb(pdbf)

    cmd.load(pdbf, "complex")
    cmd.load_traj(dcdf, "complex")
    cmd.color("green", "polymer")
    cmd.intra_fit("(polymer)", 1)
    cmd.color(ligand_color, "organic")
    cmd.create("ligand", "organic", 1)
    #cmd.orient("ligand")
    cmd.zoom("ligand")
    cmd.move("z", -70)
    cmd.clip("slab", 70)
    cmd.delete("ligand")
    cmd.mset("1-2000")
    cmd.set("ray_trace_frames", 0)
    cmd.mclear()
    cmd.set("antialias",1)

    vidprefix = os.path.join(temp_dir, "traj")
    cmd.mpng(vidprefix)

    vidname = prefix + '.mp4'
    os.system(f"ffmpeg -framerate 30 -i {vidprefix}%04d.png -c:v libx264 -pix_fmt yuv420p {vidname} -y")
    print(f"Done. Video written to {vidname}")


if __name__=="__main__":
    args = parse_arguments()
    main(args.ncfile, args.topfile, args.prefix)

elif __name__=="make_vid_py":
    # envoked using the following command within Pymol GUI:
    # run make_vid.py, module
    import os
    assert "NCFILE" in os.environ
    assert "TOPFILE" in os.environ
    assert "VIDPREFIX" in os.environ
    main(os.environ["NCFILE"], os.environ["TOPFILE"], os.environ["VIDPREFIX"])

elif __name__=="pymol":
    # envoked using the following command within Pymol GUI:
    # run make_vid.py, local
    import os
    #ncnames = open("nc_fnames.txt", "r").readlines()
    ncnames = open("fnames2.txt", "r").readlines()
    for i, ncname in enumerate(ncnames):
        ncname = ncname.strip()
        topname = ncname.replace(".nc", ".top")
        assert os.path.exists(topname)
        #if "hit" in ncname:
        #    ligand_color="magenta"
        #elif "miss" in ncname:
        #    ligand_color="cyan"

        #else:
        #    ligand_color="green"
        ligand_color="magenta"
        prefix=os.path.join(os.path.dirname(ncname), "vid")
        main(ncname, topname, prefix, ligand_color)
