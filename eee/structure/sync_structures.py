"""
Functions for taking raw RCSB output with several structures and creating input
for an EEE calculation. 
"""

from eee.io.write_pdb import write_pdb
from eee.io.read_structure import read_structure
from eee._private import logger
from eee.structure.clean_structure import clean_structure
from eee.structure.align_structure_seqs import align_structure_seqs
from eee.structure.align_structures import align_structures

import os
import glob
import shutil

def _create_unique_filenames(files):
    """
    This wacky block of code trims back filenames, right to left, until they
    are unique. This solves edge case where someone puts in files with same
    name from different directory (like 1stn.pdb and ../test/1stn.pdb). This
    loop would create output files "1stn.pdb" and "test__1stn.pdb". 
    """

    found_filenames = False
    counter = -1
    while not found_filenames:
        name_mapper = []

        if len(list(set(files))) != len(files):
            err = "structure_files must have unique filenames!"
            raise ValueError(err)

        for i in range(len(files)):    
            real_path = "__".join(files[i].split(os.path.sep)[counter:])
            if real_path not in name_mapper:
                name_mapper.append(real_path)
                found_filenames = True
            else:
                counter -= 1
                found_filenames = False
                break

    name_mapper = dict([(files[i],name_mapper[i]) for i in range(len(files))])

    return name_mapper


def sync_structures(structure_files,
                    out_dir,
                    overwrite=False,
                    verbose=False,
                    keep_temporary=False):
    """
    Take a set of structures, clean up, align, and figure out which sites are
    shared among all structures. Output is a directory with pdb files and a
    report describing structures. The residue numbers are replaced with their 
    sites in the alignment (meaning residue numbers compare between structures).
    The b-factor column of each pdb file has the fraction of structures in which
    that specific site is seen. The occupancy column is 1 if the amino acids are
    same at the site for all structures, 0 if the amino acids are different. 
    (Note: at sites with a mix of Cys and Ser across structures, the Ser 
    residues are mutated to Cys). HETATM entries will always have 0 occupancy
    and b-factors. 

    Parameters
    ----------
    structure_files : list
        list of structure files to use for the calculation. These files should 
        be in RCSB cif (preferred) or pdb format.
    out_dir : str
        output directory for the cleaned up files in pdb format. This directory
        should either not exist or be empty. 
    overwrite : bool, default=False
        overwrite an existing output directory
    verbose : bool, default=False
        write out all output to standard output
    keep_temporary : bool, default=False
        do not delete temporary files
    """
    
    # See if the output directory exists
    exists = False
    if os.path.exists(out_dir):
        if os.path.isdir(out_dir):
            if len(glob.glob(os.path.join(out_dir,"*"))) > 0:
                exists = True
        else:
            exists = True

    if exists:
        if not overwrite:
            err = f"output directory {out_dir} already exists.\n"
            raise FileExistsError(err)
        else:
            shutil.rmtree(out_dir)
    
    # Make new directory. 
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # Load the specified structure files
    dfs = []
    for f in structure_files:
        dfs.append(read_structure(f))

    # Clean up structures --> build missing atoms or delete residues with
    # missing backbone atoms. 
    logger.log("Cleaning up structures with FoldX.")
    
    for i in range(len(dfs)):
        dfs[i] = clean_structure(dfs[i],
                                 verbose=verbose,
                                 keep_temporary=keep_temporary)

    # Figure out which residues are shared between what structures
    logger.log("Aligning sequences using muscle.")
    dfs = align_structure_seqs(dfs,verbose=verbose,keep_temporary=keep_temporary)

    # Align structures in 3D
    logger.log("Aligning structures using lovoalign.")
    dfs = align_structures(dfs,verbose=verbose,keep_temporary=keep_temporary)

    # Create a unique output name for each structure file
    name_mapper = _create_unique_filenames(structure_files)

    # Write out file names. 
    logger.log(f"Writing output to '{out_dir}'.")
    for i in range(len(structure_files)):

        f = f"{name_mapper[structure_files[i]]}_clean.pdb"
        f = os.path.join(out_dir,f)
        
        write_pdb(dfs[i],
                  f,
                  bfactor_column="shared_fx",
                  occ_column="identical_aa")
    

    return dfs
    