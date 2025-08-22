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
from eee.io.incorporate_models import incorporate_models
from eee.structure.reassign_chains import reassign_chains

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
                    all_models_necessary:bool,
                    align_seqs= False,
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
    align_seqs : bool
        align the sequences in all the structures. will add negative numbers to the residues, which can cause a problem with ThermoMPNN.
    overwrite : bool, default=False
        overwrite an existing output directory
    verbose : bool, default=False
        write out all output to standard output
    keep_temporary : bool, default=False
        do not delete temporary files
    all_models_necesary: bool
        if true, all models should be used in syncing
    """
    
    logger.log('Starting with '+str(len(structure_files))+' structure files.')

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

    # Load the specified structure files, incorporate models if needed, and add name of file into pdb
    raw_dfs = []

    if all_models_necessary==True:
        for f in structure_files:
            try:
                read_df=read_structure(f,remove_multiple_models=False)
                read_df['name']=str(f)
                raw_dfs.append(incorporate_models(read_df))
            except:
                logger.log('PDB '+str(f)+' failed when reading structure and incorporating models.')
        if all_models_necessary==False:
            for f in structure_files:
                try:
                    read_df=read_structure(f)
                    read_df['name']=str(f)
                    raw_dfs.append(read_df)
                except:
                    logger.log('PDB '+str(f)+' failed when reading structure.')


    
    # Clean up structures --> build missing atoms or delete residues with
    # missing backbone atoms. 
    #this shouldn't have issues with name column
    logger.log("Cleaning up structures with FoldX.")
    
    dfs=[]

    for i in range(len(raw_dfs)):

        try:
            dfs.append(clean_structure(raw_dfs[i],
                                verbose=verbose,
                                keep_temporary=keep_temporary,
                                name_in_df=True))
        except:
            failed_name=raw_dfs[i].name[0]
            logger.log('PDB '+failed_name+' failed foldx cleaning.')

            
    if len(dfs)!=len(structure_files):
        logger.log("After cleaning up structures with FoldX, there are "+str(len(dfs))+" protein dataframes.")

    #Align chains and make sure the same chains share the same chain IDs between structures
    logger.log("Changing chains to ensure chain IDs are the same between structures.")
    dfs = reassign_chains(dfs,ensemble=out_dir)


    if len(dfs)!=len(structure_files):
        logger.log("After chain reassignment, there are "+str(len(dfs))+" protein dataframes")

    # Figure out which residues are shared between what structures
    if align_seqs == True:
        logger.log("Aligning sequences using muscle.")
        dfs = align_structure_seqs(dfs,verbose=verbose,keep_temporary=keep_temporary)

        if len(dfs)!=len(structure_files):
            logger.log("After aligning seqs, there are "+str(len(dfs))+" protein dataframes.")


    # Align structures in 3D
    logger.log("Aligning structures using lovoalign.")
    dfs = align_structures(dfs,verbose=verbose,keep_temporary=keep_temporary)

    if len(dfs)!=len(structure_files):
        logger.log("After aligning structures with lovoalign, there are "+str(len(dfs))+" protein dataframes.")

        
    #create a list of files then create a unique output name for each structure file
    filenames=[]
    for i in range(len(dfs)):
        df=dfs[i]
        pdb_name=df.name[0]
        filenames.append(pdb_name)


    name_mapper = _create_unique_filenames(filenames)   

    # Write out file names. 
    logger.log(f"Writing output to '{out_dir}'.")
    for i in range(len(dfs)):

        f = f"{name_mapper[filenames[i]]}_clean.pdb"
        f = os.path.join(out_dir,f)
        
        if align_seqs == True:
            write_pdb(dfs[i],
                    f,
                    bfactor_column="shared_fx",
                    occ_column="identical_aa")
        if align_seqs == False:
            write_pdb(dfs[i],
                    f)


    return dfs
    