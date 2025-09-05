"""
Write a pdb file from a dataframe with atom information.
"""

import os

def write_pdb(df,
              pdb_file,
              overwrite=False,
              bfactor_column=None,
              occ_column=None,
              write_with_models=False):
    """
    Write a pdb file given a pandas dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe with structural data (generally created using read_structure)
    pdb_file : str
        name of pdb file to write
    overwrite : bool, default=False
        overwrite the pdb file if it exists
    bfactor_column : str, optional
        name of a column in df to use for the bfactor data in the pdb
    occ_column : str, optional 
        name of a column in df to use for the occupancy data in the pdb
    """
    
    if os.path.exists(pdb_file):
        if not overwrite:
            err = f"pdb_file {pdb_file} already exists.\n"
            raise FileExistsError(err)
        else:
            if os.path.isfile(pdb_file):
                os.remove(pdb_file)
            else:
                err = f"pdb_file {pdb_file} exists but is not a regular file.\n"
                err += "Cannot overwrite.\n"
                raise FileExistsError(err)

    df = df.copy()
    
    last_chain = None
    last_class = None
    last_model=None
    with open(pdb_file,'w') as f:

        if write_with_models==True:
            f.write('MODEL        1\n')

        counter = 1
        for i in df.index:
            
            row = df.loc[i,:]
            
            chain = row['chain']
            if last_chain is None:
                last_chain = chain

            model=row['model']
            if last_model is None:
                last_model = model
                
            atom_class = row['class']
            if last_class is None:
                last_class = atom_class
            
            if chain != last_chain and last_class == "ATOM":
                f.write("TER\n")
                last_chain = chain 

            if write_with_models==True and model != last_model:
                f.write("ENDMDL\nMODEL        ")
                f.write(str(model))
                f.write('\n')
                last_model = model
                
            f.write(f"{row['class']:6s}{counter:5d}")
        
            atom = row["atom"]
            if len(atom) < 5:
                atom = f" {atom:<4s}"

            resid_num=str(row['resid_num'])
            if resid_num.isdigit():
                resid_num=f"{resid_num:>4s} "
            else:
                resid_num=f"{resid_num:>5s}"

            f.write(f" {atom}{row['resid']:3s} {row['chain']}{resid_num}")
            f.write(f"   {row['x']:8.3f}{row['y']:8.3f}{row['z']:8.3f}")

            if bfactor_column is None:
                b = row["b"]
            else:
                b = row[bfactor_column]
            
            if occ_column is None:
                occ = row["occ"]
            else:
                occ = row[occ_column]


            f.write(f"{occ:6.2f}{b:6.2f}{row['elem']:>12s}\n")
        
            last_class = atom_class
            counter += 1
            
        f.write("END\n")
