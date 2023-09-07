
import subprocess


def ddg_calc(pdb_file:str):
    """
    Runs the calculation for ddg using ACDC-NN.

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest
        
    Returns
    -------
    None
    """
    pdb_id=pdb_file.split('.')[0]
    tsv_file=pdb_id+'_ddg_input.tsv'
    output_file=pdb_id+'_ddg_output.txt'
    
    
    verbose= True

    cmd=['acdc-nn','batch',tsv_file]

    
    popen = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True)
    
    f = open(output_file, 'w')
    
    for line in popen.stdout:
        f.write(line)
    f.close()

    for line in popen.stdout:
        if verbose:
            print(line,end="",flush=True)
        
    # Check for success
    return_code = popen.wait()
    if return_code != 0:
        err = "Program failed.\n"
        raise RuntimeError(err)