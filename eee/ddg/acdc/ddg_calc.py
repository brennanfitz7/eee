
import subprocess


def ddg_calc(muts_file:str):
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
    pdb_id=muts_file.split('_')[0]
    output_file=pdb_id+'_acdc_raw_ddgs.txt'
    
    
    verbose= True

    cmd=['acdc-nn','batch',muts_file]

    
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