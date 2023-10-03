import subprocess

def ddg_calc(muts_file:str,pdb_file:str):

    #open text file into a string
    text_file = open(muts_file, "r")
    mut_string = text_file.read()
    text_file.close()
    
    verbose= True

    cmd=['foldx','--command=PositionScan','--pdb='+pdb_file,'--positions='+mut_string]

    
    popen = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True)

    for line in popen.stdout:
        if verbose:
            print(line,end="",flush=True)
        
    # Check for success
    return_code = popen.wait()
    if return_code != 0:
        err = "Program failed.\n"
        raise RuntimeError(err)