import pandas as pd

def get_ddg_std_dev(ens_ddg_file:str,save_file=True):

    """
    This function finds the standard deviation for each mutation across different structures

    ens_ddg_file : str
        a file of ddg values for all structures in the ensemble
    
    save_file : bool(True)
        An argument that determines whether the std deviation csv outfile is saved. Set to True as default.
    """
    
    df=pd.read_csv(ens_ddg_file)
    
    non_struct_columns=['site','mut','unfolded']
    struct_columns=[]
    

    for column in df.columns:
        if column not in non_struct_columns:
            struct_columns.append(column)
            
    std_devs=df[struct_columns].std(axis=1)
    df['std_devs']=std_devs

    non_struct_columns.append('std_devs')

    for column in df.columns:
        if column not in non_struct_columns:
            df=df.drop(column,axis=1)
    
    df=df.drop('unfolded',axis=1)
     
    out_file=ens_ddg_file.split('.')[0]+'_std_dev.csv'
    
    if save_file == True:
        df.to_csv(out_file,index=False)
        
    return df

    
