import pandas as pd

def get_ddg_diff(ens_ddg_file:str,save_file=True):

    """
    This function finds the standard deviation for each mutation across different structures

    ens_ddg_file : str
        a file of ddg values for two structures in the ensemble
    
    save_file : bool(True)
        An argument that determines whether the std deviation csv outfile is saved. Set to True as default.
    """
    
    df=pd.read_csv(ens_ddg_file)
    
    non_struct_columns=['site','mut','unfolded']
    struct_columns=[]
    

    for column in df.columns:
        if column not in non_struct_columns:
            struct_columns.append(column)
            
    #create pdb1 and pdb2 columns and rename the ddg columns
    df['pdb1']=struct_columns[0]
    df['pdb2']=struct_columns[1]
    
    df.rename(columns={struct_columns[0]:'ddg1',struct_columns[1]:'ddg2'}, inplace=True)
    
    #get diff in ddg column
    df['diff_in_ddg']=abs(df.ddg1-df.ddg2)
    
    #get other info columns
    wt_res=[i[0] for i in df.mut.to_list()]
    df['wt_res']=wt_res
    
    mut_res=[i[-1] for i in df.mut.to_list()]
    df['mut_res']=mut_res

    df=df.drop('unfolded',axis=1)
     
    out_file=ens_ddg_file.split('.')[0]+'_least_sim_diff.csv'
    
    if save_file == True:
        df.to_csv(out_file,index=False)
        
    return df

    
