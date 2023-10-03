from eee.data import AA_3TO1

import pandas as pd

def convert_to_df(ddg_output:str):
#figure out how to add header without deleting first row--it's okay for now bc the first line is mutating the residue to itself
    #foldx is able to identify the protein tag alone--no path
    pdb_id=ddg_output.split('_')[1]
    #establishing the name of the output file to find--should be where I ran it regardless of where pdb files are
    
    ddg_df_raw=pd.read_csv(ddg_output, sep='\t')
    
    column_names=['Mutation','DDG']
    ddg_df=ddg_df_raw.set_axis(column_names, axis=1)
    
    #going to actually make new lists rather than deleting from old one to keep things from getting sketchy with indices
    old_mut_list=ddg_df['Mutation'].tolist()
    new_mut_list=[]
    old_ddg_list=ddg_df['DDG'].tolist()
    new_ddg_list=[]
    
    for item in old_mut_list:
        one_letter_code=AA_3TO1.get(item[0:3])
        #for mutations to a charged form of histidine
        if item[-1] in 'oef' and one_letter_code!='H':
            new_mut_list.append(one_letter_code+item[4:len(item)-1]+'H')
            new_ddg_list.append(old_ddg_list[old_mut_list.index(item)])
        #for regular mutations that are not mutate to self
        elif one_letter_code!=item[-1] and item[-1].isupper()==True:
            new_mut_list.append(one_letter_code+item[4:len(item)])
            new_ddg_list.append(old_ddg_list[old_mut_list.index(item)])
        #in case there is another weird foldx lowercase residue
        elif item[-1].isupper==False and item[-1] not in 'oef':
            err=f'Unrecognized residue. Consult FoldX residues.\n'
            raise ValueError(err)
        else:
            continue
        
    df_data={'Mutation':new_mut_list, 'DDG': new_ddg_list}
        
    ddg_df=pd.DataFrame(df_data)
    
    output_df=pdb_id.split('.')[0]+'_foldx_ddg_df.csv'
    #this creates the output_df wherever the program is run
    ddg_df.to_csv(output_df, index=False)
    
    return ddg_df