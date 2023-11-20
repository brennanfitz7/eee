from eee.data import AA_3TO1

import pandas as pd
import os


def test_convert_to_df(test_foldx_ddgs):
    for ddg_output in test_foldx_ddgs:
        pdb_id=ddg_output.split('_')[1]

        ddg_df_raw=pd.read_csv(ddg_output, sep='\t')

        column_names=['Mutation','DDG']
        ddg_df=ddg_df_raw.set_axis(column_names, axis=1) 

        #test that the dataframe has these two columns and they have the same length
        assert 'Mutation' in ddg_df.columns
        assert 'DDG' in ddg_df.columns
        assert len(ddg_df['Mutation'])==len(ddg_df['DDG'])

        old_mut_list=ddg_df['Mutation'].tolist()
        new_mut_list=[]
        old_ddg_list=ddg_df['DDG'].tolist()
        new_ddg_list=[]

        #make sure both of these are lists
        assert issubclass(type(old_mut_list),list)
        assert issubclass(type(old_ddg_list),list)
        assert len(old_mut_list)==len(old_ddg_list)

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

        assert len(new_mut_list)==len(new_ddg_list) 
        #old mut list should include all mutations but new mut list will cut mutating to self
        #I think this will work but maybe it won't?
        #adding 1 to old mut list bc iirc the first mutation gets cut off when it's converted from csv to df
        #this is fine because first mutation is a mutate to self situation
        assert len(new_mut_list)/19==(len(old_mut_list)+1)/20

        #test that final character in new_muts is a capital letter
        for item in new_mut_list:
            assert item[-1].isupper()
            assert item[-1].isalpha()


        df_data={'Mutation':new_mut_list, 'DDG': new_ddg_list}
        
        ddg_df=pd.DataFrame(df_data)

        #test that the dataframe index length is correct
        assert len(ddg_df.index)==len(new_mut_list)
        assert len(ddg_df.index)==len(new_ddg_list)

        output_df=pdb_id.split('.')[0]+'_foldx_ddg_df.csv'
        #this creates the output_df wherever the program is run
        ddg_df.to_csv(output_df, index=False)

        assert os.path.isfile(output_df)==True
    

        