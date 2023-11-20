from eee.structure.sync_structures import sync_structures
from eee.ddg import acdc
from eee.ddg import foldx
#cannot test without test csv files--create with test_pdbs

import pandas as pd
import glob
import shutil
import os


def test_run_ensemble(test_csvs):
   
   module_list=['acdc','foldx'] 
   for module in module_list:   
      for csv in test_csvs:
         pdb_df=pd.read_csv(csv)
         assert issubclass(type(pdb_df),pd.DataFrame)

         pdb_list=pdb_df['PDB']
         name_list=pdb_df['NAME']
         name_dict = dict(zip(pdb_list, name_list))
         assert issubclass(type(name_dict),dict)

         module_dict={'acdc':acdc,'foldx':foldx}
         calculator=module_dict.get(module)

         if calculator==acdc:
            acdc_info=pdb_df['ACDC'].tolist()
            assert len(acdc_info)==2
            hhblits_path=acdc_info[0]
            uniref_path=acdc_info[1]
    

         prot_name=csv.split('.')[0]
         sync_structures(structure_files=pdb_list, out_dir=prot_name)
         synced_pdbs=glob.glob('*.pdb',root_dir=prot_name)
         assert issubclass(type(synced_pdbs), list)

         if calculator==foldx:
            for pdb in synced_pdbs: 
               shutil.move(prot_name+'/'+pdb,pdb)
         
         #if using foldx, test that pdbs were moved successfully
         if calculator==foldx:
            test_move=glob.glob('*.pdb')
            for pdb in synced_pdbs:
               assert pdb in test_move==True

         #test just_a_test? how do I do that? tests seem like they don't actually run the program but maybe they do
         #have all the tests do just_a_test otherwise it will take too long--ask Mike about this
         
         df_list=[]

         for pdb in synced_pdbs:
            if calculator==foldx:
               pdb_file=pdb
               pdb_id=pdb_file.split('.')[0]
               calculator.generate_input(pdb_file)
               assert os.path.isfile(pdb_file.split('.')[0]+'_foldx_muts.txt')==True
               calculator.ddg_calc(muts_file=pdb_id+'_'+module+'_muts.txt', pdb_file=pdb_file)   
               assert os.path.isfile('PS_'+pdb_file.split('/')[-1][0:-4]+'_scanning_output.txt')==True
               ddg_df=calculator.convert_to_df('PS_'+pdb_file[0:-4]+'_scanning_output.txt')

            if calculator==acdc:
               pdb_file=str(prot_name+'/'+pdb)
               pdb_id=pdb_file.split('.')[0]
               calculator.generate_input(pdb_file,hhblits_path=hhblits_path,uniref_path=uniref_path)
               assert os.path.isfile(pdb_file.split('.')[0]+'_acdc_muts.tsv')==True
               calculator.ddg_calc(pdb_id+'_'+module+'_muts.tsv')
               assert os.path.isfile(pdb_file.split('.')[0]+'_acdc_raw_ddgs.txt')==True
               ddg_df=calculator.convert_to_df(pdb_id+'_'+module+'_raw_ddgs.txt')


            assert issubclass(type(ddg_df),pd.DataFrame)
            ddg_df.name=pdb.split('_')[0]
            df_list.append(ddg_df)

         #test that there are the correct number of dfs in in the df_list
         assert len(df_list)==len(synced_pdbs)

         mut_sets=[]

         for df in df_list:
            mut_sets.append(set(df.Mutation))
         #test that all dfs were added to mut_sets
         assert len(df_list)==len(mut_sets)
         #test that all the items in mut_sets are sets
         for item in mut_sets:
            assert issubclass(type(item),set)==True

         shared_muts=mut_sets[0].intersection(*mut_sets[1:])
         #testing that all mutations in shared_muts are also in all of the mut sets
         for mut in shared_muts:
            for pdb_set in mut_sets:
               assert mut in pdb_set==True

         clean_df_list=[]
         for df in df_list:
            mask=df.Mutation.isin(shared_muts)
            clean_df=df.loc[mask,:]
            #testing clean_df and shared_muts have the same number of mutations
            assert len(clean_df.index)==len(shared_muts)
            clean_df.name=df.name
            clean_df_list.append(clean_df)

         #testing that all dfs made it into clean_df_list
         assert len(clean_df_list)==len(df_list)

         mut_source=clean_df_list[0]
         combined_df=mut_source[['Mutation']].copy()
         #test that mutation is now a column in combined_df and has the same length as the shared_muts length
         assert len(combined_df.Mutation)==len(shared_muts)

         for df in clean_df_list:
            DDG_col = df['DDG']
            #test that mask with shared muts worked on ddg
            assert len(DDG_col)==len(shared_muts)
            ##use name_dict to get the correct name for the pdb file
            name_of_df=name_dict.get(df.name)
            combined_df[name_of_df] = DDG_col

         #test that the combined df has the correct number of columns
         assert len(combined_df.columns)==len(synced_pdbs)+1

         combined_df.to_csv(prot_name+'_'+module+'_ddgs.csv', index=False)

         #test that the csv is produced
         assert os.path.isfile(prot_name+'_'+module+'_ddgs.csv', index=False)==True

         return(combined_df)
    
    
         