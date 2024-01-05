from eee.ddg.make_mutation_file import make_mutation_file
from eee.core.data import AA_3TO1


def generate_input(pdb_file:str, remove_multiple_models:bool,just_a_test=True):
    if remove_multiple_models==True:
        mut_file_long=make_mutation_file(pdb_file,remove_multiple_models=True)
    if remove_multiple_models==False:
        mut_file_long=make_mutation_file(pdb_file,remove_multiple_models=False)
    else: 
        print('No arg for remove_multiple_models. Must enter True or False')

    if just_a_test==True:
        mut_file_long=mut_file_long.iloc[0:10]

    #making a list with all the different mutations    
    mut_file=mut_file_long.drop_duplicates(subset=["res_pos"], keep='first')
    mutation_col=[]
    for index, row in mut_file.iterrows():
        mutation_col.append(AA_3TO1.get(row[1])+row[0]+row[2]+'a')

    #making a string with items separated by commas
    mut_string=','.join(mut for mut in mutation_col)

    #writing this string out to a txt file
    txt_file=pdb_file.split('.')[0]+'_foldx_muts.txt'
    mut_txt = open(txt_file, "w")
    mut_txt.write(mut_string)
    mut_txt.close()


