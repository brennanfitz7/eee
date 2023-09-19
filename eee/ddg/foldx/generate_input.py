##Generate input could make a list of all mutations that can be put on the command line later in the ddgcalc point
from eee.ddg import make_mutation_file
from eee.data import AA_3TO1

def generate_input(pdb_file:str,just_a_test=True):
    if just_a_test==True:
        mut_file_long=make_mutation_file(pdb_file).iloc[0:10]
    elif just_a_test==False:
        mut_file_long=make_mutation_file(pdb_file)
    else:
        print('error in just_a_test arg in generate input')

    #making a list with all the different mutations    
    mut_file=mut_file_long.drop_duplicates(subset=["res_pos"], keep='first')
    mutation_col=[]
    for index, row in mut_file.iterrows():
        mutation_col.append(AA_3TO1.get(row[1])+row[0]+row[2]+'a')

    #making a string with items separated by commas
    mut_string=','.join(mut for mut in mutation_col)

    #writing this string out to a txt file
    txt_file=pdb_file.split('.')[0]+'.txt'
    mut_txt = open(txt_file, "w")
    n = mut_txt.write(mut_string)
    mut_txt.close()


