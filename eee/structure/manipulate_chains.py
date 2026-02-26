import pandas as pd

from Bio import Align

from eee.core.data import AA_3TO1

def chain_reindex(df,prev_chain,new_chain,pdb=None,write_pdb=False):
    
    #first change chain names
    
    #if there is a previous chain with the new chain name I want, change it to a placeholder chain
    df.loc[df.chain == new_chain, 'chain'] = 'ph'
    
    #now change my prev_chain residues to new_chain
    df.loc[df.chain == prev_chain,'chain'] = new_chain
    
    #change the placeholder chain to the prev_chain (keeping it simple)
    df.loc[df.chain == 'ph','chain'] = prev_chain
    
    #then sort by chain (sort by atom number secondarily so as to not mess with the order of the atoms)
    
    df=df.sort_values(['chain', 'atom_num'], ascending=[True, True])
    
    if write_pdb==True:
        write_pdb(df, pdb, overwrite=True,write_with_models=True)
    
    return df

def get_chain_df(df):
    #select for atoms and drop duplicates
    redundant_prot_seq=df[df['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')

    #getting rid of redundant chains
    #making a dataframe of Chain ID and sequence
    chain_df=pd.DataFrame(columns=['Chain_ID','seq'])
    for item in redundant_prot_seq.chain.unique():
        chain=str(item)
        temp_list=[]
        for index,row in redundant_prot_seq.iterrows():
            if chain==row.iloc[0]:
                temp_list.append((AA_3TO1.get(row.iloc[1])))

        test_string=''.join(temp_list)
        chain_df.loc[len(chain_df.index)]=[chain,test_string]

    return chain_df


def get_unique_chains(df):

    #get df made up of chains
    chain_df=get_chain_df(df)
        
    #finding original chains 
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    #now time to do alignment between each different chain
    redundant=[]
    unique_chains=[]
    for i in range(0,len(chain_df)):

        #this line is to make sure we don't run chains that have already been grouped with other chains
        if chain_df.Chain_ID[i] not in redundant:
            main_chain=chain_df.Chain_ID[i]
            #now we do a pairwise alignment of this chain with every other chain (that hasn't been grouped already)
            for n in range(i+1,len(chain_df)):
                score=aligner.score(chain_df.seq[i],chain_df.seq[n])
                match=score/min(len(chain_df.seq[i]),len(chain_df.seq[n]))
                #if the match is sufficient we append the chain to the similar seqs list
                if match >= 0.95:
                    redundant.append(chain_df.Chain_ID[n])
            unique_chains.append(main_chain)


    return unique_chains