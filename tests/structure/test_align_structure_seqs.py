
import pytest

from eee.structure.align_structure_seqs import _run_muscle
from eee.structure.align_structure_seqs import align_structure_seqs
from eee.io.read_structure import read_structure
from eee.core.data import AA_1TO3

import numpy as np

import os
import glob
import shutil

def test__run_muscle():
    pass

def test_align_structure_seqs(structure_ensembles,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    for k in structure_ensembles:
        
        # Load the specified structure files
        dfs = []
        for f in structure_ensembles[k]:
            dfs.append(read_structure(f))
        
        # Get lists of all CA atoms and residues
        seqs = []
        for df in dfs:
            mask = np.logical_and(df.atom == "CA",df["class"] == "ATOM")
            this_df = df.loc[mask,:]
            seqs.append(list(this_df["resid"]))

        aligned_dfs = align_structure_seqs(dfs,keep_temporary=True)

        tmp_dir = list(glob.glob("calculation_*"))[0]


        # Make sure right output is being created
        assert len(aligned_dfs) == len(dfs)
        all_rows = np.sum([len(out_df) for out_df in aligned_dfs])
        for out_df in aligned_dfs:
            assert "shared_fx" in out_df.columns
            assert np.min(out_df["shared_fx"]) >= 0
            assert np.max(out_df["shared_fx"]) <= 1

        assert len(glob.glob(os.path.join(tmp_dir,"*.fasta"))) == 2
        
        # Make sure the input fasta is correct
        input_fasta = os.path.join(tmp_dir,"tmp-align_input.fasta")
        input_seqs = []
        with open(input_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    input_seqs.append([])
                else:
                    input_seqs[-1].extend(list(line.strip()))
        
        assert len(input_seqs) == len(structure_ensembles[k])
        for i in range(len(input_seqs)):
            assert np.array_equal([AA_1TO3[s] for s in input_seqs[i]],
                                seqs[i])

        # Make sure the output fasta is reasonable
        output_fasta = os.path.join(tmp_dir,"tmp-align_output.fasta")
        output_seqs = []
        with open(output_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    output_seqs.append([])
                else:
                    output_seqs[-1].extend(list(line.strip()))

        # Only one length allowed; same as length of column_contents
        lengths = np.unique([len(s) for s in output_seqs])
        assert len(lengths) == 1

        # make sure muscle_binary is interpreted correctly
        with pytest.raises(RuntimeError):
            align_structure_seqs(dfs,muscle_binary="not_real")

        # Delete temporary files
        for f in glob.glob("calculation_*"):
            shutil.rmtree(f)

        # Make sure keep_temporary flag is interpreted correctly
        #column_contents, column_indexes 
        aligned_dfs= align_structure_seqs(dfs,keep_temporary=False)
        assert len(glob.glob("calculation_*")) == 0

    os.chdir(current_dir)