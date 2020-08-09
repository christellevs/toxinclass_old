import Bio as bio
import pandas as pd

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from typing import List

# local imports
import config as cfg
import protein as prot
import utils as u

# -----------------------------------------------------------------------------


def parse_fasta(path_fasta:str, is_toxic:int) -> List[str]:
    """
    Parses .fasta files into a list of Protein objects.
    """
    sequences = []
    with open(path_fasta) as fasta_file:
        for title, sequence in SimpleFastaParser(fasta_file):
            if not any(ele in sequence for ele in cfg.INVALID_AMINO):
                sequences.append( prot.Protein(title.split(None, 1)[0],
                                                is_toxic, len(sequence), sequence.split()))
    return sequences


def crop_sequences(proteins:prot.Protein) -> List[prot.Protein]:
    """
    Returns protein sequences that are less than the specified length.
    """
    return [p for p in proteins if p.length <= cfg.MAX_SEQ_LEN]


    # MATRICES VALUES
# ---------------------------------------

def _get_matrix_values(self, seq_dict):
    """
    Returns matrix from input dictionary.
    """
    return np.array([seq_dict[i] for i in seq_dict.keys()])


def update_matrices(self, protein_seq_list):
    """
    Updates both matrices.
    """
    for protein in protein_seq_list:
        protein.matrix_raw = get_matrix_values(protein.seq_dict_raw)
        protein.matrix_diff = get_matrix_values(protein.seq_dict_diff)

def append_proteins(proteins):
    """
    Appends Atchley values to the Protein object matrices.
    """
    append_atchley_values(proteins)
    update_matrices(proteins)
    u.pickle_method(cfg.f_train_proteins, 'wb', proteins)


def proteins_to_df(proteins:prot.Protein) -> pd.DataFrame:
    """
    Returns a DataFrame from a list of proteins.
    """
    df = pd.DataFrame.from_records([p.to_dict() for p in proteins])
    u.pickle_method(cfg.f_train_df, 'wb', df)
    return df


def get_atchley_values_list(aminos:List[str], idx:int) -> List[int]:
    """
    Returns all atchley values in a list for a specific amino acid.
    """
    
    return [float(cfg.DICT_ATCHLEY.get(i)[idx]) for i in aminos]


def get_atchley_values_raw(sequence, seq_dict_r):
    """
    Returns raw atchley values in a dictionary.
    """
    seq_dict_r['f1'] = get_atchley_values_list(sequence, 0)
    seq_dict_r['f2'] = get_atchley_values_list(sequence, 1)
    seq_dict_r['f3'] = get_atchley_values_list(sequence, 2)
    seq_dict_r['f4'] = get_atchley_values_list(sequence, 3)
    seq_dict_r['f5'] = get_atchley_values_list(sequence, 4)

# -----------------------------------------------

def get_change_list(atchley_list):
    """
    Calculates sequential change for single atchley value.
    """
    atchley_list.insert(0, 0)
    change_list = [i for i in (np.diff(atchley_list))]
    atchley_list.pop(0)
    return change_list


def get_atchley_diff(seq_dict_r, seq_dict_d):
    """
    Returns the sequential change for each atchley value as a dictionary.
    """
    seq_dict_d['f1_d'] = get_change_list(seq_dict_r.get('f1'))
    seq_dict_d['f2_d'] = get_change_list(seq_dict_r.get('f2'))
    seq_dict_d['f3_d'] = get_change_list(seq_dict_r.get('f3'))
    seq_dict_d['f4_d'] = get_change_list(seq_dict_r.get('f4'))
    seq_dict_d['f5_d'] = get_change_list(seq_dict_r.get('f5'))


def append_atchley_values(proteins_list):
    """
    Appends the atchley values to the dictionary of a ProteinSequence object.
    """
    for protein in proteins_list:
        get_atchley_values_raw(protein.sequence, protein.seq_dict_raw)
        get_atchley_diff(protein.seq_dict_raw, protein.seq_dict_diff)
    
# def split_seq(sequence:str) -> List[str]:
#     """
#     Splits string qeuence returns list.
#     """
#     return [char for char in sequence]    

# -----------------------------------------------------------------------------