import io
import itertools as ite
import numpy as np
import os.path
import pandas as pd
import pprint as p
import re
import statistics

from sklearn.utils import resample
from tabulate import tabulate

import sys
sys.path.insert(1, 'modules')

# module imports
import config as cfg
import data_references as dr
import helper_funcs as hf
import visual_funcs as vf
import protein as prot

# -----------------------------------------------------------------------------

# PRE-PROCESSING TRAINING SET
# PARSING FASTA
toxic_proteins = hf.parse_fasta(cfg.f_train_toxic_fasta, 1)
atoxic_proteins = hf.parse_fasta(cfg.f_train_atoxic_fasta, 0)

# print test
print(f'Total toxic sequences: {len(toxic_proteins)}')
print(f'Total atoxic sequences: {len(atoxic_proteins)}')


# DOWNSAMPLING
# -----------------------------------------------------------------------------

toxic_proteins = hf.crop_sequences(toxic_proteins)
atoxic_proteins = hf.crop_sequences(atoxic_proteins)

total_toxic_seqs = len(toxic_proteins)
print(f'Total toxic sequences ( <= {cfg.MAX_SEQ_LEN}): {total_toxic_seqs}')
print(f'Total toxic sequences ( <= {cfg.MAX_SEQ_LEN}): {len(atoxic_proteins)}')


atoxic_proteins = resample(atoxic_proteins, replace=False, n_samples=total_toxic_seqs, random_state=cfg.RANDOM_SEED)
print(f'Total protein sequences in atoxic list post-downsampling: {len(atoxic_proteins)}')


# COMBINING & APPENDING
# -----------------------------------------------------------------------------

# TODO


print(f'Total overall protein sequences: {len(proteins)}')
print(proteins[0].matrix_diff)

df_proteins = hf.proteins_to_df(proteins)

# print test
print('\nProtein training set info:')
df_proteins.info()

print('\nChecking value counts for each class in df combined:\n1 == toxic\n0 == atoxic\n----------')
print(df_proteins['toxic'].value_counts())


# -----------------------------------------------------------------------------
