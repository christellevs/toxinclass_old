import Bio as bio
import csv
import importlib
import io
import numpy as np
import os
import pandas as pd
import pickle
import random
import re
import statistics

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from typing import List, Dict

# local imports
import config as cfg
import protein as prot

# -----------------------------------------------------------------------------

# FILE RELATED FUNCTIONS

def pickle_method(fname:str, method:str, context=''):
  """
  Pickles/unpickles a file to a binary file.
  """
  if method == 'wb':
      return pickle.dump(context, open(fname, method))
  elif method == 'rb':
      return pickle.load(open(fname, method))
    
    
def write_to_file(fname:str, context, function) -> None:
  """
  Writes to file.
  """
  with open(fname, 'w') as f:
    function(f, context)
  f.close()
  

def df_to_csv(df:pd.DataFrame, fname:str, sep:str) -> pd.DataFrame:
  """
  Converts a DataFrame to a .csv file.
  """
  return df.to_csv(fname, sep, encoding='utf-8')


def cvs_to_df(fname:str, col_idx:int):
  """
  Converts a .csv file to a DataFrame.
  """
  return pd.read_csv(fname, index_col=col_idx, encoding='utf-8')


def split_seq(sequence:str) -> List[str]:
  """
  Splits string qeuence returns list.
  """
  return [char for char in sequence]


# PRE-PROCESSING FUNCIONS
# ----------------------------------------------------------

def crop_sequences(proteins: prot.Protein) -> List[prot.Protein]:
  """
  Returns protein sequences that are less than the specified length.
  """
  return [p for p in proteins if p.length <= cfg.MAX_SEQ_LEN]


def parse_fasta(path_fasta:str, is_toxic:int) -> List[str]:
  """
  Parses .fasta files into a list of Protein objects.
  """
  sequences = []
  with open(path_fasta) as fasta_file:
    for title, sequence in SimpleFastaParser(fasta_file):
      if not any(ele in sequence for ele in cfg.INVALID_AMINO):
        sequences.append( prot.Protein(title.split(None, 1)[0],
                                        is_toxic, len(sequence), split_seq(sequence)) )
  return sequences
    


# SPLTITING FUNCIONS
# -----------------------------------------------------------------------------

def get_kfold_splits(fname:str, x_features, y_labels):
  """
  Returns stratified shuffle cross validation splits from training dataset.
  """
  fold_dict = {}
  i = 1
  sss = StratifiedShuffleSplit(n_splits=cfg.K_FOLDS, test_size=cfg.VAL_SIZE,
                               random_state=cfg.RANDOM_SEED)
  print(sss)
  print(f'Number of k-fold splits: {sss.get_n_splits()}')
  for train_idx, val_idx in sss.split(x_features, y_labels):
    x_train, x_val = x_features.iloc[train_idx], x_features.iloc[val_idx]
    y_train, y_val = y_labels.iloc[train_idx], y_labels.iloc[val_idx]
    fold_dict[str(i)] = [x_train, x_val, y_train, y_val]
    i += 1
  pickle_method(fname, 'wb', fold_dict)
  return fold_dict


# -----------------------------------------------------------------------------