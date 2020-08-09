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

from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from typing import List, Dict

# local imports

import config as cfg
import protein as prot

# -----------------------------------------------------------------------------

# FILE RELATED FUNCTIONS

def pickle_method(filename:str, method:str, context=''):
  """
  Pickles/unpickles a file to a binary file.
  """
  if method == 'wb':
      return pickle.dump(context, open(filename, method))
  elif method == 'rb':
      return pickle.load(open(filename, method))
    
    
def write_to_file(filename:str, context, function) -> None:
  """
  Writes to file.
  """
  with open(filename, 'w') as f:
    function(f, context)
  f.close()
  
  
def directify_name(name:str) -> str:
    """
    In case name has spaces, returns string with directory divider
    """
    return str(cfg.DIR_DIV.join(name.split()))

def add_fname_ext(fname:str, ext:str) -> str:
    """
    Appends extension to a filename.
    """
    return fname + ext
  

def df_to_csv(df:pd.DataFrame, filename:str, sep:str) -> pd.DataFrame:
  """
  Converts a DataFrame to a .csv file.
  """
  return df.to_csv(filename, sep, encoding='utf-8')


def cvs_to_df(filename:str, sep:str='.', col_idx:int=0) -> pd.DataFrame:
  """
  Converts a .csv file to a DataFrame.
  """
  return pd.read_csv(filename, sep=sep, index_col=col_idx, encoding='utf-8')


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