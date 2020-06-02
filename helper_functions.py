import csv
import pickle
import pandas as pd


# writes or reads data to/from a binary file
def pickle_method(fname, method, context):
    if method == 'wb':
        return pickle.dump(context, open(fname, method))
    elif method == 'rb':
        return pickle.load(open(fname, method))

# saves df to a csv file
def df_to_csv(df, fname, sep):
  df.to_csv(fname, sep, encoding='utf-8')

# read in csv, returns df
def cvs_to_df(fname, col_idx):
  return pd.read_csv(fname, index_col=col_idx, encoding='utf-8')

# splits sequence, returns a list
def split_seq(sequence):
  return [char for char in sequence]