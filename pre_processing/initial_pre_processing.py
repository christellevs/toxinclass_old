#TOXIFY: https://github.com/tijeco/toxify
#ToxClassifier: https://github.com/rgacesa/ToxClassifier


# -----------------------------------------------------------------------------
# IMPORTS
import Bio as bio
import csv
import io
import itertools as ite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import os.path
import random
import re
import seaborn as sns
import sklearn
import statistics

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import OrderedDict

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.utils import resample


# # loading atchley csv data and turning to dict
# df_atchley = cvs_to_df(f_atchley, 0)
# df_atchley.rename(columns={'amino.acid': 'amino_acid'}, inplace=True)
# df_atchley.set_index('amino_acid', inplace=True)

# for col in df_atchley['f1': 'f5']:
#   df_atchley[col] = df_atchley[col].apply(lambda x: re.sub(r'[^\x00-\x7F]+','-', x)).astype(float)

# pickle_method(f_atchley_table, 'wb', df_atchley)

# dict_atchley = df_atchley.T.to_dict('list')
# pickle_method(f_atchley_dict, 'wb', dict_atchley)

# loading pickled Atchley dictionary
dict_atchley = pickle_method(f_atchley_dict, 'rb', '')
dict_atchley


# parsing toxic and atoxic data from fasta files
toxic_list = parse_fasta(f_toxic_fasta, 1)
atoxic_list = parse_fasta(f_atoxic_fasta, 0)

total_toxic_seqs = len(toxic_list)
print('Total toxic sequences: ', total_toxic_seqs)
print('Total atoxic sequences: ', len(atoxic_list))

"""# **4. Downsampling**"""

# downsamplig atoxic
atoxic_downsampled = resample(atoxic_list, replace=False, n_samples=total_toxic_seqs, random_state=RANDOM_SEED)
print('Total protein sequences in atoxic list post-downsampling: ', len(atoxic_downsampled))

"""**Comparing pre and post downsmaple distrbutions for sequence length**"""

# # downsampled distribution of atoxic protein sequence lengths
# dist_graph(df_atoxic_downsampled, 18, 8, 'length',
#            'Distribution of atoxic protein sequences lengths post downsampling',
#            'Sequence Length', '', 'seagreen', 80)

# # boxplots
# x = 2
# y = 12
# vert = True
# # atoxic pre downsample
# df_atoxic[df_atoxic.columns[2:].to_list()].boxplot(vert=vert)
# plt.gcf().set_size_inches(x,y)
# plt.title('Pre-downsampling atoxic sequences distribution')
# plt.show()
# print(df_atoxic['length'].describe().T[['mean', 'std', 'max','min', '25%', '50%', '75%']].round(decimals=2))
# print('\n')

# # atoxic post downsample
# df_atoxic_downsampled[df_atoxic_downsampled.columns[2:].to_list()].boxplot(vert=vert)
# plt.gcf().set_size_inches(x,y)
# plt.title('Post-downsampling atoxic sequences distribution')
# plt.show()
# print(df_atoxic_downsampled['length'].describe().T[['mean', 'std', 'max','min', '25%', '50%', '75%']].round(decimals=2))
# print('\n')

# # toxic
# df_toxic[df_toxic.columns[2:].to_list()].boxplot(vert=vert)
# plt.gcf().set_size_inches(x,y)
# plt.title('Toxic sequences distribution')
# plt.show()
# print(df_toxic['length'].describe().T[['mean', 'std', 'max','min', '25%', '50%', '75%']].round(decimals=2))

"""# **5. Combining Datasets & Appending Values**"""

# appending toxic and atoxic 'ProteinSequence' lists
proteins = toxic_list + atoxic_downsampled
print(len(proteins))

# -----------------------------------------------
# FUNCTIONS APPENDING ATCHLEY VALUES

# get atchley values - returns list of values
def get_atchley_values_list(letters, idx):
  return [float(dict_atchley.get(i)[idx]) for i in letters]

# matches each atchley value to each amino acid
def get_atchley_values_raw(sequence, seq_dict_r):
  seq_dict_r['f1'] = get_atchley_values_list(sequence, 0)
  seq_dict_r['f2'] = get_atchley_values_list(sequence, 1)
  seq_dict_r['f3'] = get_atchley_values_list(sequence, 2)
  seq_dict_r['f4'] = get_atchley_values_list(sequence, 3)
  seq_dict_r['f5'] = get_atchley_values_list(sequence, 4)

# -----------------------------------------------
# calculates sequential change for single atchley value
def get_change_list(atchley_list):
  atchley_list.insert(0, 0)
  change_list = [i for i in (np.diff(atchley_list))]
  atchley_list.pop(0)
  return change_list

# calculates sequential change for each atchley value
def get_atchley_diff(seq_dict_r, seq_dict_d):
  seq_dict_d['f1_d'] = get_change_list(seq_dict_r.get('f1'))
  seq_dict_d['f2_d'] = get_change_list(seq_dict_r.get('f2'))
  seq_dict_d['f3_d'] = get_change_list(seq_dict_r.get('f3'))
  seq_dict_d['f4_d'] = get_change_list(seq_dict_r.get('f4'))
  seq_dict_d['f5_d'] = get_change_list(seq_dict_r.get('f5'))

# appends the atchley values to the dictionary of a ProteinSequence object
def append_atchley_values(proteins_list):
  for protein in proteins_list:
    get_atchley_values_raw(protein.sequence, protein.seq_dict_raw)
    get_atchley_diff(protein.seq_dict_raw, protein.seq_dict_diff)

# appending atchley values
append_atchley_values(proteins)

"""**Appending Matrices**"""

# -----------------------------------------------
# FUNCTIONS UPDATING MATRICES

# returns matrix from input dictionary
def get_matrix_values(seq_dict):
  return np.array([seq_dict[i] for i in seq_dict.keys()])

# updates both matrices
def update_matrices(protein_seq_list):
  for protein in protein_seq_list:
    protein.matrix_raw = get_matrix_values(protein.seq_dict_raw)
    protein.matrix_diff = get_matrix_values(protein.seq_dict_diff)

# updating matrices
update_matrices(proteins)

# pickling protein list objects to file for later use
pickle_method(f_train_proteins_list, 'wb', proteins)

# unpickling protein list objects
proteins = pickle_method(f_train_proteins_list, 'rb', '')

# adding proteins list to dataframes
df_raw = pd.DataFrame.from_records([p.to_dict_raw() for p in proteins])
df_diff = pd.DataFrame.from_records([p.to_dict_diff() for p in proteins])
df_combined = pd.DataFrame.from_records([p.to_dict_combined() for p in proteins])

# pickling df files
pickle_method(f_train_atchley_raw, 'wb', df_raw)
pickle_method(f_train_atchley_diff, 'wb', df_diff)
pickle_method(f_train_atchley_combined, 'wb', df_combined)

# unpickling df files
df_raw = pickle_method(f_train_atchley_raw, 'rb', '')
df_diff = pickle_method(f_train_atchley_diff, 'rb', '')
df_combined = pickle_method(f_train_atchley_combined, 'rb', '')

df_raw.head(5)

df_diff.head(5)

df_combined.head(5)

print((type(df_raw['atchley_raw_avg'][0])))
print((type(df_diff['atchley_diff_avg'][0])))
print((type(df_combined['atchley_diff_avg'][0])))
print((type(df_combined['matrix_raw'][0])))
print((type(df_combined['matrix_diff'][0])))

print('Raw training set info:\n')
df_raw.info()

print('\Diff training set info:\n')
df_diff.info()

print('\nCombined training set info:\n')
df_combined.info()

print('Checking value counts for each class in df raw:\n1 == toxic\n0 == atoxic\n----------')
print(df_raw['toxic'].value_counts())

print('\nChecking value counts for each class in df diff:\n1 == toxic\n0 == atoxic\n----------')
print(df_diff['toxic'].value_counts())

print('\nChecking value counts for each class in df combined:\n1 == toxic\n0 == atoxic\n----------')
print(df_combined['toxic'].value_counts())

"""---


# **6. Training vs Validation splitting**
"""

# splitting variables
K_FOLDS = 5
VAL_SIZE = 0.20

# get splits from k fold
def get_kfold_splits(fname, cv_splits, val_size):
  fold_dict = {}
  i = 1
  sss = StratifiedShuffleSplit(n_splits=cv_splits, test_size=val_size, random_state=random_seed)
  print(sss)
  print('Number of k-fold splits: ', sss.get_n_splits())
  for train_idx, val_idx in sss.split(x_features, y_labels):
    x_train, x_val = x_features.iloc[train_idx], x_features.iloc[val_idx]
    y_train, y_val = y_labels.iloc[train_idx], y_labels.iloc[val_idx]
    fold_dict[str(i)] = [x_train, x_val, y_train, y_val]
    i += 1
  pickle_method(fname, 'wb', fold_dict)
  return fold_dict

# split into training and validation sets
y_labels = df_means['toxic']
x_features = df_means.drop(['toxic', 'identifier', 'letter'], axis=1)

# checking length
print(len(y_labels))
print(len(x_features))

# getting k-fold splits
# indexes:
# 0 -> x_train_vectorised
# 1 -> x_val_vectorised
# 2 -> y_train
# 3 -> y_val

# k_folds_dict = get_kfold_splits(f_5_fold_p20, K_FOLDS, VAL_SIZE)

# loading k-fold files

fold_5_10p = pickle_method(f_5_fold_p10, 'rb', '')
print('k-folds: ', len(fold_5_10p))
fold_5_15p = pickle_method(f_5_fold_p15, 'rb', '')
print('k-folds: ', len(fold_5_15p))
fold_5_20p = pickle_method(f_5_fold_p20, 'rb', '')
print('k-folds: ', len(fold_5_20p))

fold_10_10p = pickle_method(f_10_fold_p10, 'rb', '')
print('k-folds: ', len(fold_10_10p))
fold_10_15p = pickle_method(f_10_fold_p15, 'rb', '')
print('k-folds: ', len(fold_10_15p))
fold_10_20p = pickle_method(f_10_fold_p20, 'rb', '')
print('k-folds: ', len(fold_10_20p))

fold_15_10p = pickle_method(f_15_fold_p10, 'rb', '')
print('k-folds: ', len(fold_15_10p))
fold_15_15p = pickle_method(f_15_fold_p15, 'rb', '')
print('k-folds: ', len(fold_15_15p))
fold_15_20p = pickle_method(f_15_fold_p20, 'rb', '')
print('k-folds: ', len(fold_15_20p))

fold_20_10p = pickle_method(f_20_fold_p10, 'rb', '')
print('k-folds: ', len(fold_20_10p))
fold_20_15p = pickle_method(f_20_fold_p15, 'rb', '')
print('k-folds: ', len(fold_20_15p))
fold_20_20p = pickle_method(f_20_fold_p20, 'rb', '')
print('k-folds: ', len(fold_20_20p))

# # print test
# print('Dictionary length: ', len(k_folds_dict))
# for fold in k_folds_dict:
#   print('\nFold: ', fold)
#   print('Total training features: ', len(k_folds_dict[fold][2]))
#   print('Total validation features: ', len(k_folds_dict[fold][3]))
#   print('Total training labels: ', len(k_folds_dict[fold][2]))
#   print('Total validation labels: ', len(k_folds_dict[fold][3]))

"""---


# **6. Modelling**
"""

# fit model to data
def fit_model_to_folds(model, k_folds_dict):
  scores = []
  i = 1
  print('\nModel: ', model)
  for fold in k_folds_dict:
    fitted_model = model.fit(k_folds_dict[fold][0], k_folds_dict[fold][2])
    score = accuracy_score(k_folds_dict[fold][3], fitted_model.predict(k_folds_dict[fold][1]))
    print('{} of KFold {}'.format(fold, len(k_folds_dict)), ' --> ROC AUC score:', score)
    scores.append(score)
    i += 1
  print("\nMean model score: %.3f" % statistics.mean(scores))
  return scores

def run_models(models, k_folds_dict):
  models_dict = {}
  for model in models:
    models_dict[str(model)] = fit_model_to_folds(model, k_folds_dict)
  return models_dict

# def run_k_folds =

from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression,SGDClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score

classifiers = [SGDClassifier(loss='hinge', random_state=random_seed)]

model_dict = run_models(classifiers, fold_5_10p)
model_dict = run_models(classifiers, fold_5_15p)
model_dict = run_models(classifiers, fold_5_20p)

model_dict = run_models(classifiers, fold_10_10p)
model_dict = run_models(classifiers, fold_10_15p)
model_dict = run_models(classifiers, fold_10_20p)

model_dict = run_models(classifiers, fold_15_10p)
model_dict = run_models(classifiers, fold_15_15p)
model_dict = run_models(classifiers, fold_15_20p)

model_dict = run_models(classifiers, fold_20_10p)
model_dict = run_models(classifiers, fold_20_15p)
model_dict = run_models(classifiers, fold_20_20p)

"""---"""

# # pickling atoxic protein data
# pickle_method(op_atoxic_atchley_final, 'wb', atoxic_list)

# # unpickling complete atoxic protein data
# atoxic_list = pickle_method(op_atoxic_atchley_complete, 'rb', '')