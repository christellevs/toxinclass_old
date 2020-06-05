import io
import itertools as ite
import numpy as np
import pandas as pd
import pickle
import os.path
import random
import statistics

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.utils import resample

import sys
sys.path.insert(1, 'modules')

random.seed(cfg.RANDOM_SEED)


# Loading Pre-Processed datasets
df_combined = pickle_method(f_train_atchley_combined, 'rb', '')
df_combined.head(5)

# -----------------------------------------------------------------------------
# SPLITTING IMPORTS

# split into training and validation sets
drop_cols = ['toxic', 'identifier', 'sequence', 'atchley_raw_avg', 'atchley_diff_avg']
y_labels = df_combined['toxic']
x_features = df_combined['matrix_raw']

# checking length
print(len(y_labels))
print(len(x_features))

# SPLITTING FUNCTIONS START
# -----------------------------------------------------------------------------
sss = StratifiedShuffleSplit(n_splits=cfg.K_FOLDS, test_size=cfg.VAL_SIZE, random_state=cfg.RANDOM_SEED)

# get splits from k fold
def get_kfold_splits(fname):
  fold_dict = {}
  i = 1
  print(sss)
  print('Number of k-fold splits: ', sss.get_n_splits())
  for train_idx, val_idx in sss.split(x_features, y_labels):
    x_train, x_val = x_features.iloc[train_idx], x_features.iloc[val_idx]
    y_train, y_val = y_labels.iloc[train_idx], y_labels.iloc[val_idx]
    fold_dict[str(i)] = [x_train, x_val, y_train, y_val]
    i += 1
  pickle_method(fname, 'wb', fold_dict)

# get splits from k fold
def get_kfold_dict():
  fold_dict = {}
  i = 1
  print(sss)
  print('Number of k-fold splits: ', sss.get_n_splits())
  for train_idx, val_idx in sss.split(x_features, y_labels):
    x_train, x_val = x_features.iloc[train_idx], x_features.iloc[val_idx]
    y_train, y_val = y_labels.iloc[train_idx], y_labels.iloc[val_idx]
    fold_dict[str(i)] = [x_train, x_val, y_train, y_val]
    i += 1
  return fold_dict

"""**Stratified Splitting Dataset**"""

# getting k-fold splits
# indexes:
# 0 -> x_train_vectorised
# 1 -> x_val_vectorised
# 2 -> y_train
# 3 -> y_val
get_kfold_splits(f_20_fold_p20)

"""**Loading back K-Folds Splits**"""

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

"""**Single Dictionary**"""

# loading k-fold files
k_fold_dicts = {'fold_5_10p': fold_5_10p,
                'fold_5_15p': fold_5_15p,
                'fold_5_20p': fold_5_20p,
                'fold_10_10p': fold_10_10p,
                'fold_10_15p': fold_10_15p,
                'fold_10_20p': fold_10_20p,
                'fold_15_10p': fold_15_10p,
                'fold_15_15p': fold_15_15p,
                'fold_15_20p': fold_15_20p,
                'fold_20_10p': fold_20_10p,
                'fold_20_15p': fold_20_15p,
                'fold_20_20p': fold_20_20p}

print('Total k-folds: ', len(k_fold_dicts))

"""# **3. Modelling**"""

# -----------------------------------------------------------------------------
# COMMON IMPORTS

import multiprocessing

from IPython.display import Image

from keras.callbacks import EarlyStopping
from keras.initializers import Constant
from keras.layers import Conv1D, Conv2D, Embedding, Flatten, GlobalMaxPooling1D
from keras.layers import GRU, MaxPooling1D, SpatialDropout1D
from keras.layers.core import Activation, Dropout, Dense
from keras.layers.embeddings import Embedding
from keras.models import Sequential
from keras.utils.vis_utils import model_to_dot

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.tree import DecisionTreeClassifier

# MODELLING FUNCTIONS START
# -----------------------------------------------------------------------------

# fit model to data
def fit_model_to_folds(model, k_folds_dict):
  scores = []
  i = 1
  for fold in k_folds_dict:
    fitted_model = model.fit(k_folds_dict[fold][0], k_folds_dict[fold][2])
    score = accuracy_score(k_folds_dict[fold][3], fitted_model.predict(k_folds_dict[fold][1]))
    print('{} of k-fold {}'.format(fold, len(k_folds_dict)), ' --> ROC AUC score:', score)
    scores.append(score)
    i += 1
  mean_score = statistics.mean(scores)
  print("\nMean model score: %.3f" % mean_score)
  return mean_score

# run model, save results to dictionary
def run_model(model, k_folds_dicts, filename):
  model_results = {}
  for key, value in k_folds_dicts.items():
    print('\nCV:', key)
    model_results[key] = fit_model_to_folds(model, value)
  pickle_method(filename, 'wb', model_results)
  return model_results

# run all models, save results to a dictionary of dictionaries
def run_all_models(models, k_fold_dicts):
  models_dict = {}
  for model in models:
    print('\nModel:', model)
    models_dict[model] = run_model(model)
  return models_dict

"""# **3.1 Complex Modelling**

**CNN Model**
"""

# Common state model variables
EPOCHS = 10
BATCH_SIZE = 128

# testing on single val fold
trial_5_20p_fold = get_kfold_dict()

print(len(trial_5_20p_fold.get('1')[0]))

# creating pre-trained embedding layer for LSTM and CNN
embedding_layer = Embedding(vocab_size, EMBEDDING_DIM, input_length=MAX_WORDS, weights=[embedding_matrix], trainable=True)

# enabling callbacks so model stops in case validation loss keeps increasing
callbacks = [EarlyStopping(monitor='val_loss', min_delta=0, patience=2, verbose=0, mode='auto')]

# defining CNN model
m_cnn = Sequential()
m_cnn.add(Conv1D(filters=32, kernel_size=8, activation='relu'))
m_cnn.add(MaxPooling1D(pool_size=2))
m_cnn.add(Flatten())
m_cnn.add(Dense(256, activation='relu'))
m_cnn.add(Dense(1, activation='sigmoid'))
m_cnn.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

m_cnn_hist = fit_model_to_folds(m_cnn, test_5_20p_fold)

# CNN model summary
m_cnn.summary()
Image(model_to_dot(m_cnn).create(prog='dot', format='png'))

"""# **3.2 Other Modelling**"""

# model linear sgd classifier: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.SGDClassifier.html
m_linear_sgd = SGDClassifier(loss='modified_huber', random_state=RANDOM_SEED)
m_linear_sgd_dict = run_model(m_linear_sgd, k_fold_dicts, f_m_linear_sgd)

# model random forest: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
m_rfc = RandomForestClassifier(random_state=RANDOM_SEED)
m_rfc_dict = run_model(m_rfc, k_fold_dicts, f_m_random_foregit sst)

# model decision tree: https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeClassifier.html
m_dt = DecisionTreeClassifier(random_state=RANDOM_SEED)
m_dt_dict = run_model(m_dt, k_fold_dicts, f_m_decision_tree)

"""# **Running all models**"""

# models list
models = [m_linear_sgd, m_rfc]

# running all models
models_dict = run_models(models, k_fold_dicts)