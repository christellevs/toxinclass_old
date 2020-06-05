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

""" PRE-PROCESSING TRAINING SET
----------------------------------------------------------------------------- """

# PARSING FASTA
toxic_proteins = hf.parse_fasta(cfg.f_train_toxic_fasta, 1)
atoxic_proteins = hf.parse_fasta(cfg.f_train_atoxic_fasta, 0)

# print test
print('Total toxic sequences: ', len(toxic_proteins))
print('Total atoxic sequences: ', len(atoxic_proteins))


"""  DOWNSAMPLING
----------------------------------------------------------------------------- """

toxic_proteins = hf.crop_sequences(toxic_proteins)
atoxic_proteins = hf.crop_sequences(atoxic_proteins)

total_toxic_seqs = len(toxic_proteins)
print('Total toxic sequences ( <=', cfg.MAX_SEQ_LEN, '):', total_toxic_seqs)
print('Total toxic sequences ( <=', cfg.MAX_SEQ_LEN, '):', len(atoxic_proteins))


atoxic_proteins = resample(atoxic_proteins, replace=False, n_samples=total_toxic_seqs, random_state=cfg.RANDOM_SEED)
print('Total protein sequences in atoxic list post-downsampling: ', len(atoxic_proteins))


""" COMBINING & APPENDING
----------------------------------------------------------------------------- """

proteins = toxic_proteins + atoxic_proteins
hf.append_proteins(proteins)

print('Total overall protein sequences', len(proteins))
print(proteins[0].matrix_diff)


df_proteins = hf.proteins_to_df(proteins)

# print test
print('\Protein training set info:')
df_proteins.info()

print('\nChecking value counts for each class in df combined:\n1 == toxic\n0 == atoxic\n----------')
print(df_proteins['toxic'].value_counts())

""" SPLITTING DATA
----------------------------------------------------------------------------- """




# # split into training and validation sets
# y_labels = df_means['toxic']
# x_features = df_means.drop(['toxic', 'identifier', 'letter'], axis=1)

# # checking length
# print(len(y_labels))
# print(len(x_features))

# # getting k-fold splits
# # indexes:
# # 0 -> x_train_vectorised
# # 1 -> x_val_vectorised
# # 2 -> y_train
# # 3 -> y_val

# # k_folds_dict = get_kfold_splits(f_5_fold_p20, K_FOLDS, VAL_SIZE)

# # loading k-fold files

# fold_5_10p = pickle_method(f_5_fold_p10, 'rb', '')
# print('k-folds: ', len(fold_5_10p))
# fold_5_15p = pickle_method(f_5_fold_p15, 'rb', '')
# print('k-folds: ', len(fold_5_15p))
# fold_5_20p = pickle_method(f_5_fold_p20, 'rb', '')
# print('k-folds: ', len(fold_5_20p))

# fold_10_10p = pickle_method(f_10_fold_p10, 'rb', '')
# print('k-folds: ', len(fold_10_10p))
# fold_10_15p = pickle_method(f_10_fold_p15, 'rb', '')
# print('k-folds: ', len(fold_10_15p))
# fold_10_20p = pickle_method(f_10_fold_p20, 'rb', '')
# print('k-folds: ', len(fold_10_20p))

# fold_15_10p = pickle_method(f_15_fold_p10, 'rb', '')
# print('k-folds: ', len(fold_15_10p))
# fold_15_15p = pickle_method(f_15_fold_p15, 'rb', '')
# print('k-folds: ', len(fold_15_15p))
# fold_15_20p = pickle_method(f_15_fold_p20, 'rb', '')
# print('k-folds: ', len(fold_15_20p))

# fold_20_10p = pickle_method(f_20_fold_p10, 'rb', '')
# print('k-folds: ', len(fold_20_10p))
# fold_20_15p = pickle_method(f_20_fold_p15, 'rb', '')
# print('k-folds: ', len(fold_20_15p))
# fold_20_20p = pickle_method(f_20_fold_p20, 'rb', '')
# print('k-folds: ', len(fold_20_20p))

# # # print test
# # print('Dictionary length: ', len(k_folds_dict))
# # for fold in k_folds_dict:
# #   print('\nFold: ', fold)
# #   print('Total training features: ', len(k_folds_dict[fold][2]))
# #   print('Total validation features: ', len(k_folds_dict[fold][3]))
# #   print('Total training labels: ', len(k_folds_dict[fold][2]))
# #   print('Total validation labels: ', len(k_folds_dict[fold][3]))

# """---


# # **6. Modelling**
# """

# # fit model to data
# def fit_model_to_folds(model, k_folds_dict):
#   scores = []
#   i = 1
#   print('\nModel: ', model)
#   for fold in k_folds_dict:
#     fitted_model = model.fit(k_folds_dict[fold][0], k_folds_dict[fold][2])
#     score = accuracy_score(k_folds_dict[fold][3], fitted_model.predict(k_folds_dict[fold][1]))
#     print('{} of KFold {}'.format(fold, len(k_folds_dict)), ' --> ROC AUC score:', score)
#     scores.append(score)
#     i += 1
#   print("\nMean model score: %.3f" % statistics.mean(scores))
#   return scores

# def run_models(models, k_folds_dict):
#   models_dict = {}
#   for model in models:
#     models_dict[str(model)] = fit_model_to_folds(model, k_folds_dict)
#   return models_dict

# # def run_k_folds =

# from sklearn.naive_bayes import MultinomialNB
# from sklearn.linear_model import LogisticRegression,SGDClassifier
# from sklearn.svm import SVC
# from sklearn.metrics import classification_report,confusion_matrix,accuracy_score

# classifiers = [SGDClassifier(loss='hinge', random_state=random_seed)]

# model_dict = run_models(classifiers, fold_5_10p)
# model_dict = run_models(classifiers, fold_5_15p)
# model_dict = run_models(classifiers, fold_5_20p)

# model_dict = run_models(classifiers, fold_10_10p)
# model_dict = run_models(classifiers, fold_10_15p)
# model_dict = run_models(classifiers, fold_10_20p)

# model_dict = run_models(classifiers, fold_15_10p)
# model_dict = run_models(classifiers, fold_15_15p)
# model_dict = run_models(classifiers, fold_15_20p)

# model_dict = run_models(classifiers, fold_20_10p)
# model_dict = run_models(classifiers, fold_20_15p)
# model_dict = run_models(classifiers, fold_20_20p)

# """---"""

# # # pickling atoxic protein data
# # pickle_method(op_atoxic_atchley_final, 'wb', atoxic_list)

# # # unpickling complete atoxic protein data
# # atoxic_list = pickle_method(op_atoxic_atchley_complete, 'rb', '')