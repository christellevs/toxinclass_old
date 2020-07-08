# File to define global variables

# links
# TOXIFY: https://github.com/tijeco/toxify
# ToxClassifier: https://github.com/rgacesa/ToxClassifier

import random
import helper_funcs as hf

# -----------------------------------------------------------------------------

RANDOM_SEED = 273

# pre-processing
INVALID_AMINO = ['B', 'J', 'O', 'U', 'X', 'Z']
MAX_SEQ_LEN = 500
MAX_PROT_NUM = 300

# splitting
K_FOLDS = 10
VAL_SIZE = 0.1


# FILEPATHS
# -----------------------------------------------------------------------------

# common paths
PATH_DATASETS = 'data/datasets'
PATH_ILP_PROCESS = 'ILP/ILP_content'

# data files
PATH_DATA_TRAIN = PATH_DATASETS + '/training_data'
PATH_DATA_BENCH = PATH_DATASETS + '/benchmark_data'
PATH_DATA_TEST =  PATH_DATASETS + '/testing_data'

PATH_DATA_REFS = 'data/data_refs'
PATH_DATA_SPLITS = 'data/data_splits'


# MAIN DATA
# -------------------------------------

# training data
f_train_toxic_fasta = PATH_DATA_TRAIN + '/pre.venom.fasta'
f_train_atoxic_fasta = PATH_DATA_TRAIN + '/pre.NOT.venom.fasta'

# benchmark data
f_test_toxic_fasta = PATH_DATA_TRAIN + '/post.venom.fasta'
f_test_atoxic_fasta = PATH_DATA_TRAIN + '/post.NOT.venom.fasta'

#testing data

# processed main data
f_train_proteins = PATH_DATA_TRAIN + '/train_proteins.pickle'
f_train_df = PATH_DATA_TRAIN + '/train_df.pickle'


# REFERENCE DATA
# -------------------------------------

# files
f_aminoacid = PATH_DATA_REFS + '/aminoacids.csv'
f_atchley = PATH_DATA_REFS + '/atchley.txt'

f_atchley_dict = PATH_DATA_REFS + '/atchley_dict.pickle'
f_atchley_df = PATH_DATA_REFS + '/atchley_df.pickle'

# references
DICT_ATCHLEY = hf.pickle_method(f_atchley_dict, 'rb', '')


# ILP FILES
# -------------------------------------
f_ILP_classes = PATH_ILP_PROCESS + '/classes.pl'
f_ILP_nclasses = PATH_ILP_PROCESS + '/nclasses.pl'
f_ILP_positions = PATH_ILP_PROCESS + '/positions.pl'
f_ILP_examples = PATH_ILP_PROCESS + '/examples.pl'


# DATASET SPLITTING
# -------------------------------------
# filename format: kN_pM
# kn == no. of splits
# pM == % taken for split

# # k-fold cross validation dictionaries

# f_5_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p10.pickle'
# f_5_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p15.pickle'
# f_5_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p20.pickle'

f_k10_p10 = PATH_DATA_SPLITS + '/k10_p10.pickle'
# f_10_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/10_fold_p15.pickle'
# f_10_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/10_fold_p20.pickle'

# f_15_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p10.pickle'
# f_15_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p15.pickle'
# f_15_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p20.pickle'

# f_20_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p10.pickle'
# f_20_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p15.pickle'
# f_20_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p20.pickle'

# MODELS
# -------------------------------------
# f_m_linear_sgd = PATH_TO_FOLDER + '/models/m_linear_sgd.pickle'
# f_m_random_forest = PATH_TO_FOLDER + '/models/m_random_forest.pickle'
# f_m_decision_tree = PATH_TO_FOLDER + '/models/m_decision_tree.pickle'