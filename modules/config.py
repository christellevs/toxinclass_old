# File to define global variables

RANDOM_SEED = 273

# common variables
INVALID_AMINO = ['B', 'J', 'O', 'U', 'X', 'Z']

# FILEPATHS
# -----------------------------------------------------------------------------

# common paths
PATH_DATA_TRAIN = 'data/datasets/training_data'
PATH_DATA_BENCH = 'data/datasets/benchmark_data'
PATH_DATA_TEST = 'data/datasets/testing_data'
PATH_DATA_REFS = 'data/data_refs'


# training data
f_train_toxic_fasta = PATH_DATA_TRAIN + '/pre.venom.fasta'
f__train_atoxic_fasta = PATH_DATA_TRAIN + '/pre.NOT.venom.fasta'

# benchmark data


# reference data
# reference values files
f_amino_acid = PATH_DATA_REFS + '/amino_acid_codes.csv'
f_atchley = PATH_DATA_REFS + '/atchley.csv'

f_atchley_dict = PATH_DATA_REFS + '/atchley_dict.pickle'
f_atchley_table = PATH_DATA_REFS + '/atchley_table.pickle'



# # processed data
# f_train_proteins_list = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_proteins_list.pickle'
# # f_train_complete = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_complete.pickle'
# # f_train_atchley_means = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_atchley_mean.pickle'

# f_train_atchley_raw = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_atchley_raw.pickle'
# f_train_atchley_diff = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_atchley_diff.pickle'
# f_train_atchley_combined = PATH_TO_FOLDER + '/pre_processing_files/dataframes/train_atchley_combined.pickle'

# # k-fold cross validation dictionaries
# f_5_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p10.pickle'
# f_5_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p15.pickle'
# f_5_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/5_fold_p20.pickle'

# f_10_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/10_fold_p10.pickle'
# f_10_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/10_fold_p15.pickle'
# f_10_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/10_fold_p20.pickle'

# f_15_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p10.pickle'
# f_15_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p15.pickle'
# f_15_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/15_fold_p20.pickle'

# f_20_fold_p10 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p10.pickle'
# f_20_fold_p15 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p15.pickle'
# f_20_fold_p20 = PATH_TO_FOLDER + '/pre_processing_files/data_splits/20_fold_p20.pickle'
