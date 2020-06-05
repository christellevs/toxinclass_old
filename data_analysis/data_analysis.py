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

import sys
sys.path.insert(1, 'modules')

# module imports
import config as cfg
import protein as prot
import helper_functions as hp
import fasta_parser

"""# **1. Atchley Value Analysis**"""
# # unpickling protein list objects
# proteins = pickle_method(f_train_proteins_list, 'rb', '')
# # unpickling df files
# df_raw = pickle_method(f_train_atchley_raw, 'rb', '')
# df_diff = pickle_method(f_train_atchley_diff, 'rb', '')
# df_combined = pickle_method(f_train_atchley_combined, 'rb', '')


df_combined = hp.pickle_method(cfg.f_train_atchley_combined, 'rb', '')
df_combined.head(5)

df_toxic = df_combined.loc[df_combined['toxic'] == 1]
df_toxic.head(5)

df_atoxic = df_combined.loc[df_combined['toxic'] == 0]
df_atoxic.head(5)


# get avgs of columns with arrays, returns arrays same size
def get_average(df, column):
    return np.average(df[column], axis=0)


def get_atchley_avgs(df, labels):
    return dict((el, get_average(df, el)) for el in labels)


raw_labels = ['f1_raw', 'f2_raw', 'f3_raw', 'f4_raw', 'f5_raw', 'atchley_raw_avg']
diff_labels = ['f1_diff', 'f2_diff', 'f3_diff', 'f4_diff', 'f5_diff', 'atchley_diff_avg']

df_toxic_raw = pd.DataFrame.from_dict(get_atchley_avgs(df_toxic, raw_labels))
df_toxic_diff = pd.DataFrame.from_dict(get_atchley_avgs(df_toxic, diff_labels))
df_atoxic_raw = pd.DataFrame.from_dict(get_atchley_avgs(df_atoxic, raw_labels))
df_atoxic_diff = pd.DataFrame.from_dict(get_atchley_avgs(df_atoxic, diff_labels))

# print test
print('Toxic raw\n', df_toxic_raw.head(5))
print('\nToxic change\n', df_toxic_diff.head(5))
print('\nAtoxic raw\n', df_atoxic_raw.head(5))
print('\nAtoxic change\n', df_atoxic_diff.head(5))

# DF LIST
df_list = [df_toxic_raw, df_toxic_diff, df_atoxic_raw, df_atoxic_diff]

# BOXPLOTS
BOXPLOT_XDIM = 20
BOXPLOT_YDIM = 10

titles_boxplot = ['Distribution of Atchley Values of a Toxic Protein Sequence',
                  'Distribution of Change in Atchley Values of a Toxic Protein Sequence',
                  'Distribution of Atchley Values of an Atoxic Protein Sequence',
                  'Distribution of Change in Atchley Values of an Atoxic Protein Sequence']


# plots boxplots
def plot_boxplot(df, title, idx1, idx2):
    print('\n')
    sns.set(style='whitegrid')
    fig = plt.subplots(figsize=(BOXPLOT_XDIM, BOXPLOT_YDIM))
    ax = sns.boxplot(data=df[df.columns[idx1:idx2].tolist()], palette='Set2', orient='h')
    ax.set_title(title, fontsize=26)
    ax.set_xlabel('Atchley Feature', fontsize=16)
    ax.set_ylabel('Atchley Value', fontsize=16)
    # print(df['length'].describe().T[['mean', 'std', 'max','min', '25%', '50%', '75%']].round(decimals=2))


# plots all boxplots
def plot_all_boxplots(df_list, titles, idx1, idx2):
    plot_boxplot(df_list[0], titles[0], idx1, idx2)
    plot_boxplot(df_list[1], titles[1], idx1, idx2)
    plot_boxplot(df_list[2], titles[2], idx1, idx2)
    plot_boxplot(df_list[3], titles[3], idx1, idx2)


# plotting all boxplots
plot_all_boxplots(df_list, titles_boxplot, 0, (-1))

"""**Plotting Atchley Values**"""

# PLOTTING ATCHLEY VALUES FUNCTIONS
# -----------------------------------------------------------------------------

# ATCHLEY PLOT VARIABLES
plot_raw_labels = ['f1-Polarity', 'f2-Secondary Structure', 'f3-Molecular Volume', 'f4-Relative Composition',
                   'f5-Electrostatic Charge', 'Average']
plot_diff_labels = ['f1-Change in Polarity', 'f2-Change in Secondary Structure', 'f3-Change in Molecular Volume',
                    'f4-Change in Relative Composition', 'f5-Change Electrostatic Charge', 'Average Change']
titles_linplot = ['Atchley Values of a Toxic Protein Sequence', 'Change in Atchley Values of a Toxic Protein Sequence',
                  'Atchley Values of an Atoxic Protein Sequence',
                  'Change in Atchley Values of an Atoxic Protein Sequence']


def plot_atchley_values(df, title, labels, xdim, ydim):
    sns.set()
    fig = plt.subplots(figsize=(xdim, ydim))
    sns.set(style='whitegrid')
    plt.ylim(top=3)
    plt.ylim(bottom=-3)
    plt.plot(df)
    plt.legend(labels, loc='upper right', title='Atchley Value')
    plt.ylabel('Atchley Value')
    plt.xlabel('Amino Acid Position in Sequence')
    plt.title(title, fontsize=22)
    plt.show()
    print('\n')


def plot_all_atchley_values(titles, labels_raw, labels_diff, xdim, ydim):
    plot_atchley_values(df_toxic_raw, titles[0], labels_raw, xdim, ydim)
    plot_atchley_values(df_toxic_diff, titles[1], labels_raw, xdim, ydim)
    plot_atchley_values(df_atoxic_raw, titles[2], labels_diff, xdim, ydim)
    plot_atchley_values(df_atoxic_diff, titles[3], labels_diff, xdim, ydim)


# plotting atchley values
plot_all_atchley_values(titles_linplot, plot_raw_labels, plot_diff_labels, 20, 10)

