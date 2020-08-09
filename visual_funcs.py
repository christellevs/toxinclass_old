import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sklearn
import statistics

from textwrap import wrap

# -----------------------------------------------------------------------------
# GRAPH FUNCTIONS

def dist_graph(df, x_dim, y_dim, x_axis_col, title, x_label, y_label, colour, bins):
  """
  Creates a distribution graph.
  """
  sns.set(style='whitegrid')
  fig = plt.subplots(figsize=(x_dim, y_dim))
  dist = np.asarray(df[[x_axis_col]].values.tolist()).ravel()
  ax = sns.distplot(df[x_axis_col],
                       bins=bins, kde=True,
                       kde_kws={"color": "k", "lw": 2, "label": "KDE"},
                       color=colour)
  ax.set_title(title, fontsize=26)
  ax.set_xlabel(x_label, fontsize=16)
  ax.set_ylabel(y_label, fontsize=16)
  

def describe_df(df, decimal):
  """
  Describes a DataFrame table.
  """
  return df.describe().T[['mean', 'std', 'max','min', '25%','50%', '75%']].round(decimals=decimal)
  

# a = [1, 3, 5, 44, 5]
# b = ['a', 'd', 't', 'o', 'p']

# for i, (vol1, vol2) in enumerate(zip(a, b):
#   print(i,vol1, vol2)


# df = df.drop('letters', 1).assign(**df['letters'].dropna().apply(pd.Series))

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