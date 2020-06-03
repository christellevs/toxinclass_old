import seaborn as sns
import sklearn
import statistics

from textwrap import wrap

# -----------------------------------------------------------------------------
# GRAPH FUNCTIONS

# make distribution plot
def dist_graph(df, x_dim, y_dim, x_axis_col, title, x_label, y_label, colour, bins):
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

# a = [1, 3, 5, 44, 5]
# b = ['a', 'd', 't', 'o', 'p']

# for i, (vol1, vol2) in enumerate(zip(a, b):
#   print(i,vol1, vol2)


# df = df.drop('letters', 1).assign(**df['letters'].dropna().apply(pd.Series))