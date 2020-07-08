import os
import re
import pprint as p

from sklearn.utils import resample
from typing import List

import sys
sys.path.insert(1, 'modules')

# module imports
import config as cfg
import helper_funcs as hf
import protein as prot

# -----------------------------------------------------------------------------
# ILP PRE-PROCESSING FUNCTIONS

def parse_identifier(identifier):
  pattern = '\|(.*?)\|'
  return re.search(pattern, identifier).group(1)

def resample_proteins(proteins):
  return resample(proteins, replace=False, n_samples=cfg.MAX_PROT_NUM, random_state=cfg.RANDOM_SEED)

def get_positions_list(fname:str, proteins):
  for p in proteins:
    seq = p.sequence
    fname.write('\nposition(p_%s,%s,%s).' % (parse_identifier(p.identifier),
                                         [(seq.index(a) + 1) for a in seq],
                                         [a.lower() for a in seq]) )


def get_positions_lines(fname:str, proteins:List[prot.Protein]) -> None:
  for p in proteins:
    for idx, amino in enumerate(p.sequence):
      fname.write('\nposition(p_%s,%s,%s).' %
                  (parse_identifier(p.identifier), (idx + 1), amino.lower() ))


def separate_proteins(proteins:List[prot.Protein]) -> List[prot.Protein]:
  """Separate proteins into different lists depending on toxicity."""
  toxic_proteins = [prot for prot in proteins if prot.toxic == 1]
  atoxic_proteins = [prot for prot in proteins if prot.toxic == 0]
  return toxic_proteins, atoxic_proteins
  

def get_examples(fname:str, proteins:List[prot.Protein]) -> None:
  for p in proteins:
    if p.toxic == 1:
      fname.write('\nexample(toxic(p_%s),1).' % (parse_identifier(p.identifier)))
    else:
      fname.write('\nexample(toxic(p_%s),-1).' % (parse_identifier(p.identifier)))
      
      
# -----------------------------------------------------------------------------



# unloading complete df
train_proteins = hf.pickle_method(cfg.f_train_proteins, 'rb', '')
print(f'Number of proteins {len(train_proteins)}')

# cropping protein sequences
train_proteins_crop = hf.crop_sequences(train_proteins)

toxic_proteins, atoxic_proteins = separate_proteins(train_proteins_crop)

print(f'Number of toxic proteins {len(toxic_proteins)}')
print(f'Number of atoxic proteins {len(atoxic_proteins)}')


# downsampling protein sequences
toxic_prot_down = resample_proteins(toxic_proteins)
atoxic_prot_down = resample_proteins(atoxic_proteins)

print(f'Number of downsampled toxic proteins {len(toxic_prot_down)}')
print(f'Number of downsampled atoxic proteins {len(atoxic_prot_down)}')

downsampled_proteins = toxic_prot_down + atoxic_prot_down

# -----------------------------------------------------------------------------
# WRITING TO FILES

# positions
hf.write_to_file(cfg.f_ILP_positions, downsampled_proteins, get_positions_lines)

# examples
hf.write_to_file(cfg.f_ILP_examples, downsampled_proteins, get_examples)
