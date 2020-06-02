# module imports
# import config as cfg
import protein as prot
import helper_functions as hp

# other imports
import Bio as bio
import csv
import io
import pickle
import os.path
import random
import re
import statistics

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

# seed setting
random.seed(cfg.RANDOM_SEED)

# -----------------------------------------------------------------------------
# START

# opens fasta file and creates dataframe
def parse_fasta(path_fasta, is_toxic):
  sequences = []
  with open(path_fasta) as fasta_file:
    for title, sequence in SimpleFastaParser(fasta_file):
      if not any(ele in sequence for ele in cfg.INVALID_AMINO):
        sequences.append( prot.Protein(title.split(None, 1)[0],
                                        is_toxic, len(sequence), hp.split_seq(sequence)) )
  return sequences




# END 
# -----------------------------------------------------------------------------