
import numpy as np

# PROTEIN PROCESSING FUNCTIONS START
# -----------------------------------------------------------------------------

# PROTEIN CLASS
class Protein:
  def __init__(self, identifier, toxic, length, sequence):
    self.identifier = identifier
    self.toxic = toxic
    self.length = length
    self.sequence = sequence
    self.seq_dict_raw = {}
    self.seq_dict_diff = {}
    self.matrix_raw = np.zeros((5, length))
    self.matrix_diff = np.zeros((5, length))


  

  def _to_dict(self):
    """
    Used to easily transform parsed protein data in Dictionary first, and then DataFrame.
    """
    return {'identifier': self.identifier,
            'toxic': self.toxic,
            'length': self.length,
            'sequence': self.sequence,
            'matrix_raw': self.matrix_raw,
            'matrix_diff': self.seq_dict_diff,
            'f1_raw': self.matrix_raw[0],
            'f2_raw': self.matrix_raw[1],
            'f3_raw': self.matrix_raw[2],
            'f4_raw': self.matrix_raw[3],
            'f5_raw': self.matrix_raw[4],
            'atchley_raw_avg': np.average(self.matrix_raw, axis=0),
            'f1_diff': self.matrix_diff[0],
            'f2_diff': self.matrix_diff[1],
            'f3_diff': self.matrix_diff[2],
            'f4_diff': self.matrix_diff[3],
            'f5_diff': self.matrix_diff[4],
            'atchley_diff_avg': np.average(self.matrix_diff, axis=0)}
    
    

  

# -----------------------------------------------------------------------------