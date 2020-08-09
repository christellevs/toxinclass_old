import pprint as p
import pytest

# import sys
# sys.path.insert(0, 'modules')

# local imports
from modules.atchley import Atchley
from modules.config import config as cfg

# -----------------------------------------------------------------------------

atchley = Atchley(f_in_atchley=cfg.f_in_atchley, f_out_atchley=cfg.f_out_atchley)




def test_atchley_to_df():
    """
    """
    df = atchley._atchley_to_df(f_in_atchley=atchley.f_in_atchley)
    p.pprint(df)
    
# -----------------------------------------------------------------------------
