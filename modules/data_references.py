""" Data References for using in pre-processing
----------------------------------------------------------------------------- """

import sys
sys.path.insert(1, 'modules')

# module imports
import config as cfg
import helper_funcs as hf

# loading pickled Atchley dictionary
dict_atchley = hf.pickle_method(cfg.f_atchley_dict, 'rb', '')