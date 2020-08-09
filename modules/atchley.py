
import pandas as pd
import pprint as p

from typing import List, Dict

# local imports
import config as cfg
import utils as u

# -----------------------------------------------------------------------------

class Atchley:    
    def __init__(self, f_in_atchley:str, f_out_atchley:str):
        self.f_in_atchley = f_in_atchley
        self.f_out_atchley = f_out_atchley
        
    def get_atchley_dict(self) -> Dict[str, List[float]]:
        """
        Returns a parsed Dictionary of the 5 Atchley values per amino acid in a list.
        """
        return u.pickle_method(self.f_out_atchley, 'rb', '')

    def parse_atchley(self) -> None:
        """
        Parses a .txt file containing Atchley values into a dataframe for easy manipulation.
        Pickles the file for later use.
        """
        df = self._atchley_to_df(f_in_atchley=self.f_in_atchley)
        df = self._fix_text(df=df)
        self._atchley_to_dict(df=df, f_atchley_dict=self.f_atchley_dict)
        
    # -------------------------------------------------

    def _atchley_to_dict(self, df:pd.DataFrame, f_out_atchley:str) -> Dict:
        """
        Takes in a DataFrame of Atchley values and returns it as a Dict.
        """
        dict_atchley = df.T.to_dict('list')
        u.pickle_method(f_out_atchley, 'wb', dict_atchley)
        return dict_atchley
    
    def _atchley_to_df(self, f_in_atchley:str) -> pd.DataFrame:
        """
        Takes in .txt file containing Atchley values.
        Returns values in a pandas DataFrame.
        """
        df = u.cvs_to_df(filename=f_in_atchley, col_idx=0)
        df.rename(columns={'amino.acid': 'amino_acid'}, inplace=True)
        df.set_index('amino_acid', inplace=True)
        return df
    
    def _fix_text(self, df:pd.DataFrame) -> pd.DataFrame:
        """
        Ensures dash in the original Atchley test is a negative sign.
        """
        for col in df['f1': 'f5']:
            df[col] = df[col].apply(lambda x: re.sub(r'[^\x00-\x7F]+','-', x)).astype(float)
        return df

# -----------------------------------------------------------------------------

