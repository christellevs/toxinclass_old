


# -----------------------------------------------------------------------------

class Atchley:
    
    def __init__(self):
        
        
    def to_dict(self):
        return 
        

    def _get_atchley_values_list(self, aminos:List[str], idx:int) -> List[int]:
        """
        Returns all atchley values in a list for a specific amino acid.
        """
        
        return [float(cfg.DICT_ATCHLEY.get(i)[idx]) for i in aminos]


    def _get_atchley_values_raw(self, sequence, seq_dict_r):
        """
        Returns raw atchley values in a dictionary.
        """
        seq_dict_r['f1'] = get_atchley_values_list(sequence, 0)
        seq_dict_r['f2'] = get_atchley_values_list(sequence, 1)
        seq_dict_r['f3'] = get_atchley_values_list(sequence, 2)
        seq_dict_r['f4'] = get_atchley_values_list(sequence, 3)
        seq_dict_r['f5'] = get_atchley_values_list(sequence, 4)

    # -----------------------------------------------

    def _get_change_list(self, atchley_list):
        """
        Calculates sequential change for single atchley value.
        """
        atchley_list.insert(0, 0)
        change_list = [i for i in (np.diff(atchley_list))]
        atchley_list.pop(0)
        return change_list


    def get_atchley_diff(self, seq_dict_r, seq_dict_d):
        """
        Returns the sequential change for each atchley value as a dictionary.
        """
        seq_dict_d['f1_d'] = get_change_list(seq_dict_r.get('f1'))
        seq_dict_d['f2_d'] = get_change_list(seq_dict_r.get('f2'))
        seq_dict_d['f3_d'] = get_change_list(seq_dict_r.get('f3'))
        seq_dict_d['f4_d'] = get_change_list(seq_dict_r.get('f4'))
        seq_dict_d['f5_d'] = get_change_list(seq_dict_r.get('f5'))
    
    
    def append_atchley_values(self, proteins_list):
        """
        Appends the atchley values to the dictionary of a ProteinSequence object.
        """
        for protein in proteins_list:
            get_atchley_values_raw(protein.sequence, protein.seq_dict_raw)
            get_atchley_diff(protein.seq_dict_raw, protein.seq_dict_diff)