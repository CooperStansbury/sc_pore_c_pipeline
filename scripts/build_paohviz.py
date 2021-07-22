import pandas as pd
import numpy as np
import sys
import os


def build_paohvis_output(df):
    """A function to format the output for paohviz
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        
    Returns:
    -----------------------------
        : pao_df (pd.DataFrame): paohvis structured CSV
    """
    pao_df = df[['read_idx', 'fragment_id', 'chrom', 'read_name']].reset_index(drop=True)
    pao_df['slot_1'] = ""
    pao_df['slot_2'] = ""
    pao_df = pao_df[['read_idx', 'fragment_id', 'slot_1', 'slot_2', 'chrom', 'read_name']]
    return pao_df 


if __name__ == "__main__":
    input_path = sys.argv[1]
    df = pd.read_csv(input_path)
  
    df = build_paohvis_output(df)
    df.to_csv(sys.stdout, index=False)