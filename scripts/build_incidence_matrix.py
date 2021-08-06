import pandas as pd
import numpy as np
import sys


def read_joiner(read_values):
    """A function to structure lists of items into a single pandas
    dataframe cell
    
    Parameters:
    -----------------------------
        :  read_values (pd.Series): column values, likely grouped
        
    Returns:
    -----------------------------
        : str_read_values (str): a string separated by a semi-colon with values
    """
    str_list = list(read_values.astype(str).to_list())
    return ";".join(str_list)


def spread(read_values):
    """A function to compute the number of base pairs between the most
    separated fragements in a read.
    
    Parameters:
    -----------------------------
        :  read_values (pd.Series): column values, likely grouped
        
    Returns:
    -----------------------------
        : max_contact_distance (float): number of base pairs between midpoints
        of the furthest contacts
    """
    return np.max(read_values) - np.min(read_values)


def build_incidence(df):
    """A function to construct the incidence matrix from the alignment 
    table
    
    
    Parameters:
    -----------------------------
        : df (pd.DataFrame): the alignment table
        
    Returns:
    -----------------------------
        : incidence (pd.DataFrame): the incidence table
    """
    incidence = df.groupby('read_name', as_index=False).agg(
        cardinality=('fragment_id', 'nunique'),
        max_contact_distance=('fragment_midpoint', spread),
        contact_midpoints=('fragment_midpoint', read_joiner),
        fragment_ids=('fragment_id', read_joiner),
        chromosome=('Chromosome', read_joiner),
    )

    incidence = incidence.sort_values(by='cardinality', ascending=False)
    return incidence

if __name__ == "__main__":
    input_path = sys.argv[1]
    df = pd.read_csv(input_path)
    df = build_incidence(df)
    df.to_csv(sys.stdout, index=False)