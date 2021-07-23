import pandas as pd
import numpy as np
import sys


def load_coverage(in_path):
    """A function to read coverage information 
    to a pandas dataframe. Note: this funcation 
    expects the output of:
    
    `genomeCoverageBed -ibam <bam_file.bam> -bga`
    
    Parameters:
    -----------------------------
        : in_path (str): path to the coverage file
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): dataframe with column names
    """
    colnames = ['RefSeq accession', 
                'position_start', 
                'position_end', 
                'depth']
    
    df = pd.read_csv(in_path, 
                     header=None, 
                     sep='\t', 
                     names=colnames)
    
    df['RefSeq accession'] = df['RefSeq accession'].astype(str).str.strip()
    return df


def get_sum_zero_reads(df):
    """A function to return the number of basepairs
    with zero reads
    
    Parameters:
    -----------------------------
        : df (pd.DataFrame): dataframe with read depth information
        
    Returns:
    -----------------------------
        : zero_depth_bp (int): the number of basepairs with no read depth
    """
    zips = df[df['depth'] == 0].reset_index()
    zips['region_size'] = zips['position_end'] - zips['position_start']
    return zips['region_size'].sum()


def get_read_coverage_per_xsome(df):
    """A function to produce a dataframe with coverage estimated by
    individual basepair vs. by each fragement of the reference genome
    
    Parameters:
    -----------------------------
        : df (pd.DataFrame): dataframe with read depth information
        
    Returns:
    -----------------------------
        : coverage_df (pd.dataframe): per chromosome read coverage
    """
    new_rows = []
    
    for xsome in df['RefSeq accession'].unique():
        tmp = df[df['RefSeq accession'] == xsome]
        
        total_len = tmp['position_end'].max()
        no_coverage = get_sum_zero_reads(tmp)
        coverage = 1 - (no_coverage / total_len)
        
        mean_depth_bp = tmp['depth'].sum() / total_len
        mean_depth = tmp['depth'].mean()
        max_depth = tmp['depth'].max()
        min_depth = tmp['depth'].min()
        
        row = {
            'rname' : xsome,
            'total_length' : total_len,
            'n_bps_no_coverage' : no_coverage,
            'perc_coverage' : coverage,
            'mean_depth_bp' : mean_depth_bp,
            'mean_depth_fragment' : mean_depth,
            'max_depth' : max_depth,
            'min_depth' : min_depth
        }
        
        new_rows.append(row)
        
    coverage = pd.DataFrame(new_rows)
    return coverage

if __name__ == "__main__":
    df = load_coverage(sys.argv[1])
    coverage = get_read_coverage_per_xsome(df)
    coverage.to_csv(sys.stdout, index=False)

    