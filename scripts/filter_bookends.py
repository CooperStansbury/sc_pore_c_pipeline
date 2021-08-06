import pandas as pd
import numpy as np
import sys
import os


def add_fragment_midpoints(df):
    """A function to add a column of fragment midpoints. These coordinates are reference genome 
    positions halfway (rounding up) along the mapped fragment. This decision was made by Dr. Can Chen
    based on the Nanopore logic in pore-C snakemake during generation of the pairs table.
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the 'raw' alignment table
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table with 1 additional column
    """
    df['fragment_mean'] = np.ceil((df['fragment_start'] + df['fragment_end']) /2 )
    df['fragment_mean'] = df['fragment_mean'].astype(int)
    return df


def drop_low_fragment_count_reads(df, n=1, verbose=True):
    """A function to drop reads that have n or fewer fragments 
        
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        :  n (int): reads with n or fewer fragments will be dropped
        :  verbose (bool): if True, print numner of rows
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table without low fragement reads
    """
    grped  = df.groupby('read_name', as_index=False).agg({
        'fragment_id' : 'count',
    })
    if verbose:
        print(f"\nn reads BEFORE striping {n}-fragment reads: {len(grped)}")
    
    # filter 
    grped = grped[grped['fragment_id'] > n]
    df = df[df['read_name'].isin(grped['read_name'])]
    
    if verbose:
        print(f"n reads AFTER striping {n}-fragment reads: {len(grped)}")
        
    return df


def per_read_filter(df, criterion, verbose=0):
    """A function to filter individual reads and remove duplicate with-in read 
    fragments
    
    Parameters:
    -----------------------------
      :  df (pd.DataFrame): the alignment table
      :  criterion (str): column name of the column to use for determining which 
      fragment (if duplciated), should be kept. The logic takes the fragment with 
      the HIGHEST value for the criterion specified. 
      :  verbose (int): one of 0, 1, 2
          0: show nothing
          1: print before and after alignment dataframe sizes (rows)
          2: print individual read filter results (debugging only) 
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the 'raw' alignment table
    """
    df_list = []
    
    if verbose > 0:
        print(f"\nalignment table total fragments: {len(df)}")
    
    for read_idx in df['read_name'].unique():
        read_df = df[df['read_name'] == read_idx].reset_index()
        
        read_df['read_min'] = read_df['fragment_start'].min() # add first fragment
        read_df['read_max'] = read_df['fragment_end'].max() # add last fragment
        read_df = read_df.sort_values(by=['fragment_id', criterion], ascending=False)
        
        before = len(read_df)
        read_df = read_df.drop_duplicates(subset=['fragment_id'], keep='first')
        after = len(read_df)
        read_df['n_fragments'] = len(read_df)
        
        first_fragment = read_df[read_df['fragment_start'] == read_df['read_min']]['fragment_id'].unique()[0]
        last_fragment = read_df[read_df['fragment_end'] == read_df['read_max']]['fragment_id'].unique()[0]
        read_df['first_fragment'] = first_fragment
        read_df['last_fragment'] = last_fragment
        
        if verbose == 2:
            print(read_idx, before, after)
        
        df_list.append(read_df)
        
    res = pd.concat(df_list, ignore_index=True)
    
    if verbose > 0:
        print(f"alignment table fragments after filtering reads: {len(res)}")
    return res


def get_maximal_reads(df, n=2):
    """A function to return the n reads that share bookending 
    fragment IDS, sorted by mean perc_of_alignment (read)
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table
        :  n (int): the number of top reads to take (biologically infeasible to take > 4)
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table with maximal reads
    """
    
    # reduce rows with unique fragments 
    # - since we perform this excersise on reads, not fragments
    grped = df.groupby('read_name', as_index=False).agg({
        'first_fragment' : 'first', 
        'last_fragment' : 'first', 
        'n_fragments' : 'first',
        'perc_of_alignment' : np.mean
    })
    
    # NOTE: we may need to do some cross-read fragment
    # decision logic here
    
    # sort by number of fragments and mapping quality
    grped = grped.sort_values(by=['first_fragment', 'last_fragment', 'n_fragments', 'perc_of_alignment'], ascending=False)
    
    # add an identifier for reads that have the same start and end fragment 
    N_TOP = grped.groupby(['first_fragment', 'last_fragment']).cumcount()
    grped.loc[:, 'N_TOP'] = N_TOP + 1
    
    # compute the maximal number of fragments for reads sharing bookended 
    # fragment IDS
    grped['matching_group_max'] = df.groupby(['first_fragment', 'last_fragment'])["n_fragments"].transform(np.max)
    
    # set a flag to take the top N reads with the highest number of fragements
    mask = (grped['n_fragments'] == grped['matching_group_max']) & (grped['N_TOP'] <= n)
    grped['SELECT'] = np.where(mask, 1, 0)

    grped = grped[grped['SELECT'] == 1]
    read_ids = grped['read_name']
    
    # filter the original data frame for only those reads
    df = df[df['read_name'].isin(read_ids)]
    return df


if __name__ == "__main__":
    input_path = sys.argv[1]
    criterion = 'perc_of_alignment'
    filter_n_fragments = 1
    n_top_reads = 2
    
    df = pd.read_csv(input_path)
    df = add_fragment_midpoints(df) # compute mid points
    df = per_read_filter(df, criterion, verbose=0) # filter duplicate within read fragments 
    df = drop_low_fragment_count_reads(df, n=filter_n_fragments, verbose=False) # drop reads with low fragment counts
    df = get_maximal_reads(df, n_top_reads) # filter for reads spanning the same fragments
    
    df.to_csv(sys.stdout, index=False)