import pandas as pd
import numpy as np
import sys
import os
import itertools


def get_pairs(df, chroms):
    """A function to arrange mutli-way contacts into pairs of contacts
    
    NOTE: only build pairs for contacts on chromosomes
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table
        :  chroms (list of str): all chromosomes (named)
        
    Returns:
    -----------------------------
        : pairs (pd.DataFrame): upper triangle formatted pairs table (such that the coordinate 
        of pos1 < pos2).
    """
    
    new_rows = []
    
    for read in df['read_name'].unique():
        tmp = df[df['read_name'] == read]
        
        if tmp.shape[0] > 1:
            # only include unique loci in pairs
            pairs = list(itertools.combinations(tmp['fragment_midpoint'].unique(), 2))
            
            for pair in pairs:
                
                pos1 = np.min(pair)
                pos2 = np.max(pair)
                
                chrom1 = tmp.loc[tmp['fragment_midpoint'] == pos1, 'chrom'].iloc[0]
                chrom2 = tmp.loc[tmp['fragment_midpoint'] == pos2, 'chrom'].iloc[0]
                
                # only add contact between loci on chromosomes
                if chrom1 in chroms and chrom2 in chroms:
                    new_row = {
                        'readID' : read,
                        'chr1' : chrom1,
                        'pos1' : pos1,
                        'chr2' : chrom2,
                        'pos2' : pos2,
                    }

                    new_rows.append(new_row)
                
    pairs = pd.DataFrame(new_rows)
    return pairs


def build_header(assemby_df, assembly_field):
    """A function bulid the pairs file header for sys.stdout
    
    Parameters:
    -----------------------------
        : assemby_df (pd.DataFrame): the assmebly file
        : assembly_field (str): the field of the assembluy file to use for 
        chromosome information
        
    Returns:
    -----------------------------
        : None: all information piped to sys.stdout
    """
    
    col = [c for c in assemby_df.columns if assembly_field in c][0]
    
    print('## pairs format v1.0')
    print('#shape: upper triangle')
    print('#command: build_pairs.py')
    
    for idx, row in assemby_df.iterrows():
        print(f"#chromsize {row[col]} {row['Total length (bp)']}")
    
    print('#columns: readID chr1 pos1 chr2 pos2')
    

def get_chroms(assemby_df, assembly_field):
    """A function get a list of chromosomes
    
    Parameters:
    -----------------------------
        : assemby_df (pd.DataFrame): the assmebly file
        : assembly_field (str): the field of the assembluy file to use for 
        chromosome information
        
    Returns:
    -----------------------------
        : chroms (list of str): all chromosomes (named)
    """
    col = [c for c in assemby_df.columns if assembly_field in c][0]
    return assemby_df[col].to_list()



if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1])
    assembly = pd.read_csv(sys.argv[2])
    assembly_field = sys.argv[3]
    
    chroms = get_chroms(assembly, assembly_field)
    build_header(assembly, assembly_field)
    pairs = get_pairs(df, chroms)
    pairs.to_csv(sys.stdout, index=False, mode='a')
    
    print(pairs.head())
    
    