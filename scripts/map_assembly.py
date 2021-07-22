import sys
import pandas as pd 

def map_chromosome_names(align, assembly, output_csv):
    """A function to map the accession to a chromosome number and other metadata
    based on the input assembly.
    
    
    Parameters:
    -----------------------------
        : align (str): path to the alignment file
        : assembly (str): path to the assembly file
        
    Returns:
    -----------------------------
        : align (stream): the alignment table with chromosome information
    """
    df = pd.read_csv(align)
    assembly = pd.read_csv(assembly)
    
    df['chrom'] = df['chrom'].astype(str).str.strip()
    assembly['RefSeq accession'] = assembly['RefSeq accession'].astype(str).str.strip()
    
    df = pd.merge(df, assembly, 
                  left_on='chrom', 
                  right_on='RefSeq accession', 
                  how='left')
    
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    map_chromosome_names(sys.argv[1], sys.argv[2], sys.stdout)