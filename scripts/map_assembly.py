import sys
import pandas as pd 

def map_chromosome_names(align, assembly, field, output_csv):
    """A function to map the accession to a chromosome number and other metadata
    based on the input assembly.
    
    Parameters:
    -----------------------------
        : align (str): path to the alignment file
        : assembly (str): path to the assembly file
        : field (str): the assembly field to use for mapping chromosome 
        information. Can be either of:
            `GenBank`
            `RefSeq`
        
    Returns:
    -----------------------------
        : align (stream): the alignment table with chromosome information
    """
    clean_field = field.strip() + " accession"
    

    df = pd.read_csv(align)
    assembly = pd.read_csv(assembly)
    
    if not clean_field in assembly.columns:
        raise ValueError(f'Must pass either `GenBank` or `RefSeq` as the assembly_field in config.yaml. Passed {clean_field}')
    
    df['chrom'] = df['chrom'].astype(str).str.strip()
    assembly[clean_field] = assembly[clean_field].astype(str).str.strip()
    
    df = pd.merge(df, assembly, 
                  left_on='chrom', 
                  right_on=clean_field, 
                  how='left')
    
    df['fragment_midpoint'] = (df['fragment_end'] + df['fragment_start']) // 2
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    map_chromosome_names(sys.argv[1], sys.argv[2], sys.argv[3], sys.stdout)