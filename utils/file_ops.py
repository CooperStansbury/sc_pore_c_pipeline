import sys
import os
from pathlib import Path


def get_samples(input_path):
    """A function to get a list of input sample reads from
    the input dir (reads_dir)
    
    Parameters:
    -----------------------------
        : input_path (str): the path to the inputs
        
    Returns:
    -----------------------------
        : samples (list): list of file base names
    """
    samples = []
    for f in os.listdir(input_path):
        base = os.path.basename(f)
        basename = os.path.splitext(base)[0]
        samples.append(basename)
    return samples
        
        
def build_output_dir(output_path):
    """A function to build the output directory
    
    Parameters:
    -----------------------------
        : output_path (str): the path to the new outputdir
        
    Returns:
    -----------------------------
        : NA: creates nested directory structure
    """
    Path(output_path).mkdir(parents=True, exist_ok=False)
    
    alignments_path = f"{output_path}/alignments"
    tables_path = f"{output_path}/tables"
    logs_path = f"{output_path}/logs"
    
    Path(alignments_path).mkdir()
    Path(tables_path).mkdir()
    Path(logs_path).mkdir()
    
    
    
    
    