import sys
import os
from pathlib import Path

def get_dir(directory):
    """A function to get a list of filepaths from a 
    directory
    
    Parameters:
    -----------------------------
        : directory (str): the path to the input directory
        
    Returns:
    -----------------------------
        : f_list (list): list of filepaths
    """
    f_list = []
    for f in os.listdir(directory):
        if directory.endswith("/"):
            f_path = f"{directory}{f}"
        else:
            f_path = f"{directory}/{f}"
        
        f_list.append(f_path)
        
    return f_list


def get_sample_list(input_dir):
    """A function to get a list of input sample reads from
    the input dir (reads_dir)
    
    Parameters:
    -----------------------------
        : input_dir (str): the path to the input directory
        
    Returns:
    -----------------------------
        : samples (list): list of file base names
    """
    samples = []
    for f in os.listdir(input_dir):
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