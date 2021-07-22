import subprocess

def load_modules():
    """A function to load GL modules
    """
    bashCommand = "module load python3.8-anaconda Bioinformatics bwa samtools"
    subprocess.run(bashCommand, shell=True, universal_newlines=True, check=True)