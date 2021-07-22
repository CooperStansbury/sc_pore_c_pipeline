import sys
import pysam


def reformat_bam(infh, outfh):
    """A function to reformat a bam file QNAME to contain
    the read index, readname, and align index separated by colons.
    
    The purpose of this operation is to pipe a BAM file in, and create 
    IDS for various downstream tasks.
    
    Notes: 
        - this script is taken from: https://github.com/nanoporetech/Pore-C-Snakemake/blob/ec40f1e0dedc6992108b9efa494d4fc1aaf89fee/scripts/reformat_bam.py#L1
        -    <read_id> -> <read_id>:<read_idx>:<align_idx>
    Where 'read_idx' is a unique integer id for each read within the file and
    'align_idx' is a unique integer id for each alignment within the file. The
    tool also adds a 'BX' tag consisting of the 'read_id' to each record.
    
    Parameters:
    -----------------------------
        : infh (stream): sys.stdin of a BAM file
        
    Returns:
    -----------------------------
        : outfh (list): sys.stdout of a reformatted BAM file
    """

    infile = pysam.AlignmentFile(infh, "r")
    outfile = pysam.AlignmentFile(outfh, "w", template=infile)
    read_indices = {}
    for align_idx, align in enumerate(infile.fetch(until_eof=True)):
        read_id = align.query_name
        read_idx = read_indices.get(read_id, None)
        if read_idx is None:
            read_idx = len(read_indices)
            read_indices[read_id] = read_idx
        align.query_name = f"{read_id}:{read_idx}:{align_idx}"
        outfile.write(align)


if __name__ == "__main__":
    reformat_bam(sys.stdin, sys.stdout)