import sys
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO


def pre_digest(infile, outfile):
    """A function to digest the reads by a cut site prior to mapping 
    against the reference.
    
    This is adapted from:
    
    https://github.com/dekkerlab/MC-3C_scripts
    
    Parameters:
    -----------------------------
        : infile (str): path to the input fastq file
        : outfile (str): path to the output fastq file
    """
    output_handle = open(outfile, "w")
    handle = open(infile, "r")
    records = list(SeqIO.parse(handle, "fastq"))

    handle.close()

    for record in records:
        restriction_site_len = len( Restriction.NlaIII.site)
        cut_sites = Restriction.NlaIII.search(record.seq)

        if len(cut_sites) == 0:
            output_handle.write(record.format("fastq"))
        else :
            subread_counter = 0

            for i in range(0, len(cut_sites)):
                if i == 0:
                    subread = record[0 : cut_sites[i]-1]
                else:
                    subread = record[cut_sites[i-1]-1 : cut_sites[i]-1]

                if len(subread.seq) <= 4:
                    continue

                subread.id += "/0" + str(subread_counter)
                subread.name += "/0" + str(subread_counter)
                subread.description += "/0" + str(subread_counter)
                subread_counter += 1
                output_handle.write(subread.format("fastq"))

            subread = record[cut_sites[i]-1 :]
            subread.id += "/0" + str(subread_counter)
            subread.name += "/0" + str(subread_counter)
            subread.description += "/0" + str(subread_counter)

            if len(subread.seq) <= 4:
                continue

            output_handle.write(subread.format("fastq"))
    output_handle.close()

if __name__ == "__main__":
    pre_digest(sys.argv[1], sys.argv[2])