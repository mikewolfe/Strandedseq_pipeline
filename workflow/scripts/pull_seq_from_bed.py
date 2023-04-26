import fasta
import bed_utils


if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser("Python script to harvest sequences corresponding to a bed")

    parser.add_argument('bed_file', type=str, help = "bed file containing locations to get sequence for")
    parser.add_argument('fasta_file', type=str, help = "fasta file to pull sequences from")
    parser.add_argument('--upstream', type=int, help = "amount upstream to pad sequence by. Default = 0", default = 0)
    parser.add_argument('--downstream', type=int, help = "amount downstream to pad sequence by. Default = 0", default = 0)
    parser.add_argument('--circular', action = "store_true", help = "treat chromosomes as circular?")

    args = parser.parse_args()
    bed_file = args.bed_file
    fasta_file = args.fasta_file
    upstream = args.upstream
    downstream = args.downstream

    bed = bed_utils.BedFile()
    if bed_file == "-":
        bed_file = 0
    bed.from_bed_file(bed_file, header = 1)

    genome = fasta.FastaFile()
    with open(fasta_file, mode = "r") as inf:
        genome.read_whole_file(inf)

    for entry in bed:
        this_chrm = entry["chrm"]
        this_start = entry["start"]
        this_end = entry["end"]
        this_strand = entry["strand"]
        if this_strand == "-":
            this_start = this_start - downstream
            this_end = this_end + upstream 
            rc = True
        else:
            this_start = this_start - upstream
            this_end = this_end + downstream 
            rc = False
        chrm = genome.pull_entry(this_chrm)
        seq = chrm.pull_seq(this_start, this_end, rc = rc, circ = args.circular)
        sys.stdout.write(seq + "\n")

