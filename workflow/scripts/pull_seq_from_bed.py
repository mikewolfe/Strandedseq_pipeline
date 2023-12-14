import fasta
import bed_utils
import operator


if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser("Python script to harvest sequences corresponding to a bed")

    parser.add_argument('bed_file', type=str, help = "bed file containing locations to get sequence for")
    parser.add_argument('fasta_file', type=str, help = "fasta file to pull sequences from")
    parser.add_argument('--upstream', type=int, help = "amount upstream to pad sequence by. Default = 0", default = 0)
    parser.add_argument('--downstream', type=int, help = "amount downstream to pad sequence by. Default = 0", default = 0)
    parser.add_argument('--circular', action = "store_true", help = "treat chromosomes as circular?")
    parser.add_argument('--filter_column', type = int, help = "0-based column number to filter on. Default = 4", default = 4)
    parser.add_argument('--filter_value', type = float, help = "cutoff to filter at. Default None (no filter)", default = None)
    parser.add_argument('--filter_direction', type = str, help = "greater or less, default = greater", default = "greater")

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

    if args.filter_value is not None:
        cutoff_value = args.filter_value
        direction_factory = {"greater": operator.ge, "less": operator.le}
        direction = direction_factory[args.filter_direction]
        filter_func = lambda entry: direction(float(entry[args.filter_column]), cutoff_value)
    else:
        filter_func = lambda entry: True


    for entry in bed:
        this_chrm = entry["chrm"]
        this_start = entry["start"]
        this_end = entry["end"]
        this_strand = entry["strand"]
        if not filter_func(entry):
            continue
        if this_strand == "-":
            this_start = this_start - downstream
            this_end = this_end + upstream 
            rc = True
        else:
            this_start = this_start - upstream
            this_end = this_end + downstream 
            rc = False
        chrm = genome.pull_entry(this_chrm)
        try:
            seq = chrm.pull_seq(this_start, this_end, rc = rc, circ = args.circular)
        except ValueError as err:
            sys.stderr.write("Skipping...\n")
            sys.stderr.write(str(err) + "\n")
            continue
        sys.stdout.write(seq + "\n")

