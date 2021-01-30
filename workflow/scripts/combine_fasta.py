import fasta as fa
import sys


if __name__ == "__main__":
    genome_name = sys.argv[1]
    fasta_files = sys.argv[2:]
    all_fastas = [fa.FastaFile() for fafile in fasta_files]
    # read in all fastas
    for fafile, fafilename in zip(all_fastas, fasta_files):
        with open(fafilename, mode = "r") as inf:
            fafile.read_whole_file(inf)

    final_fasta = fa.FastaFile()
    for fafile in all_fastas:
        for entry in fafile:
            final_fasta.add_entry(entry)
    lengths = {}
    for entry in final_fasta:
        chrm_name = entry.chrm_name()
        lengths[chrm_name] = len(entry)

    with open(genome_name + "_contig_sizes.tsv", mode = "w") as outf:
        for chrm, length in lengths.items():
            outf.write("%s\t%s\n"%(chrm, length))
    with open(genome_name + ".fa", mode = "w") as outf:
        final_fasta.write(outf)




