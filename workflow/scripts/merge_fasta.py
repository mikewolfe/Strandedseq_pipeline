import fasta
import argparse
# TODO more sophisticated duplicate handling

def read_file(ff):
    this_fasta = fasta.FastaFile()
    with open(ff) as inf:
        this_fasta.read_whole_file(inf)
    return this_fasta

def merge_fasta(existing_fasta, new_fasta):
    for entry in new_fasta:
        orig_name = entry.chrm_name()
        dup_num = 0
        if orig_name in existing_fasta.chrm_names():
            new_name = orig_name + "_dup_%d"%(dup_num)
            while new_name in existing_fasta.chrm_names():
                dup_num += 1
                new_name = orig_name + "_dup_%d"%(dup_num)
            entry.set_header(">" + new_name)
            existing_fasta.add_entry(entry)
        else:
            existing_fasta.add_entry(entry)
    return existing_fasta
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge fastas")
    parser.add_argument("outfile", type=str, help="path to output file")
    parser.add_argument("infiles", type=str, nargs = "+", help="path to input fasta")
    args = parser.parse_args()
    
    all_files = [read_file(this_file) for this_file in args.infiles]

    outfasta = fasta.FastaFile()
    for this_file in all_files:
        outfasta = merge_fasta(outfasta, this_file)
    with open(args.outfile, mode = "w") as outf:
        outfasta.write(outf)
