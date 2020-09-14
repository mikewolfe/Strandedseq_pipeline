from Bio import SeqIO
import sys
import argparse

#
#gi|48994873|gb|U00096.2|mod|ATCC.47076| - 1..4639676
#233 RNAs
#Location     

# example ptt file
# Pseudomonas aeruginosa PAO1 chromosome, complete genome. - 0..6264404
# 5572 proteins
# Location    Strand  Length  PID Gene    Synonym Code    COG Product
# 483..2027   +   514 -   dnaA    PA0001  -   chromosome replication initiator DnaA
# 2056..3159  +   367 -   dnaN    PA0002  -   DNA polymerase III subunit beta
# 3169..4278  +   369 -   recF    PA0003  -   DNA replication and repair protein RecF
# 4275..6695  +   806 -   gyrB    PA0004  -   DNA gyrase subunit B
# 7018..7791  -   257 -   lptA    PA0005  -   lysophosphatidic acid acyltransferase
# 7803..8339  -   178 -   -   PA0006  -   D-glycero-beta-D-manno-heptose-1,7-bisphosphate 7-phosphatase
# 8671..10377 +   568 -   -   PA0007  -   hypothetical protein
# example rnt file
# Bacteria name - 1..6264404
# 106 RNAs (number sRNAs)
# Location Strand Length PID Gene Synonym Code COG Product
# 298816..298892 - 77 110645304 - PA0263.1 - - Arg tRNA

def convert_strandval(strandval):
    if strandval == 1:
        out = "+"
    elif strandval == -1:
        out = "-"
    else:
        out = "NA"
    return out

def out_ptt_rnt(gb, ftype = "ptt", chrm = None):
    # rockhopper only parses the genome name if it has >=5 fields split on
    # pipes. The 0-indexed [3] field is the name given to *_transcripts when
    # split on a "."
    # The [4] field is another name which I don't know where it goes
    if chrm:
        outchrm = str(chrm)
    else:

        outchrm = str(gb.id)
    if ftype == "rnt":
        features = ["tRNA", "rRNA", "sRNA"]
        name =  "RNAs"
        unknown_name = "NA_RNA"
        length_func = lambda start, end: end - start
    elif ftype == "ptt":
        features = ["CDS"]
        name = "proteins"
        unknown_name = "NA_Prot"
        length_func = lambda start, end: int(((end - start)/3) - 1)
    else:
        raise ValueError("Filetype %s not recognized"%ftype)

    all_features = []
    for NA_num, feature in enumerate(gb.features):
        outstr = ""
        if feature.type in features:
        # in one-based coordinates need to adjust
            outstr += str(feature.location.start + 1) + ".." + str(feature.location.end) + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += str(length_func(feature.location.start,feature.location.end)) + "\t"
            outstr += str(feature.qualifiers.get("protein_id", ["-"])[0]) + "\t"
            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = str(feature.qualifiers.get("locus_tag", ["%s_%d"%(unknown_name,NA_num)])[0])
            outstr += str(feature.qualifiers.get("locus_tag", [unique_name])[0]) + "\t"
            outstr += str(feature.qualifiers.get("gene", [unique_name])[0]) + "\t"
            outstr += str(feature.qualifiers.get("code", ["-"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("cog", ["-"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("product", ["-"])[0]) 
            all_features.append(outstr)
    sys.stdout.write("%s - 1..%d\n"%(outchrm, len(gb.seq)))
    sys.stdout.write("%d %s\n"%(len(all_features), name))
    sys.stdout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
    for feature in all_features:
        sys.stdout.write(feature + "\n")



def out_tsv(gb, chrm = None):
    if chrm:
        outchrm = str(chrm)
    else:
        outchrm = str(gb.name)
    sys.stdout.write("chr\tstart\tend\tstrand\tlocus_tag\tprotein_id\tproduct\tseq\ttranslation\n")
    for feature in gb.features:
        outstr = ""
        if feature.type == "CDS":
            outstr += outchrm + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += str(feature.qualifiers.get("locus_tag", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("protein_id", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("product", ["NA"])[0]) + "\t"
            # pull the nucleotide sequence
            nuc_seq = feature.extract(gb.seq)
            this_translation = nuc_seq.translate(table = "Bacterial", cds = True)
            translation = str(feature.qualifiers['translation'][0])
            # check that the translation matches the expected product. Have to drop the stop codon in
            # the nucleotide translation
            if not (this_translation == translation):
                raise ValueError("Nucleotide sequence doesn't match translation\n%s\n%s"%(this_translation, translation))
            outstr += str(nuc_seq) + "\t"
            outstr += translation + "\n"
            sys.stdout.write(outstr)


def out_bed(gb, features = ["CDS"], chrm = None):
    if chrm:
        name = str(chrm)
    else:
        name = str(gb.name)
    for NA_num, feature in enumerate(gb.features):
        if feature.type in features:
            outstr = ""
            outstr += name + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = str(feature.qualifiers.get("locus_tag", ["%s_%d"%("NA",NA_num)])[0])
            outstr += str(feature.qualifiers.get("locus_tag", [unique_name])[0]) + "\t"
            outstr += "." + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\n"
            sys.stdout.write(outstr)

def out_fasta(gb, features = ["CDS"]):
    import fasta as fa
    out_fasta = fa.FastaFile()
    for NA_num, feature in enumerate(gb.features):
        if feature.type in features:
            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = str(feature.qualifiers.get("locus_tag", ["%s_%d"%("NA",NA_num)])[0])
            header = ">" + str(feature.qualifiers.get("locus_tag", [unique_name])[0]) + \
            " " +  \
            str(feature.qualifiers.get("product", ["NA"])[0])
            seq = str(feature.extract(gb.seq))
            entry = fa.FastaEntry(header = header, seq = seq)
            out_fasta.add_entry(entry)
    out_fasta.write(sys.stdout)

def out_fna(gb, chrm = None):
    import fasta as fa
    out_fasta = fa.FastaFile()
    if chrm:
        header = ">" + str(chrm)
    else:
        header = ">" + str(gb.name)
    seq = str(gb.seq)
    entry = fa.FastaEntry(header = header, seq = seq)
    out_fasta.add_entry(entry)
    out_fasta.write(sys.stdout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a genbank file for features")
    parser.add_argument("infile", type=str, help="path to input file")
    parser.add_argument("--outfmt", type=str, help='''output format. default = bed. 
                                                    Supported: [fa, fasta, bed, tsv, ptt, rnt, fna]
                                                    fa or fasta give each features nucleotide sequence as a fasta record.
                                                    ptt and rnt give format files for rockhopper input.
                                                    bed gives each feature in bed format.
                                                    fna pulls the full nucleotide sequence.
                                                    Note that tsv is a special format
                                                    that only pulls CDS's and their translation''',
                        default = "bed")
    parser.add_argument('--features', type=str, nargs="+", help='''
    feature types to parse for. Can specify multiple. default = "CDS"
    ''',
            default = "CDS")

    parser.add_argument('--chrm', type=str, help='''
    specify the chromosome name for the output
    ''', default = None)
    args = parser.parse_args()

    fname = args.infile
    ftype = args.outfmt
    features = args.features
    chrm = args.chrm
    with open(fname, mode = "r") as inf:
        gb = SeqIO.read(inf, "genbank")
    if ftype == "fasta" or ftype == "fa":
        out_fasta(gb, features)
    elif ftype == "bed":
        out_bed(gb, features, chrm = chrm)
    elif ftype == "ptt" or ftype == "rnt":
        out_ptt_rnt(gb, ftype, chrm = chrm)
    elif ftype == "fna":
        out_fna(gb, chrm = chrm)
    else:
        out_tsv(gb, chrm = chrm)
