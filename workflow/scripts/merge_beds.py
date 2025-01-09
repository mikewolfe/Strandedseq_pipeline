import bed_utils
import scipy.cluster
import scipy.spatial

# functions to augment BedEntry class
def parse_name(bedentry):
    split_name = bedentry["name"].split(";")
    genes = split_name[0]
    duplicate_number = split_name[1:]
    bedentry.genes = genes

# functions to deal with merging entries
def parse_all_names(bedfile):
    for entry in bedfile:
        parse_name(entry)

def group_entries_by_genes(bedfiles):
    gene_dict_plus = {}
    gene_dict_minus = {}
    gene_dict_unsure = {}
    for this_bed in bedfiles:
        for entry in this_bed:
            if entry["strand"] == "+":
                entry_list = gene_dict_plus.get(entry.genes, [])
                entry_list.append(entry)
                gene_dict_plus[entry.genes] = entry_list
            elif entry["strand"] == "-":
                entry_list = gene_dict_minus.get(entry.genes, [])
                entry_list.append(entry)
                gene_dict_minus[entry.genes] = entry_list
            else:
                entry_list = gene_dict_unsure.get(entry.genes, [])
                entry_list.append(entry)
                gene_dict_unsure[entry.genes] = entry_list
    return gene_dict_plus, gene_dict_minus, gene_dict_unsure

def determine_overlap(entry1, entry2):
    """ 
    Doing interval distances using the Hausdroff Distance. Which for the
    1-dimensional case is
    D(a, b) = max(abs(a1 - b1), abs(a2 - b2))
    http://www.math.u-bordeaux.fr/~mchave100p/wordpress/wp-content/uploads/2012/12/slides_IFCS04.pdf
    https://en.wikipedia.org/wiki/Hausdorff_distance
    https://math.stackexchange.com/questions/41269/distance-between-two-ranges
    """
    return max([abs(entry1[0] - entry2[0]), abs(entry1[1] - entry2[1])])

def split_to_similar(entries, min_diff = 100):
    isoform_dict = {}
    # if there is only one entry then don't bother to do anything
    if len(entries) < 2:
        isoform_dict[1] = entries
        return isoform_dict
    # convert entries into an array
    entry_array = []
    for entry in entries:
        entry_array.append([entry["start"], entry["end"]])
    # make a distance matrix
    d = scipy.spatial.distance.pdist(entry_array, lambda x, y: determine_overlap(x, y))
    # get the linkage
    z = scipy.cluster.hierarchy.linkage(d, metric = 'single')
    # make flat clusters based on the linkage
    c = scipy.cluster.hierarchy.fcluster(z, min_diff, criterion='distance')
    # make and fill dictionary of clusters
    for entry, cluster_val in zip(entries, c):
        clust_list = isoform_dict.get(cluster_val, [])
        clust_list.append(entry)
        isoform_dict[cluster_val] = clust_list
    return isoform_dict
        
def merge_isoforms(isoform_dict):
    new_entries = []
    for key in isoform_dict:
        new_entry = bed_utils.BedEntry()
        chrms = [entry["chrm"] for entry in isoform_dict[key]]
        # check that all the chrms are the same
        if not all(x == chrms[0] for x in chrms):
            raise ValueError("Tried to merge isoforms on different chromosomes")
        starts = [entry["start"] for entry in isoform_dict[key]]
        ends = [entry["end"] for entry in isoform_dict[key]]
        strands = [entry["strand"] for entry in isoform_dict[key]]
        if not all(x == strands[0] for x in strands):
            raise ValueError("Tried to merge isoforms on different strands %s"%isoform_dict)
        new_entry["chrm"] = chrms[0]
        new_entry["start"] = min(starts)
        new_entry["end"] = max(ends)
        # assumes entries had a name of the form gene1,gene2;1
        if strands[0] is not "+" and strands[0] is not "-":
            this_strand = "."
        else:
            this_strand = strands[0]
        new_entry["name"] = isoform_dict[key][0].genes + ";%d"%(key)
        new_entry["strand"] = this_strand
        new_entries.append(new_entry)
    return new_entries

def determine_isoforms(gene_dict, min_diff = 100):
    out_bed = bed_utils.BedFile()
    for genes in gene_dict:
        isoforms = split_to_similar(gene_dict[genes], min_diff = min_diff)
        new_entries = merge_isoforms(isoforms)
        for entry in new_entries:
            out_bed.add_entry(entry)
    return out_bed
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge beds")
    parser.add_argument("outfile", type=str, help="path to output file")
    parser.add_argument("infiles", type=str, nargs = "+", help="path to input fasta")
    parser.add_argument("--min_diff", type=int, default=100, help="minimum difference in basepairs to consider unique")
    args = parser.parse_args()

    beds = []
    # read in files
    for this_file in args.infiles:
        this_bed = bed_utils.BedFile()
        this_bed.from_bed_file(this_file)
        beds.append(this_bed)

    # parse the names to add genes
    for bed in beds:
        parse_all_names(bed)


    # determine groupings by names
    gene_dict_plus, gene_dict_minus, gene_dict_unsure = group_entries_by_genes(beds)

    # merge into isoforms
    out_bed_plus = determine_isoforms(gene_dict_plus, min_diff = args.min_diff)
    out_bed_minus = determine_isoforms(gene_dict_minus, min_diff = args.min_diff)
    out_bed_unsure = determine_isoforms(gene_dict_unsure, min_diff = args.min_diff)
    for entry in out_bed_minus:
        out_bed_plus.add_entry(entry)
    for entry in out_bed_unsure:
        entry["strand"] = "."
        out_bed_plus.add_entry(entry)
    out_bed_plus.sort()
    out_bed_plus.write_bed_file(args.outfile)
