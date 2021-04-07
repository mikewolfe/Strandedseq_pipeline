# great module from deeptools that helps deal with bigwigs. I am only 
# adding support functions here for the module
import pyBigWig
import numpy as np
import arraytools

def write_arrays_to_bigwig(outfilename, arrays, chrm_dict, res = 1):
    """
    Convert a set of arrays for each contig to a bigwig file
    
    Args:
        outfilename - (str) output file name
        arrays - (dict) a dictionary of numpy arrays
        bw - (dict) a dictionary of chromosome lengths
        res - resolution that the array is in
    Returns:
        Writes a bigwig file
    """
    with pyBigWig.open(outfilename, "w") as outf:
        chrm_list = [(key, chrm_dict[key]) for key in chrm_dict.keys()]
        outf.addHeader(chrm_list)
        for key in chrm_dict.keys():
            this_array = arrays[key]
            # make sure nans don't get added to the bw file
            the_finite = np.isfinite(this_array)
            starts = np.arange(0, chrm_dict[key], res, dtype=np.int64)
            ends = np.arange(res, chrm_dict[key], res, dtype=np.int64)
            if len(ends) < len(starts):
                ends = np.append(ends, chrm_dict[key])
            names = np.array([key]*len(starts))
            outf.addEntries(names[the_finite], starts[the_finite], 
                    ends = ends[the_finite], values=this_array[the_finite])

def bigwig_to_arrays(bw, res = None):
    """
    Convert a bigwig to a dictionary of numpy arrays, one entry per contig
    
    Args:
        bw - pyBigWig object
        res - resolution that you want the array in
    Returns:
        arrays (dict) - a dictionary of numpy arrays
    """
    arrays = {}
    for chrm in bw.chroms().keys():
        arrays[chrm] = contig_to_array(bw, chrm, res)
    return arrays


def contig_to_array(bw, chrm, res = None):
    """
    Convert single basepair bigwig information to a numpy array
    

    Args:
        bw - a pyBigWig object
        chrm - name of chromosome you want 
        res - resolution you want data at in bp. 

    Returns:
        outarray - numpyarray at specified resolution

    """
    chrm_length = bw.chroms(chrm)
    # makes an array at 1 bp resolution
    out_array = bw.values(chrm, 0, chrm_length, numpy=True)
    if res:
        out_array = out_array[::res]
    return out_array


# Normalization/Transformation functions

def RobustZ_transform(arrays):
    """
    Transform data by the median and MAD of all contigs
    

    Args:
        arrays (dict) - dictionary of numpy arrays
    Returns:
        outarrays - dictionary of numpy arrays

    """
    from scipy import stats
    one_array = np.hstack(list(arrays.values()))
    median = np.nanmedian(one_array)
    mad = stats.median_abs_deviation(one_array, nan_policy='omit', scale='normal')
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], median, mad)
    return arrays

def Median_norm(arrays):
    """
    Normalize data by the median value at all contigs
    

    Args:
        arrays (dict) - dictionary of numpy arrays
    Returns:
        outarrays - dictionary of numpy arrays

    """
    one_array = np.hstack(list(arrays.values()))
    median = np.nanmedian(one_array)
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, median)
    return arrays

def background_subtract(arrays, background_regions = None, res = 1):
    if background_regions:
        import bed_utils
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(background_regions)
        background_vals = []
        for region in inbed:
            background_vals.extend(arrays[region["chrm"]][region["start"]//res:region["end"]//res])
        background_vals = np.array(background_vals)
        subtract_val = np.mean(background_vals[np.isfinite(background_vals)])
        for chrm in arrays.keys():
            arrays[chrm] = arraytools.normalize_1D(arrays[chrm], subtract_val, 1)
    return arrays

def scale_max(arrays, num_top, res = 1):
    one_array = np.hstack(list(arrays.values()))
    one_array = one_array[np.isfinite(one_array)]
    ind = np.argpartition(one_array, -num_top//res)[-num_top//res:]
    scale_factor = np.nanmedian(one_array[ind])

    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, scale_factor)
    return arrays



def manipulate_main(args):

    operation_dict={"RobustZ": RobustZ_transform, 
            "Median_norm": Median_norm,
            "background_subtract": lambda x: background_subtract(x, args.background_regions, args.res),
            "scale_max": lambda x: scale_max(x, 1000, args.res)}

    # read in file 
    inf = pyBigWig.open(args.infile)
    
    # convert to a dictionary of arrays
    arrays = bigwig_to_arrays(inf, res = args.res)
    
    # perform operation on arrays
    arrays = operation_dict[args.operation](arrays)

    # write out file
    write_arrays_to_bigwig(args.outfile, arrays, inf.chroms(), res = args.res)
    inf.close()

def read_multiple_bws(bw_files, res = 1):
    all_bws = {}
    open_fhandles = [pyBigWig.open(fname) for fname in bw_files]
    for fname, handle in zip(bw_files, open_fhandles):
        print(fname)
        all_bws[fname] = bigwig_to_arrays(handle, res = res)
    return (all_bws, open_fhandles)
        
def query_main(args):
    import bed_utils

    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.regions)
    res = args.res
    all_bws, open_fhandles = read_multiple_bws(args.infiles, res = res)
    
    if args.samp_names:
        samp_to_fname = {samp_name : fname for fname, samp_name in zip(args.infiles, args.samp_names)}
        samp_names = args.samp_names
    else:
        samp_to_fname = {fname : fname for fname in args.infiles}
        samp_names = args.infiles
    
    outvalues = {fname: [] for fname in args.infiles}
    region_names = []
    coordinates = []
    for region in inbed:
        for fname in args.infiles:
            these_arrays = all_bws[fname]
            array_len = len(these_arrays[region["chrm"]])
            if region["strand"] == "-":
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
    
            these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
            outvalues[fname].extend(these_values)
        these_coordinates = list(range((left_coord//res)*res,\
                (right_coord//res)*res, res))
        coordinates.extend(these_coordinates)
        region_names.extend([region["name"]]*len(these_coordinates))

    with open(args.outfile, mode = "w") as outf:
        header = "region\tcoord\t%s"%("\t".join(samp_names))
        outf.write(header + "\n")
        for i, (region, coord) in enumerate(zip(region_names, coordinates)):
            values = "\t".join([str(outvalues[samp_to_fname[samp]][i]) for samp in samp_names])
            outf.write("%s\t%s\t%s\n"%(region, coord, values))
    for fhandle in open_fhandles:
        fhandle.close()

 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Tools to work with BigWig files")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    # manipulate verb
    parser_manipulate = subparsers.add_parser("manipulate", help = "manipulate BigWig files")
    parser_manipulate.add_argument('infile', type=str, 
            help="file to convert from")
    parser_manipulate.add_argument('outfile', type=str, help="file to convert to")
    parser_manipulate.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_manipulate.add_argument('--operation', type=str, default=None,
            help="operation to perform before writing out file. \
            All operations, neccesitate conversion to array internally \
            options {'RobustZ', 'Median_norm', 'background_subtract', 'scale_max'}")
    parser_manipulate.add_argument('--background_regions', type=str, default=None,
            help="bed file containing known regions of background.")
    parser_manipulate.set_defaults(func=manipulate_main)


    # query verb
    parser_query = subparsers.add_parser("query", help = "Lookup values in BigWig files")
    parser_query.add_argument("outfile", type=str, help = "file to output to")
    parser_query.add_argument("infiles", nargs = "+", type=str, help = "input files")
    parser_query.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_query.add_argument('--regions', type=str, help = "regions to grab data from")
    parser_query.add_argument('--upstream', type=int, help = "bp upstream to add")
    parser_query.add_argument('--downstream', type = int, help = "bp downstream to add")
    parser_query.add_argument('--samp_names', type = str, nargs = "+", help = "sample names for each file")
    parser_query.set_defaults(func=query_main)
    
    args = parser.parse_args()
    args.func(args) 
