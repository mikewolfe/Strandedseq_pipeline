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
            is_a_nan = np.isnan(this_array)
            starts = np.arange(0, chrm_dict[key], res, dtype=np.int64)
            ends = np.arange(res, chrm_dict[key], res, dtype=np.int64)
            if len(ends) != len(starts):
                ends = np.append(ends, chrm_dict[key])
            names = np.array([key]*len(starts))
            outf.addEntries(names[~is_a_nan], starts[~is_a_nan], 
                    ends = ends[~is_a_nan], values=this_array[~is_a_nan])

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
        subtract_val = np.nanmean(background_vals)
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

     
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Manipulate BigWig files")
    parser.add_argument('infile', type=str, help="file to convert from")
    parser.add_argument('outfile', type=str, help="file to convert to")
    parser.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser.add_argument('--operation', type=str, default=None,
            help="operation to perform before writing out file. \
            All operations, neccesitate conversion to array internally \
            options {'RobustZ', 'Median_norm', 'background_subtract', 'scale_max'}")
    parser.add_argument('--background_regions', type=str, default=None,
            help="bed file containing known regions of background.")
    
    args = parser.parse_args()
    
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
