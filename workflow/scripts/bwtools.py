# great module from deeptools that helps deal with bigwigs. I am only 
# adding support functions here for the module
import pyBigWig
import numpy as np
import arraytools

def write_arrays_to_bigwig(outfilename, arrays, chrm_dict, res = None):
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
            outf.addEntries(key, 0, values=arrays[key], span = res, step = res)


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
    one_array = np.hstack(arrays.values())
    median = np.nanmedian(one_array)
    mad = stats.median_absolute_deviation(one_array, nan_policy='omit', scale=1.4826)
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
    one_array = np.hstack(arrays.values())
    median = np.nanmedian(one_array)
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, median)
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
            options {'RobustZ', 'Median_norm'}")
    
    args = parser.parse_args()
    
    operation_dict={"RobustZ": RobustZ_transform, 
            "Median_norm": Median_norm}

    # read in file 
    inf = pyBigWig.open(args.infile)
    
    # convert to a dictionary of arrays
    arrays = bigwig_to_arrays(inf, res = args.res)
    
    # perform operation on arrays
    arrays = operation_dict[args.operation](arrays)

    # write out file
    write_arrays_to_bigwig(args.outfile, arrays, inf.chroms(), res = args.res)
    inf.close()
