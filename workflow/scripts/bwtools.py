# great module from deeptools that helps deal with bigwigs. I am only 
# adding support functions here for the module
import pyBigWig
import numpy as np

def bigwig_to_array(bw, chrm, res = None):
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

def array_to_bigwig(array, bw, chrm_list, res=1):
    """
    Convert a 1D numpy array to a bigwig. Support only for 1 chromosome files
    to start
    
    Args:
        array - numpy array of values
        bw - open pyBigWig object
        chrm_list - list of (chrm_name, length) tuples
        res - resolution that the array is in
    Returns:
        nothing, edits pyBigWig object in place
    """
    bw.addHeader(chrm_list)
    chrm_name, chrm_length = chrm_list[0]
    # determine starts
    starts = np.arange(0, chrm_length, res, dtype=np.int64)
    # determine ends
    ends = np.arange(res, chrm_length, res, dtype=np.int64)
    if len(ends) != len(starts):
        # yes this is inefficient, but fixes end effect
        ends = np.append(ends, chrm_length)
    # make array for chromosome names
    names = np.array([chrm_name]*len(starts))
    bw.addEntries(names, starts, ends=ends, values=array)
     
if __name__ == "__main__":
    import argparse
    import arraytools
    parser = argparse.ArgumentParser("Convert to and from bigwig files")
    parser.add_argument('infile', type=str, help="file to convert from")
    parser.add_argument('outfile', type=str, help="file to convert to")
    parser.add_argument('--fr', type=str, default="bigwig", help="file format convert from")
    parser.add_argument('--to', type=str, default="numpy", help="file format converting to")
    parser.add_argument('--chrm_length', type=int, default=None,
            help="length of chromosome. Needed for conversions to bigwigs")
    parser.add_argument('--res', type=int, default=1,
            help="resolution. Needed for array operations")
    parser.add_argument('--chrm_name', type=str, default=None,
            help="name of chromosome. Needed for conversions to bigwigs")
    parser.add_argument('--operation', type=str, default=None,
            help="operation to perform before writing out file. \
            All operations, neccesitate conversion to array internally \
            options {'RobustZ', 'Median_norm'}")
    
    args = parser.parse_args()
    
    operation_dict={"RobustZ": arraytools.robustz_1D, 
            "Median_norm": lambda a : arraytools.normalize_1D(a, 0, np.nanmedian(a))}

    # read in file 
    if args.fr == "numpy":
        outarray = np.load(args.inf)
    elif args.fr == "bigwig":
        inf = pyBigWig.open(args.infile)
        outarray = bigwig_to_array(inf, args.chrm_name, args.res)
        inf.close()
 
    if args.operation:
        outarray = operation_dict[args.operation](outarray)

    # write out file
    if args.to == "bigwig":
        out = pyBigWig.open(args.outfile, "w")
        array_to_bigwig(outarray, out, [(args.chrm_name, args.chrm_length)], args.res)
        out.close()
    elif args.to == "numpy":
        np.save(args.outfile, outarray)
    else:
        raise ValueError("--to %s filetype not supported"%(args.to))
