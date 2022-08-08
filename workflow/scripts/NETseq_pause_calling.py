import bwtools
import numpy as np
import scipy.stats as sci_stats
from statsmodels.stats.multitest import fdrcorrection
import multiprocessing as mp
import sys

class GenomicWindow(object):
    '''
    Store information about a window in the genome

    Attributes:
        coord (int): Absolute genomic coordinate of window center
        window (np.array) : array of int values
        center_value (int): value at the center of the array
    '''
    def __init__(self, coord, array):
        self.coord = coord
        self.window = array
        self.center_value = self.window[int(len(array)/2)]

    def larson_method(self):
        """ Calculate rarity of seeing center value 
        
        Expects that the values in the window are normally distributed

        Returns:
            A tuple of information about the test including:
            coord, coord + 1,
            center_value, mean of window, stdev of window, pvalue = (1-cdf
            normal(mean, stdev))
        """
        mean = np.nanmean(self.window)
        stdev = np.nanstd(self.window)
        pvalue = sci_stats.norm.sf(self.center_value, mean, stdev)
        return(self.coord, self.coord + 1, self.center_value, mean, stdev,  pvalue)

    def bootstrap_method(self, n_boot = 10000, rng = np.random.default_rng()):
        """ Calculate rarity of seeing center value 

        Assumes that the total count of reads is distributed according to a 
        multinomial under the null.

        Args:
            n_boot (int): number of random samples from a multinomial to draw
            rng (np.random.Generator): state needed for random number generation

        Returns:
            A tuple of information about the test including:

            coord, coord + 1, center_value, expected count, total count in
            window, total random samples where max in sample >= actual
            center_value, pvalue = (total_greater_eq + 1) / (n_boot + 1)
        
        """
        # get the total number of counts
        total = self.window.sum()
        # redistribute equally using a multinomial for n_boot times
        samples = rng.multinomial(total, [1/len(self.window)]*len(self.window), size = n_boot)
        # check the number of samples where the max is higher than the actual
        total_greater_eq = np.sum(np.max(samples, axis=1) >= self.center_value)
        # compute the p value. Use a prior of one sample being greater than 
        # actual
        pvalue = (total_greater_eq + 1)/ (n_boot + 1)
        return (self.coord, self.coord+1, self.center_value, total/len(self.window), total, total_greater_eq, pvalue)

    def poisson_method(self):
        """ Calculate rarity of seeing center value

        Assumes that the reads are distributed according to a poisson
        distribution under the null. Lamba is determined as the average read
        count in the window

        Returns:
            A tuple of information about the test including:

            coord, coord + 1, center_value, expected count, total count in
            window, pvalue = 1 - cdf poisson(avg_read_count)

        """
        total = self.window.sum()
        mu = total/len(self.window)
        pvalue = sci_stats.poisson.sf(self.center_value, mu)
        return (self.coord, self.coord+1, self.center_value, mu, total, pvalue)

    def max_filter(self):
        """
        Tests if the center value of the window is the max of the window

        Returns:
            True if center value is the max value
        """
        return self.window.max() <= self.center_value

    def coverage_filter(self, cutoff):
        """
        Tests if the center value of the window is greater than a cutoff

        Args:
            cutoff (int): coverage cutoff

        Returns:
            True if greater than cutoff
        """
        return self.center_value > cutoff

    def avg_coverage_filter(self, cutoff):
        """
        Tests if the average value of the window is greater than a
        cutoff

        Args:
            cutoff (int): coverage cutoff

        Returns:
            True if greater than cutoff
        """
        return np.nanmean(self.window) > cutoff

    def zero_filter(self, cutoff):
        """
        Tests if the total fraction of zeros is less than a specified
        cutoff

        Args:
            cutoff (float): fraction of zeros that is no longer tolerable

        Returns:
            True if greater than cutoff
        """
        return np.sum(self.window == 0)/len(self.window) < cutoff

def window_generator(array, wsize, circ = True, offset = 0):
    """ Take an array and generate windows from it

    Args:
        array (1D np.array): array to create windows from
        wsize (int): half the size of the window. Placed on either side of
                     center
        circ (bool): Treat the array as circular?

    Yields:
        GenomicWindow instance
    """
    # Augment size of the array by going around the horn
    if circ:
        out_array = np.concatenate((array[-wsize:],array,array[:wsize]))
        # want wsize on either side of the center so need to add one on the right end.
        # This centers the window at center w/ wsize on either side. Total window
        # size should be 2*wsize + 1
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center -wsize + offset, out_array[center-wsize:center+wsize+1])
    else:
        out_array = array
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center + offset, out_array[center-wsize:center+wsize+1])

class RegionWindow(object):
    '''
    Store information about a region in the genome

    Attributes:
        chrm (str): chromsome region is on
        start (int): absolute start coordinate
        end (int): absolute end coordinate
        strand (str) : strand of region
        name (str): name of region
        array (np.array): array of values for the region
    '''
    def __init__(self, bed_entry):
        self.chrm = bed_entry["chrm"]
        self.start = bed_entry["start"]
        self.end = bed_entry["end"]
        self.strand = bed_entry["strand"]
        self.name = bed_entry["name"]

    def __iter__(self):
        """
        Iterates over the region from start to end.

        Returns:
            tuple(array index, absolue coordinate, relative coordinate)
        """
        for index in np.arange(0, len(self.array)):
            if self.strand == "-":
                rel_coord = len(self.array) - index - 1
            else:
                rel_coord = index
            abs_coord = self.start + index
            yield((index, abs_coord, rel_coord))

    def add_array(self, arrays):
        """
        Adds an array to the attributes of the region

        Args:
            array (np.array): np array of values typically a view
        """
        self.array = arrays[self.chrm][self.start:self.end]

    def larson_method(self):
        """ Calculate rarity of seeing each value in the region
        
        Expects that the values in the region are normally distributed

        Returns:
            A list of tuples, one for each bp in the region.

            Each tuple is in bed format with extra fields with information about
            the test. These include relative coordinate, mean, stdev
            a pvalue (1-cdf normal(mean, stdev))
        """
        out = []
        mean = np.nanmean(self.array)
        stdev = np.nanstd(self.array)
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            pvalue = sci_stats.norm.sf(this_value, mean, stdev)
            out.append((self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, mean, stdev, pvalue))
        return out

    def poisson_method(self): 
        """ Calculate rarity of seeing each value in the region
        
        Expects that the values in the region are poisson distributed

        Returns:
            A list of tuples, one for each bp in the region. 

            Each tuple is in bed format with extra fields with information about
            the test. These include relative coordinate, poisson mu, and 
            a pvalue (1-cdf poisson(mu))
        """
        out = []
        total = self.array.sum()
        mu = total/len(self.array) 
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            pvalue = sci_stats.poisson.sf(this_value, mu)
            out.append((self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, mu, pvalue))
        return out
 
    def bootstrap_method(self, n_boot = 10000, rng = np.random.default_rng()):
        """ Calculate rarity of seeing each value in the region
        
        Expects that the values in the region are randomly distributed 
        according to a multinomial model

        Returns:
            A list of tuples, one for each bp in the region.

            Each tuple is in bed format with extra fields with information about
            the test. These include relative coordinate, exp_count, total
            bootstrap samples with values >= the bp a pvalue (total_greater
            + 1 / n_boot + 1)
        """
        out = []
        # get the total number of counts
        total = self.array.sum()
        # redistribute equally using a multinomial for n_boot times
        samples = rng.multinomial(total, [1/len(self.array)]*len(self.array), size = n_boot)
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            total_greater_eq = np.sum(np.max(samples, axis=1) >= this_value) 
            pvalue = (total_greater_eq + 1)/ (n_boot + 1)
            out.append((self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, total/len(self.array), total_greater_eq, pvalue))
        return out

    def avg_coverage_filter(self, cutoff):
        """
        Check if the region passes an average coverage cutoff

        Args:
            cutoff (float): average the region bust be over
        Returns:
            True if region average is > cutoff
        """
        return np.nanmean(self.array) > cutoff

    def zero_filter(self, cutoff):
        """
        Check if the region has an acceptable number of zeros

        Args:
            cutoff (float): fraction of zeros that is no longer acceptable
        Returns:
            True if region frac_zeros is < cutoff
        """
        return np.sum(self.array == 0)/len(self.array) < cutoff

def bed_regions(arrays, bed_object):
    """ Generate regions using a dictionary of arrays and a bed file

    Args:
        arrays (dict of 1D np.array): arrays to find regions in
        bed_object (bed_utils.BedFile): A set of regions

    Yields:
        RegionWindow instance
    """

    for region in bed_object:
        out = RegionWindow(region)
        out.add_array(arrays)
        yield out

def augment_with_random_state(iterable, starting_seed):
    """ Zip an iterator with a np random seed

    This is needed to ensure that each process does not inherit the same
    randomstate during multiprocessing

    Args:
        iterable: anything you need to iterate
        starting_seed (int): starting seed for the rng

    Returns
        a list of (rng, iterable)
    """
    return [(entry, np.random.default_rng(starting_seed + index)) for index, entry in enumerate(iterable)]

def multiprocess_bootstrap(window, rng):
    return window.bootstrap_method(rng = rng)

def multiprocess_poisson(window, rng):
    return window.poisson_method()

def multiprocess_larson(window, rng):
    return window.larson_method()

def genomewide_main(args):
    
    # read in bw files
    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
    minus_strand = bwtools.pyBigWig.open(args.infile_minus)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1)

    # function factory for different methods
    methods = {'bootstrap': multiprocess_bootstrap,
            'Poisson': multiprocess_poisson,
            'Larson': multiprocess_larson}

    headers = {'bootstrap': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tsamp_gr_eq\tpvalue\tqvalue\n",
               'Poisson': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tpvalue\tqvalue\n",
               'Larson': "#chrm\tstart\tend\tname\tcount\tstrand\tmean_count\tstd_count\tpvalue\tqvalue\n"}

    # function factory for different filters
    filters = {'zero': lambda window: window.zero_filter(args.zero_cutoff),
               'avg_cov': lambda window: window.avg_coverage_filter(args.avg_cov_cutoff),
               'cov': lambda window: window.coverage_filter(args.cov_cutoff),
               'max': lambda window: window.max_filter()}

    # combine the filters into a filter list
    if len(args.filters) > 0:
        all_filters = [filters[val] for val in args.filters]
    else:
        all_filters = [lambda window: True]
    
    sys.stdout.write(headers[args.method])
    for chrm, array in arrays_plus.items():
        pool = mp.Pool(args.p)
        # allow the windows to be distributed over the processes
        output_plus = pool.starmap(methods[args.method], 
                # add a new random seed for each process. Otherwise will use the
                # same seed for every process
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    window_generator(array, args.wsize, circ = args.circular)), args.seed))
        pool.close()
        pool.join()

    for chrm, array in arrays_minus.items():
        pool = mp.Pool(args.p)
        # allow the windows to be distributed over the processes
        output_minus = pool.starmap(methods[args.method], 
                # add a new random seed for each process. Otherwise will use the
                # same seed for every process
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    window_generator(array, args.wsize, circ = args.circular)), args.seed))
        pool.close()
        pool.join()

        total_plus = len(output_plus)
        output = output_plus + output_minus
        pvals = [val[-1] for val in output]
        _, qvals = fdrcorrection(pvals, method = "i")
        strand = "+"
        for i, (val, qval) in enumerate(zip(output, qvals)):
            if i >= total_plus:
                strand = "-"
            if qval < args.alpha:
                sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\t%s\n"%(chrm, val[0], val[1], val[2],strand, "\t".join([str(entry) for entry in val[3:]]), qval))
        plus_strand.close()
        minus_strand.close()


def region_main(args):
    import bed_utils
    
    # read in bw files
    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
    minus_strand = bwtools.pyBigWig.open(args.infile_minus)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1)

    # function factory for different methods
    methods = {'bootstrap': multiprocess_bootstrap,
            'Poisson': multiprocess_poisson,
            'Larson': multiprocess_larson}

    headers = {'bootstrap': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\texp_count\tsamp_gr_eq\tpvalue\tqvalue\n",
               'Poisson': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\texp_count\tpvalue\tqvalue\n",
               'Larson': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\tmean\tstdev\tpvalue\tqvalue\n"}

    # function factory for different filters
    filters = {'zero': lambda window: window.zero_filter(args.zero_cutoff),
               'avg_cov': lambda window: window.avg_coverage_filter(args.avg_cov_cutoff),
               'max': lambda window: window.max_filter()}

    # combine the filters into a filter list
    if len(args.filters) > 0:
        all_filters = [filters[val] for val in args.filters]
    else:
        all_filters = [lambda window: True]
    
    sys.stdout.write(headers[args.method])
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.bedfile) 
    pool = mp.Pool(args.p)
    output_plus = pool.starmap(methods[args.method], 
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    bed_regions(arrays_plus, inbed.filter_entries(lambda entry: entry["strand"] == "+"))), args.seed))
    pool.close()
    pool.join()

    pool = mp.Pool(args.p)
    output_minus = pool.starmap(methods[args.method], 
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    bed_regions(arrays_minus, inbed.filter_entries(lambda entry: entry['strand'] == "-"))), args.seed))
    pool.close()
    pool.join()

    output = [loc for gene in output_plus for loc in gene] + [loc for gene in output_minus for loc in gene]
    # pvalue is always the last thing in the output of the tests
    pvals = [val[-1] for val in output]
    _, qvals = fdrcorrection(pvals, method = "i")
    for entry, qval in zip(output, qvals):
        if qval < args.alpha:
            sys.stdout.write("%s\t%s\n"%("\t".join([str(val) for val in entry]), qval))

    plus_strand.close()
    minus_strand.close()


 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Python script to call pauses in NET-seq data")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    parser_genome = subparsers.add_parser("Genomewide", help = "Call pauses across the entire genome")
    parser_genome.add_argument('infile_plus', type=str, help = "bigwig file containing raw counts on the plus strand")
    parser_genome.add_argument('infile_minus', type=str, help = "bigwig file containing raw counts on the minus strand")
    parser_genome.add_argument('wsize', type = int, help = "Half the window size. Adds to each side of the query bp. I.e. 100 would give you a \
            total of a 201 bp window with the query bp in the center. Default = 100", default = 100)
    parser_genome.add_argument('--p', type = int, help = "Number of processors to use. Default = 1", default = 1)
    parser_genome.add_argument('--method', type = str, help = "Must be one of Larson (gaussian), Poisson, or bootstrap. Default = Poisson",
            default = "Poisson")
    parser_genome.add_argument('--n_boot', type = int, help = "Number of bootstraps to do for bootstrap method. Default = 10000", default = 10000)
    parser_genome.add_argument('--seed', type = int, help = "Random number generator starting seed for reproducibility. Default = 42", default = 42)
    parser_genome.add_argument('--filters', type = str, nargs = "+",
            help = "List of filters to apply before considering windows. Possibilities include: 'zero', 'max', 'avg_cov', 'cov'. Default = None",
            default = [])
    parser_genome.add_argument('--avg_cov_cutoff', type = float, help = "Cutoff for average coverage filter. Default = 10",
            default = 10.0)
    parser_genome.add_argument('--cov_cutoff', type = float, help = "Cutoff for center coverage filter. Default = 1",
            default = 1.0)
    parser_genome.add_argument('--zero_cutoff', type = float, help = "Cutoff for fraction of zeros in a window. Default = .10",
            default = 0.1)
    parser_genome.add_argument('--circular', type = bool, help = "Consider the genome as circular? Wraps edge windows around end of genome. Default = True",
            default = True)
    parser_genome.add_argument('--alpha', type = float, help = "q-value cutoff for significance. Default = 0.05", default = 0.05)
    parser_genome.set_defaults(func=genomewide_main)


    parser_region = subparsers.add_parser("Region", help = "Call pauses across based on defined regions")
    parser_region.add_argument('infile_plus', type=str, help = "bigwig file containing raw counts on the plus strand")
    parser_region.add_argument('infile_minus', type=str, help = "bigwig file containing raw counts on the minus strand")
    parser_region.add_argument('bedfile', type=str, help = "bed file containing regions to consider")
    parser_region.add_argument('--p', type = int, help = "Number of processors to use. Default = 1", default = 1)
    parser_region.add_argument('--method', type = str, help = "Must be one of Larson (gaussian), Poisson, or bootstrap. Default = Poisson",
            default = "Poisson")
    parser_region.add_argument('--n_boot', type = int, help = "Number of bootstraps to do for bootstrap method. Default = 10000", default = 10000)
    parser_region.add_argument('--seed', type = int, help = "Random number generator starting seed for reproducibility. Default = 42", default = 42)
    parser_region.add_argument('--filters', type = str, nargs = "+",
            help = "List of filters to apply before considering windows. Possibilities include: 'zero', 'avg_cov'. Default = None",
            default = [])
    parser_region.add_argument('--avg_cov_cutoff', type = float, help = "Cutoff for average coverage filter. Default = 20",
            default = 20)
    parser_region.add_argument('--zero_cutoff', type = float, help = "Cutoff for fraction of zeros in a window. Default = .10",
            default = 0.1)
    parser_region.add_argument('--alpha', type = float, help = "q-value cutoff for significance. Default = 0.05", default = 0.05)
    parser_region.set_defaults(func=region_main)

    args = parser.parse_args()
    args.func(args)
