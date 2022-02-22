import bwtools
import numpy as np
import scipy.stats as sci_stats
import multiprocessing as mp

class GenomicWindow(object):
    def __init__(self, coord, array):
        self.coord = coord
        self.window = array
        self.center_value = self.window[int(len(array)/2)]

    def larson_method(self):
        mean = np.nanmean(self.window)
        stdev = np.nanstd(self.window)
        pvalue = sci_stats.norm.sf(self.center_value, mean, stdev)
        return(self.coord, self.coord + 1, self.center_value, mean, stdev,  pvalue)

    def bootstrap_method(self, n_boot = 10000, rng = np.random.default_rng()):
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
        total = self.window.sum()
        mu = total/len(self.window)
        pvalue = sci_stats.poisson.sf(self.center_value, mu)
        return (self.coord, self.coord+1, self.center_value, mu, total, pvalue)

    def max_filter(self):
        return self.window.max() <= self.center_value

    def coverage_filter(self, cutoff):
        return self.center_value > cutoff

    def avg_coverage_filter(self, cutoff):
        return np.nanmean(self.window) > cutoff

    def zero_filter(self, cutoff):
        return np.sum(self.window == 0)/len(self.window) < cutoff

def window_generator(array, wsize, circ = True, offset = 0):
    # Augment size of the array by going around the horn
    if circ:
        out_array = np.concatenate((array[-wsize:],array,array[:wsize]))
        # want wsize on either side of the center so need to add one on the left end.
        # This centers the window at center w/ wsize on either side. Total window
        # size should be 2*wsize + 1
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center -wsize + offset, out_array[center-wsize:center+wsize+1])
    else:
        out_array = array
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center + offset, out_array[center-wsize:center+wsize+1])

class RegionWindow(object):
    def __init__(self, bed_entry):
        self.chrm = bed_entry["chrm"]
        self.start = bed_entry["start"]
        self.end = bed_entry["end"]
        self.strand = bed_entry["strand"]
        self.name = bed_entry["name"]

    def __iter__(self):
        for index in np.arange(0, len(self.array)):
            if self.strand == "-":
                rel_coord = self.end - index
            else:
                rel_coord = index
            abs_coord = self.start + index
            yield((index, abs_coord, rel_coord))

    def add_array(self, arrays):
        self.array = arrays[self.chrm][self.start:self.end]

    def larson_method(self):
        mean = np.nanmean(self.array)
        stdev = np.nanstd(self.array)
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            pvalue = sci_stats.norm.sf(this_value, mean, stdev)
            yield (self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, pvalue)

    def poisson_method(self): 
        total = self.array.sum()
        mu = total/len(self.array) 
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            pvalue = sci_stats.poisson.sf(this_value, mu)
            yield (self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, pvalue)
 
    def bootstrap_method(self, n_boot = 10000, rng = np.random.default_rng()):
        # get the total number of counts
        total = self.array.sum()
        # redistribute equally using a multinomial for n_boot times
        samples = rng.multinomial(total, [1/len(self.array)]*len(self.array), size = n_boot)
        for index, abs_coord, rel_coord in self:
            this_value = self.array[index]
            total_greater_eq = np.sum(np.max(samples, axis=1) >= this_value) 
            pvalue = (total_greater_eq + 1)/ (n_boot + 1)
            yield (self.chrm, abs_coord, abs_coord + 1, self.name, this_value, self.strand, rel_coord, total_greater_eq, pvalue)

    def avg_coverage_filter(self, cutoff):
        return np.nanmean(self.array) > cutoff

    def zero_filter(self, cutoff):
        return np.sum(self.array == 0)/len(self.array) < cutoff

def bed_regions(arrays, bed_object):
    for region in bed_object:
        out = RegionWindow(region)
        out.add_array(arrays)
        yield out

def augment_with_random_state(iterable, starting_seed):
    return [(entry, np.random.default_rng(starting_seed + index)) for index, entry in enumerate(iterable)]

def region_multiprocess_helper(region, rng):
    out = []
    for val in region.poisson_method():
        out.append(val)
    return(out)

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

    headers = {'bootstrap': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tsamp_gr_eq\tpvalue\n",
               'Poisson': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tpvalue\n",
               'Larson': "#chrm\tstart\tend\tname\tcount\tstrand\tmean_count\tstd_count\tpvalue\n"}

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
    for chrm, array in arrays_plus.items():
        pool = mp.Pool(args.p)
        # allow the windows to be distributed over the processes
        output = pool.starmap(methods[args.method], 
                # add a new random seed for each process. Otherwise will use the
                # same seed for every process
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    window_generator(array[0:100000], args.wsize, circ = args.circular)), args.seed))
        pool.close()
        pool.join()
        for val in output:
            sys.stdout.write("%s\t%s\t%s\t.\t%s\t+\t%s\n"%(chrm, val[0], val[1], val[2], "\t".join([str(entry) for entry in val[3:]])))
        plus_strand.close()
        minus_strand.close()


#def region_main(args):
#    
#    # read in bw files
#    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
#    minus_strand = bwtools.pyBigWig.open(args.infile_minus)
#
#    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1)
#    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1)
#
#    # function factory for different methods
#    methods = {'bootstrap': lambda(window, rng): window.bootstrap_method(rng = rng),
#            'Poisson': lambda(window, rng): window.poisson_method(),
#            'Larson': lambda(window, rng): window.larson_method()}
#    headers = {'bootstrap': "chrm\tstart\tend\t.\tcount\tstrand\texp_count\ttotal_count\tsamp_gr_eq\tpvalue\n",
#               'Poisson': "chrm\tstart\tend\t.\tcount\tstrand\texp_count\ttotal_count\tpvalue\n",
#               'Larson': "chrm\tstart\tend\t.\tcount\tstrand\tmean_count\tstd_count\tpvalue\n"}
#
#    # function factory for different filters
#    filters = {'zero': lambda window: window.zero_filter(args.zero_cutoff),
#               'avg_cov': lambda window: window.avg_coverage_filter(args.avg_cov_cutoff),
#               'max': lambda window: window.max_filter()}
#    # combine the filters into a filter list
#    all_filters = [filters[val] for val in args.filters]
#    
#    sys.stdout.write(headers[args.method])
#    for chrm, array in arrays_plus.items():
#        pool = mp.Pool(args.p)
#        # allow the windows to be distributed over the processes
#        output = pool.starmap(methods[args.method], 
#                # add a new random seed for each process. Otherwise will use the
#                # same seed for every process
#                augment_with_random_state(
#                    # super filter goes here. Checks each window to make sure it
#                    # passes all filters
#                    filter(lambda window: all([f(window) for f in all_filters]),
#                    # This generates a window for each possible location in the array
#                    window_generator(array[0:100000], args.wsize, circ = args.circular)), args.seed))
#        pool.close()
#        pool.join()
#        for val in output:
#            sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\n"%(chrm, val[0], val[1], val[2], "\t".join([str(entry) for entry in val[3:]]))
 

if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser("Python script to call pauses in NET-seq data")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    parser_genome = subparsers.add_parser("Genomewide", help = "Call pauses across the entire genome")
    parser_genome.add_argument('infile_plus', type=str, help = "bigwig file containing raw counts on the plus strand")
    parser_genome.add_argument('infile_minus', type=str, help = "bigwig file containing raw counts on the plus strand")
    parser_genome.add_argument('outfile', type = str, help = "output file prefix to write, extension will be tsv.gz format")
    parser_genome.add_argument('wsize', type = int, help = "Half the window size. Adds to each side of the query bp. I.e. 100 would give you a \
            total of a 201 bp window with the query bp in the center. Default = 100", default = 100)
    parser_genome.add_argument('--p', type = int, help = "Number of processors to use. Default = 1", default = 1)
    parser_genome.add_argument('--method', type = str, help = "Must be one of Larson (gaussian), Poisson, or bootstrap. Default = Poisson",
            default = "Poisson")
    parser_genome.add_argument('--n_boot', type = int, help = "Number of bootstraps to do for bootstrap method. Default = 10000", default = 10000)
    parser_genome.add_argument('--seed', type = int, help = "Random number generator starting seed for reproducibility. Default = 42", default = 42)
    parser_genome.add_argument('--filters', type = str, nargs = "+",
            help = "List of filters to apply before considering windows. Possibilities include: 'zero', 'max', 'avg_cov'. Default = 'avg_cov' 'max'",
            default = [])
    parser_genome.add_argument('--avg_cov_cutoff', type = float, help = "Cutoff for average coverage filter. Default = 20",
            default = 20)
    parser_genome.add_argument('--zero_cutoff', type = float, help = "Cutoff for fraction of zeros in a window. Default = .10",
            default = 0.1)
    parser_genome.add_argument('--circular', type = bool, help = "Consider the genome as circular? Wraps edge windows around end of genome. Default = True",
            default = True)
    parser_genome.set_defaults(func=genomewide_main)

    args = parser.parse_args()
    args.func(args)

    



#    plusfile = sys.argv[1]
#    minusfile = sys.argv[2]
#    if len(sys.argv) > 3:
#        bedfile = sys.argv[3]
#    else:
#        bedfile = None
#
#    plus_strand = bwtools.pyBigWig.open(plusfile)
#    minus_strand = bwtools.pyBigWig.open(minusfile)
#
#    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1)
#    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1)
#
#    if bedfile:
#        import bed_utils
#        inbed = bed_utils.BedFile()
#        inbed.from_bed_file(bedfile) 
#        pool = mp.Pool(3)
#        output = pool.starmap(region_multiprocess_helper, 
#                    augment_with_random_state(
#                        filter(lambda region: region.avg_coverage_filter(10),
#                        filter(lambda region: region.zero_filter(0.10), 
#                        bed_regions(arrays_plus, inbed))), 42))
#        pool.close()
#        pool.join()
#        for region in output:
#            for val in region:
#                sys.stdout.write("%s\n"%("\t".join([str(entry) for entry in val])))
#    else:
#        for chrm, array in arrays_plus.items():
#            pool = mp.Pool(3)
#            output = pool.starmap(genome_multiprocess_helper, 
#                    augment_with_random_state(
#                        filter(lambda window: window.max_filter(),
#                        filter(lambda window: window.coverage_filter(20), 
#                        window_generator(array[0:100000], 100, circ = False))), 42))
#            pool.close()
#            pool.join()
#            for val in output:
#                sys.stdout.write("%s\t%s\t%s\n"%(chrm, "+", "\t".join([str(entry) for entry in val])))
