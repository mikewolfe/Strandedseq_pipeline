import bwtools
import numpy as np
import scipy.stats as sci_stats

class GenomicWindow(object):
    def __init__(self, coord, array):
        self.coord = coord
        self.window = array
        self.center_value = self.window[int(len(array)/2)]

    def larson_method(self):
        mean = np.nanmean(self.window)
        stdev = np.nanstd(self.window)
        if self.center_value > mean + 4*stdev:
            peak = True
        else:
            peak = False
        return(self.center_value, mean, stdev,  peak)

    def bootstrap_method(self, n_boot = 10000):
        total = self.window.sum()
        samples = np.random.multinomial(total, [1/len(self.window)]*len(self.window), size = n_boot)
        total_greater_eq = np.sum(np.max(samples, axis=1) >= self.center_value)
        pvalue = (total_greater_eq + 1)/ (n_boot + 1)
        return (self.center_value, total_greater_eq, pvalue)

    def poisson_method(self):
        total = self.window.sum()
        mu = total/len(self.window)
        pvalue = sci_stats.poisson.sf(self.center_value, mu)
        return (self.center_value, mu, total, pvalue)

    def max_filter(self):
        if self.window.max() <= self.center_value:
            return True
        else:
            return False

    def coverage_filter(self, cutoff):
        if self.center_value > cutoff:
            return True
        else:
            return False


def window_generator(array, wsize, circ = True):
    # Augment size of the array by going around the horn
    if circ:
        out_array = np.concatenate((array[-wsize:],array,array[:wsize]))
        # want wsize on either side of the center so need to add one on the left end.
        # This centers the window at center w/ wsize on either side. Total window
        # size should be 2*wsize + 1
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center -wsize, out_array[center-wsize:center+wsize+1])
    else:
        out_array = array
        for center in np.arange(wsize,len(array)-wsize+1, 1):
            yield GenomicWindow(center, out_array[center-wsize:center+wsize+1])

if __name__ == "__main__":
    import sys

    plusfile = sys.argv[1]
    minusfile = sys.argv[2]

    plus_strand = bwtools.pyBigWig.open(plusfile)
    minus_strand = bwtools.pyBigWig.open(minusfile)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1)

    for chrm, array in arrays_plus.items():
        for window in filter(lambda window: window.max_filter(),
                filter(lambda window: window.coverage_filter(20), window_generator(array, 100, circ = True))):
            sys.stdout.write("%s\t%s\t%s\n"%(chrm, window.coord, "\t".join([str(val) for val in window.larson_method()])))
