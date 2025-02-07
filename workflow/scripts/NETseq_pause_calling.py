import bwtools
import numpy as np
import scipy.stats as sci_stats
from statsmodels.stats.multitest import fdrcorrection
import multiprocessing as mp
import sys

def create_final_filter(array, wsize, circular, filter_keep, filters):
    """ Create a filter array the size of the genome of places you want to keep
    
    Args:
        array - array under consideration
        wsize - half the wsize to consider around each location
        circular - consider this array as circular?
        filter_keep - an initial filter to start out with. Must be same length
                      as the input array. Typically np.ones(array, dtype = bool)
        filter - list of functions to apply to each window

    Returns:
        A new filter array after application of each filter function to each
        window in a nested manner. I.e. once a window filter returns false that
        window continues to stay false in the output filter array
    """

    for f in filters:
        for loc, arr in generate_windows(array, wsize, filter_keep, circular):
            filter_keep[loc] = f(arr)
    return filter_keep
            
def max_filter(array, wsize):
    """ Filter that returns true when the center value of a window is the max
    value """
    return array[wsize] >= np.max(array)

def local_max_filter(array, wsize):
    """ Filter that returns true when the center value of a window is greater
    than the values next to it"""
    return (array[wsize] >= array[wsize-1]) & (array[wsize] >= array[wsize+1])

def coverage_filter(array, wsize, cutoff):
    """ Filter that returns true when the center value of a window is greater
    than some cutoff"""
    return array[wsize] >= cutoff

def avg_coverage_filter(array, wsize, cutoff):
    """ Filter that returns true when the average value of a window is greater
    than some cutoff"""
    return np.mean(array) >= cutoff

def zero_filter(array, wsize, cutoff):
    """ Filter that returns true when the total number of zeros is less than
    some cutoff"""
    return (np.sum(array == 0)/len(array)) < cutoff


def generate_windows(array, wsize, filter_keep, circular = False):
    """ A generator for created windows of a specified size across an array

    Options to consider arrays that are circular
    
    Args:
        array - array under consideration
        wsize - half the wsize to consider around each location
        filter_keep - Locations to consider. Must be same length as the input
                      array. For all loctations this can be 
                      np.ones(array, dtype = bool)
        circular - consider this array as circular?

    Yields:
        a tuple with two items:
            1. the location of the center of the window
            2. an array of size wsize*2 + 1 centered at each location in the
               array. 
    """
    array_length = len(array)
    for loc in np.arange(0, array_length)[filter_keep]:
        start = loc - wsize
        end = loc + wsize + 1
        if start < 0:
            if circular:
                yield (loc, np.concatenate((array[start:], array[:end])))
            else:
                continue
        elif end >= array_length:
            if circular:
                yield (loc, np.concatenate((array[start:], array[:(end - array_length)])))
            else:
                continue
        else:
            yield (loc, array[start:end])

def generate_windows_nofilter(array, wsize, circular = False):
    """ A generator for created windows of a specified size across an array

    Options to consider arrays that are circular
    
    Args:
        array - array under consideration
        wsize - half the wsize to consider around each location
        circular - consider this array as circular?

    Yields:
        an array of size wsize*2 + 1 centered at each location in the array. 
    """
    array_length = len(array)
    for loc in np.arange(0, array_length):
        start = loc - wsize
        end = loc + wsize + 1
        if start < 0:
            if circular:
                yield (loc, np.concatenate((array[start:], array[:end])))
            else:
                continue
        elif end >= array_length:
            if circular:
                yield (loc, np.concatenate((array[start:], array[:(end - array_length)])))
            else:
                continue
        else:
            yield array[start:end]



def larson_method(window, wsize, *extra, no_center = False):
    """ Calculate rarity of seeing center value 
    
    Expects that the values in the window are normally distributed
    
    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        rng (unused) - a np.random.generator

    Returns:
        A list of information about the test including:
        coord, coord + 1,
        center_value, mean of window, stdev of window, pvalue = (1-cdf
        normal(mean, stdev))
    """
    loc, array = window
    mean = np.mean(array)
    stdev = np.std(array)

    median = np.median(array)
    num_zeros = np.count_nonzero(array==0)

    center_value = array[wsize]
    pvalue = sci_stats.norm.sf(center_value, mean, stdev)
    Zscore = (center_value - mean)/stdev
    return [loc, loc+1, center_value, mean, stdev, Zscore, median, num_zeros, pvalue]


def larson_nocenter_method(window, wsize, *extra):
    """ Calculate rarity of seeing center value 
    
    Expects that the values in the window are normally distributed
    
    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        rng (unused) - a np.random.generator

    Returns:
        A list of information about the test including:
        coord, coord + 1,
        center_value, mean of window, stdev of window, pvalue = (1-cdf
        normal(mean, stdev))
    """

    loc, array = window
    bkg_array = np.concatenate((array[:wsize], array[wsize+1:])) 
    mean = np.mean(bkg_array)
    stdev = np.std(bkg_array)
    median = np.median(bkg_array)
    num_zeros = np.count_nonzero(bkg_array==0)

    center_value = array[wsize]
    pvalue = sci_stats.norm.sf(center_value, mean, stdev)
    Zscore = (center_value - mean)/stdev
    return [loc, loc+1, center_value, mean, stdev, Zscore, median, num_zeros, pvalue]


def larson_robust_method(window, wsize, *extra):
    """ Calculate rarity of seeing center value 
    
    Expects that the values in the window are normally distributed
    
    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        rng (unused) - a np.random.generator

    Returns:
        A list of information about the test including:
        coord, coord + 1,
        center_value, mean of window, stdev of window, pvalue = (1-cdf
        normal(mean, stdev))
    """

    loc, array = window
    bkg_array = np.concatenate((array[:wsize], array[wsize+1:])) 
    median = np.median(bkg_array)
    mad = sci_stats.median_abs_deviation(bkg_array)
    num_zeros = np.count_nonzero(bkg_array==0)

    center_value = array[wsize]
    pvalue = sci_stats.norm.sf(center_value, median, mad*1.4826)
    Zscore = (center_value - median)/(mad*1.4826)
    return [loc, loc+1, center_value, median, mad, Zscore, num_zeros,pvalue]

def estimate_nbinom(array):
    mu = np.mean(array)
    std = np.std(array)
    p = mu/(std**2)
    n = (mu**2)/(std**2 - mu)
    return n,p

def churchman_method(window, wsize, *extra):
    """ Calculate rarity of seeing center value 
    
    Expects that the values in the window are normally distributed
    
    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        rng (unused) - a np.random.generator

    Returns:
        A list of information about the test including:
        coord, coord + 1,
        center_value, mean of window, stdev of window, pvalue = (1-cdf
        normal(mean, stdev))
    """
    loc, array = window
    bkg_array = np.concatenate((array[:wsize], array[wsize+1:])) 
    n, p = estimate_nbinom(bkg_array)
    mean = sci_stats.nbinom.mean(n, p)
    std = sci_stats.nbinom.std(n, p)
    center_value = array[wsize]
    pvalue = sci_stats.nbinom.sf(center_value, n, p)
    Zscore = (center_value - mean)/std
    return [loc, loc+1, center_value, mean, std, Zscore, pvalue]


def poisson_method(window, wsize, *extra):
    """ Calculate rarity of seeing center value

    Assumes that the reads are distributed according to a poisson
    distribution under the null. Lamba is determined as the average read
    count in the window

    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        rng (unused) - a np.random.generator

    Returns:
        A list of information about the test including:

        coord, coord + 1, center_value, expected count, total count in
        window, pvalue = 1 - cdf poisson(avg_read_count)

    """
    loc, array = window
    total = int(array.sum())
    mu = total/len(array)
    center_value = int(array[wsize])
    pvalue = sci_stats.poisson.sf(center_value, mu)
    return [loc, loc+1, center_value, mu, total, pvalue]


def bootstrap_method(window, wsize, rng, n_boot = 10000):
    """ Calculate rarity of seeing center value 

    Assumes that the total count of reads is distributed according to a 
    multinomial under the null.

    Args:
        window - a 1D numpy array
        wsize - half the size of the window. The center location
        n_boot (int): number of random samples from a multinomial to draw
        rng (np.random.Generator): state needed for random number generation

    Returns:
        A list of information about the test including:

        coord, coord + 1, center_value, expected count, total count in
        window, total random samples where max in sample >= actual
        center_value, pvalue = (total_greater_eq + 1) / (n_boot + 1)
    
    """
    loc, array = window
    # get the total number of counts
    total = int(np.sum(array))
    arr_size = len(array)
    center_value = array[wsize]
    # redistribute equally using a multinomial for n_boot times
    samples = rng.multinomial(total, [1/arr_size]*arr_size, size = n_boot)
    # check the number of samples where the max is higher than the actual
    total_greater_eq = np.sum(np.max(samples, axis=1) >= center_value)
    #total_greater_eq = np.sum(np.array([np.max(row) for row in samples]) >= self.center_value)
    # compute the p value. Use a prior of one sample being greater than 
    # actual
    pvalue = (total_greater_eq + 1)/ (n_boot + 1)
    return [loc, loc+1, center_value, total/arr_size, total, total_greater_eq, pvalue]

def gini_coefficient(data, unbiased = False):
    '''
    Calculate the gini coefficient for a 1D set of data
    Arguments:
        data - 1D numpy array of data
        unbiased - boolean to determine if unbiased estimator should be used
    References:
        https://en.wikipedia.org/wiki/Gini_coefficient
        https://www.statsdirect.com/help/nonparametric_methods/gini_coefficient.htm
        https://stackoverflow.com/questions/48999542/more-efficient-weighted-gini-coefficient-in-python
    Formula:
        Standard:
        $G = \frac{1}{n} (n + 1 - 2 \frac{\sum_{i=1}^n (n + 1 - i)y_i}{\sum_{i=1}^n y_i})$
        Unbiased:                                                                
        $G = \frac{1}{n-1} (n + 1 - 2 \frac{\sum_{i=1}^n (n + 1 - i)y_i}{\sum_{i=1}^n y_i})$
    '''
    data_sorted = np.sort(data)
    n = len(data)
    # the sum of the cumulative sum vector gives the numerator when the data
    # is sorted in ascending order. The float here prevents integer overflows
    cumulative = np.cumsum(data_sorted, dtype = float)
    # the last value in a cumulative vector is the total sum
    stat = (n + 1 - 2 * np.sum(cumulative)/ cumulative[-1])
    if unbiased:
        stat = stat/(n-1)
    else:
        stat = stat/n
    return stat

def windowed_gini(data, wsize=100):
    '''
    Calculated rolling window of ginis for each position in array
    '''
    out = []
    for window in generate_windows_nofilter(data, wsize):
        out.append(gini_coefficient(window))
    out = np.array(out)
    return out


def pause_per_kb(data, wsize=100, Zscore_cutoff = 4):
    '''
    Calculate pauses per kb of a region
    Uses Larson_nocenter method to call pauses

    '''
    total_pauses = 0
    for loc, array in generate_windows(data, wsize, np.ones(len(data), dtype=bool)):
        bkg_array = np.concatenate((array[:wsize], array[wsize+1:])) 
        mean = np.mean(bkg_array)
        stdev = np.std(bkg_array)
        center_value = array[wsize]
        Zscore = (center_value - mean)/stdev
        if Zscore > Zscore_cutoff:
            total_pauses += 1
    stat = (total_pauses/(len(data)/1000))

    return stat

def pause_per_kb_multi(data, Zscore_cutoff = 4):
    '''
    Calculate pauses per kb of a region
    Estimates pause calling from multinomial

    '''
    total_pauses = 0
    total = np.sum(data)
    length = len(data)
    expected = total/length
    sigma = np.sqrt(expected*(1-(1/length)))
    stat = np.sum((data - expected)/sigma > Zscore_cutoff)
    stat = (stat/length)*1000
    return stat

def gini_boot_method(window, wsize, rng, n_boot = 10000):
    loc, array = window
    stat = gini_coefficient(array)
    boot = sci_stats.bootstrap((array,), gini_coefficient, method = 'percentile',
            vectorized = False, n_resamples = n_boot)
    num_zeros = np.count_nonzero(array==0)
    total = np.nansum(array)
    return [loc, loc+1, stat, boot.standard_error, boot.confidence_interval[0], boot.confidence_interval[1], num_zeros, total]

def stat_bootstrapped(array, true_val, stat_func, rng, n_boot=1000, alpha = 0.05):
    total = int(np.sum(array))
    arr_size = len(array)
    sample_data = rng.multinomial(total, [1/arr_size]*arr_size, size = n_boot)
    samples = np.apply_along_axis(stat_func, 1, sample_data)
    total_greater_eq = np.sum(samples >= true_val)
    median_boot = np.median(samples)
    low, high  = np.quantile(samples, [alpha, 1-alpha])
    mad_boot = sci_stats.median_abs_deviation(samples)
    pvalue = (total_greater_eq +1)/(n_boot + 1)
    return [median_boot, mad_boot, low, high, pvalue]


def stat_bootstrapped_compare(array, true_val, stat_func, compare_func, rng, n_boot=1000, alpha = 0.05):
    total = int(np.sum(array))
    arr_size = len(array)
    sample_data = rng.multinomial(total, [1/arr_size]*arr_size, size = n_boot)
    samples = np.apply_along_axis(stat_func, 1, sample_data)
    total_greater_eq = np.sum(samples >= true_val)
    mean_boot = np.mean(samples)
    std_boot = np.std(samples)
    Z_boot = (true_val - mean_boot)/std_boot
    median_boot = np.median(samples)
    low, high  = np.quantile(samples, [alpha, 1-alpha])
    mad_boot = sci_stats.median_abs_deviation(samples)
    pvalue = (total_greater_eq +1)/(n_boot + 1)
    compared_samples = compare_func(true_val, samples)
    mean_compare = np.mean(compared_samples)
    std_compare = np.std(compared_samples)
    median_compare = np.median(compared_samples)
    low_compare, high_compare = np.quantile(compared_samples, [alpha, 1-alpha])
    mad_compare = sci_stats.median_abs_deviation(compared_samples)
    return [mean_boot, std_boot, Z_boot,median_boot, mad_boot, low, high, pvalue, mean_compare, std_compare, median_compare, mad_compare, low_compare, high_compare]


def downsample_array(array, total_count, rng, replace = False):
    # Create a list of indices based on the counts in the array
    indices = np.repeat(np.arange(len(array)), array)

    # Sample without replacement from the indices
    sampled_indices = rng.choice(indices, total_count, replace=replace)
 
    # Create the downsampled array using np.bincount for efficiency
    downsampled_array = np.bincount(sampled_indices, minlength=len(array))

    return downsampled_array

def multiarray_compare_stat_bootstrapped(arrays, stat_func, baseline_sample, total, rng, n_boot = 1000, alpha = 0.05, pseudocount = 1, include_null_columns = True):
    arr_size = arrays[0].size
    sampled_stats = []
    for array in arrays:
        this_array = array + pseudocount
        sampled_array = rng.multinomial(total, this_array/np.sum(this_array), size = n_boot)
        sampled_stats.append(np.apply_along_axis(stat_func, 1, sampled_array))
    if include_null_columns:
        null_sampled_array = rng.multinomial(total, [1/arr_size]*arr_size, size = n_boot)
        null_sample = np.apply_along_axis(stat_func, 1, null_sampled_array)

    output_vector = []

    ref_sample = sampled_stats[baseline_sample]

    for i, sampled_stat in enumerate(sampled_stats):
        # uncertainty on gini
        output_vector.append(np.mean(sampled_stat))
        output_vector.append(np.std(sampled_stat))
        low, high  = np.quantile(sampled_stat, [alpha, 1-alpha])
        output_vector.append(low)
        output_vector.append(high)
        # uncertainty relative to reference sample
        ratio = sampled_stat/ref_sample
        output_vector.append(np.mean(ratio))
        output_vector.append(np.std(ratio))
        low, high  = np.quantile(ratio, [alpha, 1-alpha])
        output_vector.append(low)
        output_vector.append(high)
        # uncertainty relative to null
        if include_null_columns:
            output_vector.append(np.mean(sampled_stat/null_sample))
            output_vector.append(np.std(sampled_stat/null_sample))
            low, high  = np.quantile(sampled_stat/null_sample, [alpha, 1-alpha])
            output_vector.append(low)
            output_vector.append(high)
        # array total and zero
        output_vector.append(np.sum(arrays[i]))
        output_vector.append(np.count_nonzero((arrays[i]) == 0))
    output_vector.append(total)
    return output_vector

def multiarray_dwnsample_compare_stat_bootstrapped(arrays, stat_func, baseline_sample, total, rng, n_boot = 1000, alpha = 0.05, pseudocount = 1, include_null_columns = True, multi_sample = True, bootstrap_est = True):
    arr_size = arrays[0].size
    sampled_stats = []
    for array in arrays: 
        sampled_array = np.zeros((n_boot, len(array)), dtype="float32")
        for bootstrap in np.arange(0, n_boot):
            this_array = downsample_array(array, total, rng)
            this_array += pseudocount
            if bootstrap_est:
                if multi_sample:
                    sampled_array[bootstrap,:] = rng.multinomial(total, this_array/np.sum(this_array), size = 1)
                else:
                    sampled_array[bootstrap,:] = downsample_array(this_array, total, rng, replace = True)
            else:
                sampled_array[bootstrap,:] = this_array
        sampled_stats.append(np.apply_along_axis(stat_func, 1, sampled_array))

    if include_null_columns:
        null_sampled_array = rng.multinomial(total, [1/arr_size]*arr_size, size = n_boot)
        null_sample = np.apply_along_axis(stat_func, 1, null_sampled_array)
    
    output_vector = []

    ref_sample = sampled_stats[baseline_sample]

    for i, sampled_stat in enumerate(sampled_stats):
        # uncertainty on gini
        output_vector.append(np.mean(sampled_stat))
        output_vector.append(np.std(sampled_stat))
        low, high  = np.quantile(sampled_stat, [alpha, 1-alpha])
        output_vector.append(low)
        output_vector.append(high)
        # uncertainty relative to reference sample
        ratio = sampled_stat/ref_sample
        #ratio = ratio[np.isfinite(ratio)]
        output_vector.append(np.mean(ratio))
        output_vector.append(np.std(ratio))
        low, high  = np.quantile(ratio, [alpha, 1-alpha])
        output_vector.append(low)
        output_vector.append(high)
        # uncertainty relative to null
        if include_null_columns:
            output_vector.append(np.mean(sampled_stat/null_sample))
            output_vector.append(np.std(sampled_stat/null_sample))
            low, high  = np.nanquantile(sampled_stat/null_sample, [alpha, 1-alpha])
            output_vector.append(low)
            output_vector.append(high)
        # array total and zero
        output_vector.append(np.sum(arrays[i]))
        output_vector.append(np.count_nonzero((arrays[i]) == 0))
    output_vector.append(total)
    return output_vector


def multiarray_stat(arrays, stat_func, baseline_sample, total_count, rng, dwnsample = False, include_extras = True):
    '''
    Calculate stats for multiple arrays for one region
    '''
    # all arrays are the same size so get the size of the first array
    arr_size = arrays[0].size
    # List to hold the output
    stats = []
    # calculate a stat for each array
    for i, array in enumerate(arrays): 
        if dwnsample:
            this_array = downsample_array(array, total_count, rng)
        else:
            this_array = array
        this_total = np.sum(this_array)
        stats.append(stat_func(this_array))

    output_vector = []
    
    # add some extra stats to the stat vector
    for i, stat in enumerate(stats):
        # stat
        output_vector.append(stat)
        if include_extras:
            # array total and zero
            output_vector.append(np.sum(arrays[i]))
            output_vector.append(np.count_nonzero((arrays[i]) == 0))
    if include_extras:
        output_vector.append(this_total)
    return output_vector

def multiarray_compare_header(sample_names, ref_sample = 0, include_null_columns = True):
    out = []
    out.append("chrm")
    out.append("start")
    out.append("end")
    out.append("strand")
    out.append("name")
    for sample in sample_names:
        out.append(sample + "_mean")
        out.append(sample + "_std")
        out.append(sample + "_low")
        out.append(sample + "_hi")
        out.append(sample + "_ovr_" + sample_names[ref_sample] + "_mean")
        out.append(sample + "_ovr_" + sample_names[ref_sample] + "_std")
        out.append(sample + "_ovr_" + sample_names[ref_sample] + "_low")
        out.append(sample + "_ovr_" + sample_names[ref_sample] + "_hi")
        if include_null_columns:
            out.append(sample + "_ovr_null" + "_mean")
            out.append(sample + "_ovr_null" + "_std")
            out.append(sample + "_ovr_null" + "_low")
            out.append(sample + "_ovr_null" + "_hi")
        out.append(sample + "_total")
        out.append(sample + "_zeros")
    out.append("sampled_count")
    header = "\t".join(out)
    return header + "\n"

def multiarray_stat_header(sample_names, ref_sample = 0, include_extras = True):
    out = []
    out.append("chrm")
    out.append("start")
    out.append("end")
    out.append("strand")
    out.append("name")
    for sample in sample_names:
        out.append(sample + "_stat")
        if include_extras:
            out.append(sample + "_total")
            out.append(sample + "_zeros")

    if include_extras:
        out.append("sampled_count")
    header = "\t".join(out)
    return header + "\n"

def gini_null_method(window, wsize, rng, n_boot = 10000):
    loc, array = window
    stat = gini_coefficient(array)
    median_boot, mad_boot, low, high, pvalue = stat_bootstrapped(array, stat, gini_coefficient, rng, n_boot)
    total = np.nansum(array)
    num_zeros = np.count_nonzero(array==0)
    return [loc, loc+1, stat, median_boot, mad_boot, low, high, pvalue, num_zeros, total]

def gini_null_ratio_method(window, wsize, rng, n_boot = 10000):
    loc, array = window
    stat = gini_coefficient(array)
    mean_boot, std_boot, Z_boot, median_boot, mad_boot, low, high, pvalue, mean_compare, std_compare, median_compare, mad_compare, low_compare, high_compare = stat_bootstrapped_compare(array, stat, gini_coefficient, np.divide, rng, n_boot)
    total = np.nansum(array)
    num_zeros = np.count_nonzero(array==0)
    return [loc, loc+1, stat, median_boot, mad_boot, low, high, pvalue, median_compare, mad_compare, low_compare, high_compare, num_zeros, total]


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


    def gini_null_ratio_method(self, n_boot = 10000, rng = np.random.default_rng()):
        out = []
        array = self.array
        stat = gini_coefficient(array)
        mean_boot, std_boot, Z_boot, median_boot, mad_boot, low, high, pvalue, mean_compare, std_compare, median_compare, mad_compare, low_compare, high_compare = stat_bootstrapped_compare(array, stat, gini_coefficient, np.divide, rng, n_boot)
        total = np.nansum(array)
        num_zeros = np.count_nonzero(array==0)
        out.append([self.chrm, self.start, self.end, self.name, stat, self.strand, mean_boot, std_boot, Z_boot, median_boot, mad_boot, low, high, mean_compare, std_compare, median_compare, mad_compare, low_compare, high_compare, num_zeros, total, pvalue])
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


class RegionMultiWindow(object):
    '''
    Store information about a region in the genome for multiple samples

    Attributes:
        chrm (str): chromsome region is on
        start (int): absolute start coordinate
        end (int): absolute end coordinate
        strand (str) : strand of region
        name (str): name of region
        array (np.array): array of values for the region
    '''
    def __init__(self, bed_entry, arrays = None, pseudocount = 1, sample_cov = None, Zscorecutoff = 4, wsize=100, start_padding = 0, end_padding = 0):
        self.chrm = bed_entry["chrm"]
        self.start = bed_entry["start"] - start_padding
        self.end = bed_entry["end"] + end_padding
        self.strand = bed_entry["strand"]
        self.name = bed_entry["name"]
        self.pseudocount = pseudocount
        self.sample_cov = sample_cov
        self.Zscorecutoff = Zscorecutoff
        self.wsize = wsize
        if arrays is not None:
            self.add_arrays(arrays)
            self.get_min_count()

    def add_arrays(self, arrays):
        """
        Adds an array to the attributes of the region

        Args:
            array (np.array): np array of values typically a view
        """
        self.arrays = [array[self.chrm][self.start:self.end].astype(int)  for array in arrays]

    def get_min_count(self):
        counts = [np.nansum(array) for array in self.arrays]
        self.min_count = np.min(counts)

    def GiniCompare(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_compare_stat_bootstrapped(self.arrays, gini_coefficient, 0, total_count, rng, n_boot, alpha, self.pseudocount)

    def PausePerKb(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_compare_stat_bootstrapped(self.arrays, pause_per_kb, 0, total_count, rng, n_boot, alpha, self.pseudocount)

    def PausePerKb_multi(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_compare_stat_bootstrapped(self.arrays, pause_per_kb_multi, 0, total_count, rng, n_boot, alpha, self.pseudocount, include_null_columns = False)

    def PausePerKb_multi_dwnsample(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_dwnsample_compare_stat_bootstrapped(self.arrays, pause_per_kb_multi, 0, total_count, rng, n_boot, alpha, self.pseudocount, include_null_columns = False)

    def PausePerKb_dwn_dwnsample(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_dwnsample_compare_stat_bootstrapped(self.arrays, pause_per_kb_multi, 0, total_count, rng, n_boot, alpha, self.pseudocount, include_null_columns = False, multi_sample = False)

    def PausePerKb_dwnsample(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        total_count = self.min_count
        # make the assumption that the first array is the reference sample always
        return multiarray_dwnsample_compare_stat_bootstrapped(self.arrays, pause_per_kb_multi, 0, total_count, rng, n_boot, alpha, self.pseudocount, include_null_columns = False, bootstrap_est = False)


    def PausePerKb_simple(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        def stat_func(x):
            return pause_per_kb_multi(x, self.Zscorecutoff)
        # make the assumption that the first array is the reference sample always
        return multiarray_stat(self.arrays, pause_per_kb_multi, 0, total_count, rng,dwnsample = False)

    def PausePerKb_dwnsample_simple(self, n_boot = 10000, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        def stat_func(x):
            return pause_per_kb_multi(x, self.Zscorecutoff)
        # make the assumption that the first array is the reference sample always
        return multiarray_stat(self.arrays, stat_func, 0, total_count, rng,dwnsample = True)

    def WindowedGini(self, n_boot = 1, rng = np.random.default_rng(), alpha = 0.05):
        if self.sample_cov is not None:
            total_count = self.sample_cov * (self.end - self.start)
        else:
            total_count = self.min_count
        def stat_func(x):
            return windowed_gini(x, self.wsize)

        return multiarray_stat(self.arrays, stat_func, 0, total_count, rng,dwnsample = True, include_extras = False)





def multiarray_bed_regions(arrays_plus, arrays_minus, bed_object, pseudocount = 1, sample_cov = None, Zscorecutoff = 4, wsize = 100, upstream = 0, downstream = 0):
    for region in bed_object:
        if region["strand"] == "-":
            start_padding = downstream
            end_padding = upstream
            out = RegionMultiWindow(region, arrays_minus, pseudocount = pseudocount, sample_cov = sample_cov, Zscorecutoff = Zscorecutoff, wsize = wsize, end_padding = end_padding)
        else:
            start_padding = upstream
            end_padding = downstream
            out = RegionMultiWindow(region, arrays_plus, pseudocount = pseudocount, sample_cov = sample_cov, Zscorecutoff = Zscorecutoff, wsize = wsize, start_padding = start_padding)
        yield out



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

def augment_with_random_state(iterable, args):
    """ Zip an iterator with a np random seed

    This is needed to ensure that each process does not inherit the same
    randomstate during multiprocessing

    Args:
        iterable: anything you need to iterate
        starting_seed (int): starting seed for the rng

    Returns
        a list of (rng, iterable)
    """
    return [(entry, np.random.default_rng(args.seed + index), args.n_boot) for index, entry in enumerate(iterable)]

def augment_with_random_state_and_wsize(iterable, args):
    """ Zip an iterator with a np random seed

    This is needed to ensure that each process does not inherit the same
    randomstate during multiprocessing

    Args:
        iterable: anything you need to iterate
        starting_seed (int): starting seed for the rng

    Returns
        a list of (rng, iterable)
    """
    return [(entry, args.wsize, np.random.default_rng(args.seed + index), args.n_boot) for index, entry in enumerate(iterable)]

def multiprocess_bootstrap(window, rng, bootstraps):
    return window.bootstrap_method(rng=rng,n_boot = bootstraps)

def multiprocess_poisson(window, rng, *extras):
    return window.poisson_method()

def multiprocess_larson(window, rng, *extras):
    return window.larson_method()

def multiprocess_gini_null_ratio(window, rng, bootstraps):
    return window.gini_null_ratio_method(rng = rng, n_boot = bootstraps)

def multiprocess_multigini(region, rng, bootstraps):
    return region.GiniCompare(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb(region, rng, bootstraps):
    return region.PausePerKb(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_multi(region, rng, bootstraps):
    return region.PausePerKb_multi(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_multi_dwnsample(region, rng, bootstraps):
    return region.PausePerKb_multi_dwnsample(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_dwn_dwnsample(region, rng, bootstraps):
    return region.PausePerKb_dwn_dwnsample(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_dwnsample(region, rng, bootstraps):
    return region.PausePerKb_dwnsample(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_simple(region, rng, bootstraps):
    return region.PausePerKb_simple(rng = rng, n_boot = bootstraps)

def multiprocess_multiPausePerKb_dwnsample_simple(region, rng, bootstraps):
    return region.PausePerKb_dwnsample_simple(rng = rng, n_boot = bootstraps)

def multiprocess_windowedgini_simple(region, rng, bootstraps):
    return region.WindowedGini(rng = rng, n_boot = bootstraps)


def run_tests_over_strand(arrays, methods, filters, initial_filters, args):
    all_output = []
    for chrm, array in arrays.items():
        filter_keep = initial_filters[chrm]
        filter_keep = create_final_filter(array, args.wsize, args.circular, filter_keep, filters)
        pool = mp.Pool(args.p)
        # allow the windows to be distributed over the processes
        output = pool.starmap(methods[args.method], 
                # add a new random seed for each process. Otherwise will use the
                # same seed for every process
                augment_with_random_state_and_wsize(
                    # This generates a window for each possible location in the array
                    generate_windows(array, args.wsize, filter_keep, circular = args.circular), args))
        pool.close()
        pool.join()
        output = [[chrm] + entry for entry in output]
        all_output += output
    return all_output

def genomewide_gini(args):

    # read in bw files
    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
    minus_strand = bwtools.pyBigWig.open(args.infile_minus)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1, nan_to_zero = True)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1, nan_to_zero = True)
    
    # initialize filters for which locations to consider
    starting_filters= {'+' : {}, '-' : {}}

    # add a way to set problem regions to zero
    if args.zero_regions is not None:
        import bed_utils
        zero_bed = bed_utils.BedFile()
        zero_bed.from_bed_file(args.zero_regions)
        for entry in zero_bed:
            strand = entry['strand']
            chrm = entry['chrm']
            start = entry['start']
            end = entry['end']
            if strand == "-":
                arrays_minus[chrm][start:end] = 0
            else:
                arrays_plus[chrm][start:end] = 0

    if args.regions is not None:
        for chrm, array in arrays_plus.items():
            starting_filters['+'][chrm] = np.zeros(len(array), dtype = bool)
        for chrm, array in arrays_minus.items():
            starting_filters['-'][chrm] = np.zeros(len(array), dtype = bool)
        import bed_utils
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(args.regions) 
        for entry in inbed:
            starting_filters[entry['strand']][entry['chrm']][entry['start']:entry['end']] = True
    else:
        for chrm, array in arrays_plus.items():
            starting_filters['+'][chrm] = np.ones(len(array), dtype = bool)
        for chrm, array in arrays_minus.items():
            starting_filters['-'][chrm] = np.ones(len(array), dtype = bool)

    # function factory for different methods
    methods = {'Gini_boot': gini_boot_method,
               'Gini_null': gini_null_method,
               'Gini_null_ratio': gini_null_ratio_method}

    # function factory for different filters
    filters = {'zero': lambda arr: zero_filter(arr, args.wsize, args.zero_cutoff),
               'avg_cov': lambda arr: avg_coverage_filter(arr, args.wsize, args.avg_cov_cutoff),
               'cov': lambda arr: coverage_filter(arr, args.wsize, args.cov_cutoff),
               'max': lambda arr: max_filter(arr, args.wsize),
               'local_max': lambda arr: local_max_filter(arr, args.wsize)}

    # combine the filters into a filter list
    if len(args.filters) > 0:
        all_filters = [filters[val] for val in args.filters]
    else:
        all_filters = [lambda arr: True]
    
    output_plus = run_tests_over_strand(arrays_plus, methods, all_filters, starting_filters['+'], args)
    output_minus = run_tests_over_strand(arrays_minus, methods, all_filters, starting_filters['-'], args)

    total_plus = len(output_plus)
    output = output_plus + output_minus
    headers = {'Gini_boot': "#chrm\tstart\tend\tname\tGini\tstrand\tboot_stderr\tboot_low\tboot_hi\tnumzeros\ttotal_count\n",
               'Gini_null': "#chrm\tstart\tend\tname\tGini\tstrand\tnull_median\tnull_mad\tnull_low\tnull_high\tpvalue\tnumzeros\ttotal_count\tqvalue\n",
               'Gini_null_ratio': "#chrm\tstart\tend\tname\tGini\tstrand\tnull_median\tnull_mad\tnull_low\tnull_high\tpvalue\tratio_median\tratio_mad\tratio_low\tratio_high\tnumzeros\ttotal_count\tqvalue\n"}
    sys.stdout.write(headers[args.method])
    if args.method in ['Gini_null', 'Gini_null_ratio']:
        # filter on pvalues when possible
        pvals = [out[8] for out in output]
        if args.no_qvalue: 
            qval = np.nan
            strand = "+"
            for i, (val, pval) in enumerate(zip(output, pvals)):
                if i >= total_plus:
                    strand = "-"
                if pval < args.alpha:
                    sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\t%s\n"%(val[0], val[1], val[2], val[3], strand, "\t".join([str(entry) for entry in val[4:]]), qval))
        else:
            _, qvals = fdrcorrection(pvals, method = "i")
            strand = "+"
            for i, (val, qval) in enumerate(zip(output, qvals)):
                if i >= total_plus:
                    strand = "-"
                if qval < args.alpha:
                    sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\t%s\n"%(val[0], val[1], val[2], val[3],strand, "\t".join([str(entry) for entry in val[4:]]), qval))
    else:
        # The bootstrap method doesn't compare to a null thus doesn't have a p-value
        for i, val in enumerate(output):
            if i >= total_plus:
                strand = "-"
            else:
                strand = "+"
            sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\n"%(val[0], val[1], val[2], val[3], strand, "\t".join([str(entry) for entry in val[4:]])))
    plus_strand.close()
    minus_strand.close()


def genomewide_fast(args):
    
    # read in bw files
    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
    minus_strand = bwtools.pyBigWig.open(args.infile_minus)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1, nan_to_zero = True)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1, nan_to_zero = True)
    
    # initialize filters for which locations to consider
    starting_filters= {'+' : {}, '-' : {}}


    # add a way to set problem regions to zero

    if args.zero_regions is not None:
        import bed_utils
        zero_bed = bed_utils.BedFile()
        zero_bed.from_bed_file(args.zero_regions)
        for entry in zero_bed:
            strand = entry['strand']
            chrm = entry['chrm']
            start = entry['start']
            end = entry['end']
            if strand == "-":
                arrays_minus[chrm][start:end] = 0
            else:
                arrays_plus[chrm][start:end] = 0

    if args.regions is not None:
        for chrm, array in arrays_plus.items():
            starting_filters['+'][chrm] = np.zeros(len(array), dtype = bool)
        for chrm, array in arrays_minus.items():
            starting_filters['-'][chrm] = np.zeros(len(array), dtype = bool)
        import bed_utils
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(args.regions) 
        for entry in inbed:
            starting_filters[entry['strand']][entry['chrm']][entry['start']:entry['end']] = True
    else:
        for chrm, array in arrays_plus.items():
            starting_filters['+'][chrm] = np.ones(len(array), dtype = bool)
        for chrm, array in arrays_minus.items():
            starting_filters['-'][chrm] = np.ones(len(array), dtype = bool)

        


    # function factory for different methods
    methods = {'bootstrap': bootstrap_method,
            'Poisson': poisson_method,
            'Larson': larson_method,
            'Larson_nocenter': larson_nocenter_method,
            'Larson_robust': larson_robust_method,
            'Churchman': churchman_method}

    headers = {
            'bootstrap': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tsamp_gr_eq\tpvalue\tqvalue\n",
               'Poisson': "#chrm\tstart\tend\tname\tcount\tstrand\texp_count\ttotal_count\tpvalue\tqvalue\n",
               'Larson': "#chrm\tstart\tend\tname\tcount\tstrand\tmean_count\tstd_count\tZscore\tmedian_count\tnum_zeros\tpvalue\tqvalue\n",
               'Larson_nocenter': "#chrm\tstart\tend\tname\tcount\tstrand\tmean_count\tstd_count\tZscore\tmedian_count\tnum_zeros\tpvalue\tqvalue\n",
               'Larson_robust': "#chrm\tstart\tend\tname\tcount\tstrand\tmedian_count\tmad_count\tZscore\tnum_zeros\tpvalue\tqvalue\n",
               'Churchman': "#chrm\tstart\tend\tname\tcount\tstrand\tmean_count\tstd_count\tZscore\tpvalue\tqvalue\n"}

    # function factory for different filters
    filters = {'zero': lambda arr: zero_filter(arr, args.wsize, args.zero_cutoff),
               'avg_cov': lambda arr: avg_coverage_filter(arr, args.wsize, args.avg_cov_cutoff),
               'cov': lambda arr: coverage_filter(arr, args.wsize, args.cov_cutoff),
               'max': lambda arr: max_filter(arr, args.wsize),
               'local_max': lambda arr: local_max_filter(arr, args.wsize)}

    # combine the filters into a filter list
    if len(args.filters) > 0:
        all_filters = [filters[val] for val in args.filters]
    else:
        all_filters = [lambda arr: True]
    
    sys.stdout.write(headers[args.method])
    output_plus = run_tests_over_strand(arrays_plus, methods, all_filters, starting_filters['+'], args)
    output_minus = run_tests_over_strand(arrays_minus, methods, all_filters, starting_filters['-'], args)

    total_plus = len(output_plus)
    output = output_plus + output_minus
    pvals = [val[-1] for val in output]
    if args.no_qvalue: 
        qval = np.nan
        strand = "+"
        for i, (val, pval) in enumerate(zip(output, pvals)):
            if i >= total_plus:
                strand = "-"
            if pval < args.alpha:
                sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\t%s\n"%(val[0], val[1], val[2], val[3],strand, "\t".join([str(entry) for entry in val[4:]]), qval))
    else:
        _, qvals = fdrcorrection(pvals, method = "i")
        strand = "+"
        for i, (val, qval) in enumerate(zip(output, qvals)):
            if i >= total_plus:
                strand = "-"
            if qval < args.alpha:
                sys.stdout.write("%s\t%s\t%s\t.\t%s\t%s\t%s\t%s\n"%(val[0], val[1], val[2], val[3],strand, "\t".join([str(entry) for entry in val[4:]]), qval))
    plus_strand.close()
    minus_strand.close()


def region_main(args):
    import bed_utils
    
    # read in bw files
    plus_strand = bwtools.pyBigWig.open(args.infile_plus)
    minus_strand = bwtools.pyBigWig.open(args.infile_minus)

    arrays_plus = bwtools.bigwig_to_arrays(plus_strand, res = 1, nan_to_zero = True)
    arrays_minus = bwtools.bigwig_to_arrays(minus_strand, res = 1, nan_to_zero = True)

    # function factory for different methods
    methods = {'bootstrap': multiprocess_bootstrap,
            'Poisson': multiprocess_poisson,
            'Larson': multiprocess_larson,
            'Gini_null_ratio': multiprocess_gini_null_ratio}

    headers = {'bootstrap': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\texp_count\tsamp_gr_eq\tpvalue\tqvalue\n",
               'Poisson': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\texp_count\tpvalue\tqvalue\n",
               'Larson': "#chrm\tstart\tend\tname\tcount\tstrand\trel_coord\tmean\tstdev\tpvalue\tqvalue\n",
               'Gini_null_ratio': "#chrm\tstart\tend\tname\tgini\tstrand\tmean_boot\tstd_boot\tZ_boot\tmedian_boot\tmad_boot\tlow\thigh\tmean_compare\tstd_compare\tmedian_compare\tmad_compare\tlow_compare\thigh_compare\tnum_zeros\ttotal\tpvalue\tqvalue\n"}

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
                    bed_regions(arrays_plus, inbed.filter_entries(lambda entry: entry["strand"] == "+"))), args))
    pool.close()
    pool.join()

    pool = mp.Pool(args.p)
    output_minus = pool.starmap(methods[args.method], 
                augment_with_random_state(
                    # super filter goes here. Checks each window to make sure it
                    # passes all filters
                    filter(lambda window: all([f(window) for f in all_filters]),
                    # This generates a window for each possible location in the array
                    bed_regions(arrays_minus, inbed.filter_entries(lambda entry: entry['strand'] == "-"))), args))
    pool.close()
    pool.join()

    output = [loc for gene in output_plus for loc in gene] + [loc for gene in output_minus for loc in gene]
    # pvalue is always the last thing in the output of the tests
    pvals = [val[-1] for val in output]
    if args.no_qvalue:
        qval = np.nan
        for entry, pval in zip(output, pvals):
            if pval < args.alpha:
                sys.stdout.write("%s\t%s\n"%("\t".join([str(val) for val in entry]), qval))
    else:
        _, qvals = fdrcorrection(pvals, method = "i")
        for entry, qval in zip(output, qvals):
            if qval < args.alpha:
                sys.stdout.write("%s\t%s\n"%("\t".join([str(val) for val in entry]), qval))

    plus_strand.close()
    minus_strand.close()


def Ginicompare_main(args):
    import bed_utils
    
    # read in bw files
    plus_strand_handles = bwtools.open_multiple_bigwigs(args.infiles_plus)
    arrays_plus = bwtools.convert_bigwigs_to_arrays(plus_strand_handles, res = 1)

    minus_strand_handles = bwtools.open_multiple_bigwigs(args.infiles_minus)
    arrays_minus = bwtools.convert_bigwigs_to_arrays(minus_strand_handles, res = 1)

    # read in regions
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.bedfile) 
    if args.metric == "PausePerKb_multi" or args.metric == "PausePerKb_multi_dwnsample" or args.metric == "PausePerKb_dwn_dwnsample" or args.metric == "PausePerKb_dwnsample":
        header = multiarray_compare_header(args.sample_names, include_null_columns = False)
    elif args.metric == "PausePerKb_simple" or args.metric == "PausePerKb_dwnsample_simple":
        header = multiarray_stat_header(args.sample_names)
    elif args.metric == "WindowedGini":
        header = multiarray_stat_header(args.sample_names, include_extras = False)
    else:
        header = multiarray_compare_header(args.sample_names, include_null_columns = True)
 
    function_factory = {"Gini": multiprocess_multigini,
                        "PausePerKb": multiprocess_multiPausePerKb,
                        "PausePerKb_multi": multiprocess_multiPausePerKb_multi,
                        "PausePerKb_multi_dwnsample": multiprocess_multiPausePerKb_multi_dwnsample,
                        "PausePerKb_dwn_dwnsample": multiprocess_multiPausePerKb_dwn_dwnsample,
                        "PausePerKb_dwnsample": multiprocess_multiPausePerKb_dwnsample,
                        "PausePerKb_simple": multiprocess_multiPausePerKb_simple,
                        "PausePerKb_dwnsample_simple": multiprocess_multiPausePerKb_dwnsample_simple,
                        "WindowedGini": multiprocess_windowedgini_simple}
 
    pool = mp.Pool(args.p)
    output = pool.starmap(function_factory[args.metric], 
                augment_with_random_state(
                    multiarray_bed_regions(arrays_plus, arrays_minus, inbed, pseudocount = args.pseudocount, sample_cov = args.samplecov, Zscorecutoff = args.Zscorecutoff, wsize = args.wsize, upstream = args.wsize, downstream = args.wsize), args))
    pool.close()
    pool.join()
    sys.stdout.write(header)
    for region, array in zip(inbed, output):
        if len(array[0]) == 1:
            region_info = "\t".join([region["chrm"], str(region["start"]), str(region["end"]), region["strand"], region["name"]])
            output_info = "\t".join(["%s"%(value) for value in array])
            sys.stdout.write(region_info + "\t" + output_info + "\n")
        elif len(array[0]) > 1:
            for pos in range(0, len(array[0])):
                region_info = "\t".join([region["chrm"], str(region["start"] + pos), str(region["start"] + pos + 1), region["strand"], region["name"]])
                output_info = "\t".join(["%s"%(value[pos]) for value in array])
                sys.stdout.write(region_info + "\t" + output_info + "\n")
        else:
            raise ValueError("Output array dimension doesn't make sense %s"%(array.ndim))

    bwtools.close_multiple_bigwigs(plus_strand_handles)
    bwtools.close_multiple_bigwigs(minus_strand_handles)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Python script to call pauses in NET-seq data")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    parser_gini = subparsers.add_parser("Genomewide_gini", help = "Test for uneveness in read distribution")
    parser_gini.add_argument('infile_plus', type=str, help = "bigwig file containing raw counts on the plus strand")
    parser_gini.add_argument('infile_minus', type=str, help = "bigwig file containing raw counts on the minus strand")
    parser_gini.add_argument('wsize', type = int, help = "Half the window size. Adds to each side of the query bp. I.e. 100 would give you a \
            total of a 201 bp window with the query bp in the center. Default = 100", default = 100)
    parser_gini.add_argument('--p', type = int, help = "Number of processors to use. Default = 1", default = 1)
    parser_gini.add_argument('--method', type = str, help = "Must be one of Gini_boot (bootstrap the Gini estimate) or Gini_null (bootstrap a null multinomial distro). Default = Gini_null",
            default = "Gini_null")
    parser_gini.add_argument('--n_boot', type = int, help = "Number of bootstraps to do for bootstrap method. Default = 10000", default = 10000)
    parser_gini.add_argument('--seed', type = int, help = "Random number generator starting seed for reproducibility. Default = 42", default = 42)
    parser_gini.add_argument('--filters', type = str, nargs = "+",
            help = "List of filters to apply before considering windows. Possibilities include: 'zero'- max fraction of zeros in window, 'max',\
                    'avg_cov'-window must be above an avg coverage, 'cov'-center of window must be above cutoff, 'local_max'- coverage must be higher than neighbors. Default = None",
            default = [])
    parser_gini.add_argument('--avg_cov_cutoff', type = float, help = "Cutoff for average coverage filter. Default = 1",
            default = 1.0)
    parser_gini.add_argument('--cov_cutoff', type = float, help = "Cutoff for center coverage filter. Default = 1",
            default = 1.0)
    parser_gini.add_argument('--zero_cutoff', type = float, help = "Cutoff for fraction of zeros in a window. Default = .10",
            default = 0.1)
    parser_gini.add_argument('--regions', type = str, help = "Regions to only consider. Ignore anything outside these regions. Default = None",
            default = None)
    parser_gini.add_argument('--zero_regions', type = str, help = "Regions to set to zero. Useful for removing known artifacts. Default = None",
            default = None)
    parser_gini.add_argument('--resolution', type = int, help = "Roll the windows at what resolution? Default = 1 bp ", default = 1)
    parser_gini.add_argument('--circular', type = bool, help = "Consider the genome as circular? Wraps edge windows around end of genome. Default = True",
            default = True)
    parser_gini.add_argument('--alpha', type = float, help = "q-value cutoff for significance. Default = 0.05", default = 0.05)
    parser_gini.add_argument('--no_qvalue', action = "store_true", help = "don't do q-value calculations. Use p-value for sig cutoff")
    parser_gini.set_defaults(func=genomewide_gini)

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
            help = "List of filters to apply before considering windows. Possibilities include: 'zero'- max fraction of zeros in window, 'max',\
                    'avg_cov'-window must be above an avg coverage, 'cov'-center of window must be above cutoff, 'local_max'- coverage must be higher than neighbors. Default = None",
            default = [])
    parser_genome.add_argument('--avg_cov_cutoff', type = float, help = "Cutoff for average coverage filter. Default = 1",
            default = 1.0)
    parser_genome.add_argument('--cov_cutoff', type = float, help = "Cutoff for center coverage filter. Default = 1",
            default = 1.0)
    parser_genome.add_argument('--zero_cutoff', type = float, help = "Cutoff for fraction of zeros in a window. Default = .10",
            default = 0.1)
    parser_genome.add_argument('--regions', type = str, help = "Regions to only consider. Ignore anything outside these regions. Default = None",
            default = None)
    parser_genome.add_argument('--zero_regions', type = str, help = "Regions to set to zero. Useful for removing known artifacts. Default = None",
            default = None)
    parser_genome.add_argument('--circular', type = bool, help = "Consider the genome as circular? Wraps edge windows around end of genome. Default = True",
            default = True)
    parser_genome.add_argument('--alpha', type = float, help = "q-value cutoff for significance. Default = 0.05", default = 0.05)
    parser_genome.add_argument('--no_qvalue', action = "store_true", help = "don't do q-value calculations. Use p-value for sig cutoff")
    parser_genome.set_defaults(func=genomewide_fast)


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
    parser_region.add_argument('--no_qvalue', action = "store_true", help = "don't do q-value calculations. Use p-value for sig cutoff")
    parser_region.set_defaults(func=region_main)

    parser_ginicompare = subparsers.add_parser("GiniCompare", help = "Compare pausing across locations")
    parser_ginicompare.add_argument('bedfile', type=str, help = "bed file containing regions to consider")
    parser_ginicompare.add_argument('--sample_names', type=str, nargs = "+", help = "names of the samples")
    parser_ginicompare.add_argument('--infiles_plus', type=str, nargs = "+", help = "bigwig files containing raw counts on the plus strand")
    parser_ginicompare.add_argument('--infiles_minus', type=str, nargs = "+", help = "bigwig files containing raw counts on the minus strand")
    parser_ginicompare.add_argument('--p', type = int, help = "Number of processors to use. Default = 1", default = 1)
    parser_ginicompare.add_argument('--n_boot', type = int, help = "Number of bootstraps to do for bootstrap method. Default = 10000", default = 10000)
    parser_ginicompare.add_argument('--seed', type = int, help = "Random number generator starting seed for reproducibility. Default = 42", default = 42)
    parser_ginicompare.add_argument('--pseudocount', type = int, help = "Pseudocount", default = 0)
    parser_ginicompare.add_argument('--samplecov', type = float, help = "How much coverage depth to sample at", default = None)
    parser_ginicompare.add_argument('--metric', type = str, help = "Which metric to use (Gini or PausePerKb) default = Gini", default = "Gini")
    parser_ginicompare.add_argument('--Zscorecutoff', type = float, help = "What Zscore cutoff to use", default = 4)
    parser_ginicompare.add_argument('--wsize', type = int, help = "Window size for windowed gini. Also pads each array by size", default = 0)
    parser_ginicompare.set_defaults(func=Ginicompare_main)

    args = parser.parse_args()
    args.func(args)
