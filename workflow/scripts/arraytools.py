import numpy as np
from scipy import stats
import logging

def call_bool_peaks(array, consolidate = 0):
    """ Take a logical array and find the start and stop of ones across
    the array. Consolidate any peaks that are within consolidate.

    TODO: Deal with peaks that go over the end of the array

    Args:
        array (np.array) - logical array to call peaks from.
        consolidate (int) - number of bins to consolidate peaks over
    Returns:
        peak_indices (list of lists) - a list of [start, stop] indices in the
                                        array
    """
    # first find all the places where there is a change from 1 to 0 or vice
    # versa, here we pad the incoming array with a zero on each end, then take
    # the difference along the array and finally take the absolute value to flip
    # the places where the difference was 0 - 1
    changes = np.abs(np.diff(np.concatenate(([0], array.view(np.int8), [0]))))
    # changes is now a logical array, I then find all the indices where changes
    # happen and reshape them into an ndarray of start, stop locations
    start_stops = np.where(changes)[0].reshape(-1, 2)
    if start_stops.size == 0:
        logging.warning("No bp was considered to be within a peak.")
        return []

    # now lets consolidate any peaks that are within consolidate
    consolidate_peaks = [[start_stops[0][0], start_stops[0][1]]]
    consolidated_peaks = 0
    removed_peaks = 0
    for start, stop in start_stops[1:]:
        if start - consolidate_peaks[-1][1] < consolidate:
            consolidate_peaks[-1][1] = stop
            consolidated_peaks += 1
        else:
            if stop-start > consolidate:
                consolidate_peaks.append([start, stop])
            else:
                removed_peaks += 1
    logging.info("Consolidated %s peaks within %s bps"%(consolidated_peaks, consolidate))
    logging.info("Removed %s peaks < %s bps"%(removed_peaks, consolidate))
    return consolidate_peaks

def threshold_peaks(array, threshold, consolidate=1, below = False):
    """ Function to call peaks by a threshold and consolidate peaks that are
    nearby

    Args:
        array - 1 dimensional numpy array
        threshold - float, a hard threshold for a peak 
        consolidate - int, number of adjacent locations to consolidate
        below - bool, consider peaks below the threshold (above if false)
    Returns:
        list - start,end tuples of peaks
    """
    if below:
        peaks = array < threshold
    else:
        peaks = array > threshold
    return call_bool_peaks(peaks, consolidate = consolidate)

def normalize_1D(array, center, scale):
    """
    Function to normalize a 1 dimensional array using the following formula:

    (array - center)/ scale

    Args:
        array - 1 dimensional numpy array
        center - a single value for the center (commonly mean(array))
        scale - a single value for the scale (commonly std(array))
    Returns:
        outarray - numpy array normalized
    """

    return (array - center)/scale

def robustz_1D(array, ignore_nan = True):
    """
    Function to Robust Z normalize an array. I.e.

    RZ = array - median(array)/(1.4826* MAD)

    Args:
        array - 1 dimensional numpy array
        ignore_nan - boolean, ignore individual nans and propogate those nans to final array
    Returns:
        outarray - numpy array robustZ normalized
    """
    if ignore_nan:
        # default scale is specified here
        MAD = stats.median_abs_deviation(array, nan_policy='omit', scale = 'normal')
        median = np.nanmedian(array)
    else:
        MAD = stats.median_absolute_deviation(array)
        median = np.median(array)
    return normalize_1D(array, median, MAD)

def smooth_1D(array, wsize, kernel_type = "flat", edge = "mirror"):
    """
    Function to smooth a 1D signal using a convolution

    Args:
        array - 1 dimensional numpy array
        wsize - size of half the window
        kernel - one of flat, gaussian
        edge - one of mirror, wrap
    Returns:
        outarray - smoothed numpy array of same size
    """
    tot_size = wsize*2 + 1
    if tot_size >= len(array):
        raise ValueError("Window must be smaller than array. Window %s, array %s"%(tot_size, len(array)))

    if kernel_type == "flat":
        kernel = np.ones(tot_size)
    elif kernel_type == "gaussian":
        x = np.arange(-wsize, wsize + 1)
        # make the weights span the kernel 
        sigma = (wsize*2)/6
        kernel = np.exp(- (x**2)/(2. *sigma**2))
    else:
        raise ValueError("kernel_type must be one of flat or gaussian. Not %s"%kernel_type)

    if edge == "mirror":
        # pad with a mirror reflection of the ends
        s = np.r_[array[wsize:0:-1], array, array[-2:-wsize-2:-1]]
    elif edge == "wrap":
        # pad with a wrap around the horn
        s = np.r_[array[(-wsize):], array, array[:wsize]]
    else:
        raise ValueError("edge must be one of mirror or wrap. Not %s"%edge)
    kernel = kernel/np.sum(kernel)
    out = np.convolve(s, kernel, mode = 'valid')
    return out

def savgol_1D(array, wsize, polyorder=5, deriv=0, delta = 1.0, edge = "mirror"):
    """
    Function to smooth a 1D signal using a Savitzky-Golay filter

    Args:
        array - 1 dimensional numpy array
        wsize - size of half the window
        edge - one of mirror, wrap, nearest, constant
    Returns:
        outarray - smoothed numpy array of same size

    See also:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
    """
    from scipy import signal
    tot_size = wsize*2 + 1

    if tot_size >= len(array):
        raise ValueError("Window must be smaller than array. Window %s, array %s"%(tot_size, len(array)))
    if polyorder >= tot_size:
        raise ValueError("polyorder shouldn't be larger than window. Window %s, polyorder %s"%(tot_size, polyorder))
    return signal.savgol_filter(array, tot_size, polyorder, deriv, delta, mode = edge)
