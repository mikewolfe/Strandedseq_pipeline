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
        MAD = stats.median_absolute_deviation(array, nan_policy='omit', scale = 1.4826)
        median = np.nanmedian(array)
    else:
        MAD = stats.median_absolute_deviation(array)
        median = np.median(array)
    return normalize_1D(array, median, MAD)

