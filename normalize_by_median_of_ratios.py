from __future__ import division
import numpy as np

def normalize_by_median_of_ratios(counts_matrix,logged=True):
    """
    Normalizes by finding the median of ratios of each sample count to the
    geometric mean of each gene across samples. 

    Parameters
    ----------
    counts_matrix : np.ndarray
        The counts matrix of all sequencing libraries. The rows are the counts
        for each gene and the columns represent the different libraries.
    logged : bool
        If true, the data will be not be logged again when calculating the
        geometric mean. 

    Returns
    -------
    np.ndarray
        A 2 x 2 matrix like counts_matrix, but normalized for median of ratios.
    """

    #make sure that all counts are floats
    if counts_matrix.dtype is not np.float32:
        counts_matrix = counts_matrix.astype(np.float32)
    #checks for logging to calculate geometric mean of each gene properly 
    if logged:
        geo_mean = np.nanmean(counts_matrix,axis=1)
    else:
        geo_mean = np.exp(np.nanmean(np.log2(counts_matrix + 1),axis=1)) - 1
    #calculate ratio of each gene count to geometric mean
    ratios = np.zeros_like(counts_matrix)
    for j in range(counts_matrix.shape[1]):
        ratios[:,j] = counts_matrix[:,j] / geo_mean
    #take the median ratio across genes
    size_factors = np.nanmedian(ratios,axis=0)
    size_factors_ratio=[np.divide(size_factors[0],size_factors[0]), np.divide(size_factors[1],size_factors[0]),np.divide(size_factors[2],size_factors[0]),np.divide(size_factors[3],size_factors[0])]
    #divide counts matrix by size factor for each sample
    normalized_matrix = counts_matrix / size_factors_ratio
    return normalized_matrix 
