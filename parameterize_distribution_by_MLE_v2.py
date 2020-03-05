import numpy as np
import scipy.stats
from scipy.optimize import minimize

def parameterize_distribution_by_MLE(data,dist_name,starting_parameters,**kwargs):
    """
     Fits a discrete distribution to data using maximum likelihood estimation.

    Parameters
    ----------
    data : np.array[int]
        A 1D or 2D array of integers to fit a distribution to.
    dist_name : str
        The name of the distribution you want to fit.
    starting_parameters: List[int]
        A list of starting values to pass through `scipy.optimize.minimize()`
    **kwargs : Optional[Dict[str, Any]]
        Kwargs to be passed to `scipy.optimize.minimize()`

    Returns
    -------
    np.array, float
        A list of paramaters that are the MLE estimates of a given distribution, 
		and a scalar corresponding to the maximum log-likelihood.
    """

    #generate lambda function of the negative log-likelihood that can be passed to `scipy.optimize.minimize()`
    #as an objective function
    if dist_name == 'poisson':
        neg_loglikelihood_fn = lambda mu : -np.sum([scipy.stats.poisson.logpmf(x,mu) for x in data])
    elif dist_name == 'nbinom':
        neg_loglikelihood_fn = lambda (r,p) : -np.sum([scipy.stats.nbinom.logpmf(x,r,p) for x in data])
    else:
        raise Warning("The distribution name is not recognized. Please set dist_name to 'poisson' or 'nbinom'")

    #Perform MLE by minimizing the negative log-likelihood with respect to 'mu' or with respect to 'r' and 'p'
    mle_parameters = minimize(neg_loglikelihood_fn, starting_parameters, **kwargs).x
	
    max_loglikelihood = -1*neg_loglikelihood_fn(mle_parameters)

    return mle_parameters, max_loglikelihood
