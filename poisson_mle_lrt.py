from scipy import stats
import numpy as np
import matplotlib.pyplot as pyplot
from parameterize_distribution_by_MLE_v2 import *
from scipy.stats import chisqprob


def poisson_mle_lrt(file_string):


    normed = np.loadtxt(file_string, usecols=(1, 2, 3, 4), dtype = 'float', skiprows=1)
    normed = normed.astype(int)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')
    pvalues = np.zeros((np.shape(normed)[0]))
    pvalues_poisson_mle_lrt = np.zeros((np.shape(normed)[0], 2))
    mle_stats = np.zeros((np.shape(normed)[0], 3))

    for i in range(0,(np.shape(normed)[0])):

        mle_poisson_es_param, ll_es = parameterize_distribution_by_MLE(normed[i,0:2],dist_name='poisson',
                                                   starting_parameters=[1],method='Nelder-Mead')
        mle_poisson_npc_param, ll_npc = parameterize_distribution_by_MLE(normed[i,2:4],dist_name='poisson',
                                                   starting_parameters=[1],method='Nelder-Mead')
        mle_poisson_null_param, ll_null = parameterize_distribution_by_MLE(normed[i,0:4],dist_name='poisson',
                                                   starting_parameters=[1],method='Nelder-Mead')

        # print normed[i,0:4], mle_poisson_es_param, (np.mean(normed[i,0:2])), mle_poisson_npc_param, np.mean(normed[i,2:4]), mle_poisson_null_param, np.mean(normed[i,0:4]), "Here is crosscheck for MLE fits - they should match sample mean for Poisson"

        mle_stats[i,0] = mle_poisson_es_param
        mle_stats[i,1] = mle_poisson_npc_param
        mle_stats[i,2] = mle_poisson_null_param

        # print ll_es, ll_npc, ll_null
        

        # Now, set up code to compute LRT test statistic 

        LRT = -2*((ll_null)-(ll_es+ll_npc))
        
        #LRT = -2*np.log2((ll_null*ll_null)/(ll_es*ll_npc))
        # print LRT

	pv_lrt = chisqprob(LRT, 1)
        pvalues[i] = pv_lrt

        pvalues_poisson_mle_lrt[i,0] = -10*np.log2(pv_lrt)
        pvalues_poisson_mle_lrt[i,1] = np.log2(np.divide((mle_stats[i,1]+1),(mle_stats[i,0]+1)))


    # Make the volcano plots to look at pvalues with respect to fold change

    pyplot.figure()
    pyplot.scatter(pvalues_poisson_mle_lrt[:,1],pvalues_poisson_mle_lrt[:,0])
    pyplot.scatter(pvalues_poisson_mle_lrt[np.where(pvalues_poisson_mle_lrt[:,0]>=-10*np.log2(0.05)),1],pvalues_poisson_mle_lrt[np.where(pvalues_poisson_mle_lrt[:,0]>=-10*np.log2(0.05)),0],color='red')
    pyplot.title('Volcano plot: ES vs. NPC - pvalues from Poisson log-likelihood ratio test code')
    pyplot.ylabel('-10*log2(pvalue)')
    pyplot.xlabel('log2 fold change')
    pyplot.show()


    return pvalues
