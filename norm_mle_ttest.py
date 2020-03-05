from scipy import stats
import numpy as np
from scipy.stats import norm
from scipy.stats import ttest_ind, ttest_ind_from_stats
import matplotlib.pyplot as pyplot


def norm_mle_ttest_volcano(file_string):


    normed = np.loadtxt(file_string, usecols=(1, 2, 3, 4), dtype = 'float', skiprows=1)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')

    pvalues_direct_ttest = np.zeros((np.shape(normed)[0], 2))
    pvalues_mle_params_ttest = np.zeros((np.shape(normed)[0], 2))
    pvalues_unbiasedvar_ttest = np.zeros((np.shape(normed)[0], 2))
    mle_stats = np.zeros((np.shape(normed)[0], 10))
    pvalues = np.zeros((np.shape(normed)[0]))

    for i in range(0,(np.shape(normed)[0])):
    	mle_norm_es_mean, mle_norm_es_var = norm.fit(normed[i,0:2])
    	mle_norm_npc_mean, mle_norm_npc_var = norm.fit(normed[i,2:4])
        mle_stats[i,0] = mle_norm_es_mean  #1st column is ES mean
        mle_stats[i,1] = np.sqrt(mle_norm_es_var)  #2nd column is standard deviation of ES
        mle_stats[i,2] = np.var(normed[i,0:2], axis=0, ddof=1)  #3rd column is sample variance of ES
        mle_stats[i,3] = np.size(normed[i,0:2])   #4th column is ES sample size
        mle_stats[i,4] = np.size(normed[i,0:2]) - 1  #5th column is ES sample size - 1
        mle_stats[i,5] = mle_norm_npc_mean  #6th column is NPC mean
        mle_stats[i,6] = np.sqrt(mle_norm_npc_var)  #7th column is NPC standard deviation
        mle_stats[i,7] = np.var(normed[i,2:4], axis=0, ddof=1)  #8th column is NPC sample variance
        mle_stats[i,8] = np.size(normed[i,2:4])  #9th column is NPC sample size
        mle_stats[i,9] = np.size(normed[i,2:4]) - 1  #10th column is NPC sample size - 1

	# compute pvalues using 2-sample t test with unequal variance (scipy.stats.ttest_ind) computes all parameters directly
	t, p = ttest_ind(normed[i,0:2], normed[i,2:4], equal_var=False)
        pvalues[i] = p
	#print("ttest_ind:            t = %g  p = %g" % (t, p))
        pvalues_direct_ttest[i,0] = -10*np.log2(p)
        pvalues_direct_ttest[i,1] = np.log2(np.divide((mle_stats[i,5]+1),(mle_stats[i,0]+1)))



    # Make the volcano plots to look at pvalues with respect to fold change

    pyplot.figure()
    pyplot.scatter(pvalues_direct_ttest[:,1],pvalues_direct_ttest[:,0])
    pyplot.scatter(pvalues_direct_ttest[np.where(pvalues_direct_ttest[:,0]>=-10*np.log2(0.05)),1],pvalues_direct_ttest[np.where(pvalues_direct_ttest[:,0]>=-10*np.log2(0.05)),0],color='red')
    pyplot.title('Volcano plot: ES vs. NPC - pvalues from direct t-test code')
    pyplot.ylabel('-10*log2(pvalue)')
    pyplot.xlabel('log2 fold change')
    pyplot.show()



    return pvalues
