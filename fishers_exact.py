from scipy import stats
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as pyplot


def fishers_exact(file_string):


    normed = np.loadtxt(file_string, usecols=(1, 2, 3, 4), dtype = 'float', skiprows=1)
    normed = normed.astype(int)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')
    column_sum_es1=np.sum(np.sum(normed[:,0]))
    column_sum_es2=np.sum(np.sum(normed[:,1]))
    column_sum_npc1=np.sum(np.sum(normed[:,2]))
    column_sum_npc2=np.sum(np.sum(normed[:,3]))
    pvalues = np.zeros((np.shape(normed)[0]))


    pvalues_fishers_exact = np.zeros((np.shape(normed)[0], 2))


    for i in range(0,(np.shape(normed)[0])):
	oddsratio, pvalue = stats.fisher_exact([[np.sum(normed[i,0:2]),np.sum(normed[i,2:4])],[((column_sum_es1+column_sum_es2)-np.sum(normed[i,0:2])),((column_sum_npc1+column_sum_npc2)-np.sum(normed[i,2:4]))]])
        #print oddsratio, pvalue, "Fisher's Exact pvalue and Odd's ratio"
        pvalues[i] = pvalue
        pvalues_fishers_exact[i,0] = -10*np.log2(pvalue)
        pvalues_fishers_exact[i,1] = np.log2(np.divide(((np.mean(normed[i,2:4]))+1),((np.mean(normed[i,0:2]))+1)))


    # Make the volcano plots to look at pvalues with respect to fold change

    pyplot.figure()
    pyplot.scatter(pvalues_fishers_exact[:,1],pvalues_fishers_exact[:,0])
    pyplot.scatter(pvalues_fishers_exact[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05)),1],pvalues_fishers_exact[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05)),0],color='red')
    pyplot.title('Volcano plot: ES vs. NPC - pvalues from Fishers Exact Test Code')
    pyplot.ylabel('-10*log2(pvalue)')
    pyplot.xlabel('log2 fold change')
    pyplot.show()

    return pvalues
