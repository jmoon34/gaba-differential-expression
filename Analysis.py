from load_RNAseq_counts import load_RNAseq_counts
from load_gene_info import load_gene_info
from filter_low_count_genes import filter_low_count_genes
from normalize_by_median_of_ratios import normalize_by_median_of_ratios
from write_RNAseq_counts import write_RNAseq_counts
from norm_mle_ttest import *
from fishers_exact import *
from ks_test import *
from poisson_mle_lrt import *
import statsmodels.stats.multitest as mtc


import numpy as np
import matplotlib.pyplot as pyplot


def main():
    #load counts file, add psudocount, log transform and normalize by median-of-ratios
    counts_matrix, gene_order, rep_list = load_RNAseq_counts('slc12a_counts.txt', delimiter=' ')
    new_gene_order, _, filtered_matrix = filter_low_count_genes(gene_order, counts_matrix, count_threshold=0)
    matrix = filtered_matrix + 1  # add pseudocount +1
    matrix = np.log(matrix)  # log transform the counts
    mor_matrix = normalize_by_median_of_ratios(matrix)
    write_RNAseq_counts(new_gene_order, rep_list, mor_matrix,'final_ESNPC_mor_normalized.txt', delimiter=' ')
    mean = np.zeros(np.shape(counts_matrix)[1])
    var = np.zeros(np.shape(counts_matrix)[1])
    std = np.zeros(np.shape(counts_matrix)[1])
    med = np.zeros(np.shape(counts_matrix)[1])
    min = np.zeros(np.shape(counts_matrix)[1])
    max = np.zeros(np.shape(counts_matrix)[1])

    # Here, create a box plot of counts across all four samples for the untransformed data
    data = [counts_matrix[:, 0], counts_matrix[:, 1], counts_matrix[:, 2], counts_matrix[:, 3]]
    pyplot.figure()
    pyplot.boxplot(data)
    pyplot.ylim(0, 70)
    pyplot.xticks([1, 2, 3, 4], ['v65es_rep1', 'v65es_rep2', 'pnpc_rep1', 'pnpc_rep2'])
    pyplot.title('Boxplot for ' + rep_list[0] + rep_list[1] + rep_list[2] + rep_list[3])
    pyplot.show()

    # Plot histogram, relative frequency distribution and boxplot for untransformed data
    for j in range(0, np.shape(counts_matrix)[1]):
        # relative frequency distribution
        pyplot.figure()
        pyplot.hist((counts_matrix[:,j]), 20, normed=True)
        pyplot.title('Relative Frequency Distribution in ' + rep_list[j])
        pyplot.show()
        # relative frequency distribution of normalized matrix
        pyplot.figure()
        pyplot.hist(mor_matrix[:, j], 20, normed=True)
        pyplot.xlabel('Counts')
        pyplot.ylabel('Relative Frequency')
        pyplot.title('Relative Frequency Ddistribution for MOR Normalized Counts in ' + rep_list[j])
        pyplot.show()

    #plots scatterplots for transformed data
    pyplot.figure()
    pyplot.scatter(mor_matrix[:, 0], mor_matrix[:, 1])
    pyplot.title('Scatterplot for Counts in ' + rep_list[0] + ' vs. ' + rep_list[1])
    pyplot.xlabel(rep_list[0])
    pyplot.ylabel(rep_list[1])
    pyplot.show()
    pyplot.figure()
    pyplot.scatter(counts_matrix[:, 2], counts_matrix[:, 3])
    pyplot.title('Scatterplot for Counts in ' + rep_list[2] + ' vs. ' + rep_list[3])
    pyplot.xlabel(rep_list[0])
    pyplot.ylabel(rep_list[1])
    pyplot.show()
    pyplot.figure()
    pyplot.scatter(counts_matrix[:, 0], counts_matrix[:, 2])
    pyplot.title('Scatterplot for Counts in ' + rep_list[0] + ' vs. ' + rep_list[2])
    pyplot.xlabel(rep_list[0])
    pyplot.ylabel(rep_list[1])
    pyplot.show()
    pyplot.figure()
    pyplot.scatter(counts_matrix[:, 1], counts_matrix[:, 3])
    pyplot.title('Scatterplot for Counts in ' + rep_list[1] + ' vs. ' + rep_list[3])
    pyplot.xlabel(rep_list[0])
    pyplot.ylabel(rep_list[1])
    pyplot.show()


    # Summary statistics for the transformed data
    for j in range(0, np.shape(mor_matrix)[1]):
        print np.mean(mor_matrix[:, j]), "here is the mean of counts in sample ", rep_list[j]
        mean[j] = np.mean(mor_matrix[:, j])
        print np.var(mor_matrix[:, j]), "here is the population variance of counts in sample ", rep_list[j]
        var[j] = np.var(mor_matrix[:, j])
        print np.median(mor_matrix[:, j]), "here is the median of counts in sample ", rep_list[j]
        med[j] = np.median(mor_matrix[:, j])
        print np.std(mor_matrix[:, j]), "here is the sample standard deviation of counts in sample ", rep_list[j]
        std[j] = np.std((mor_matrix[:, j]))
        print np.amin(mor_matrix[:, j]), "here is the minimum of counts in sample", rep_list[j]
        min[j] = np.amin(mor_matrix[:, j])
        print np.amax(mor_matrix[:, j]), "here is the maximum of counts in sample", rep_list[j]
        max[j] = np.amax(mor_matrix[:, j])

    p_ttest = norm_mle_ttest_volcano('final_ESNPC_mor_normalized.txt')
    p_KS = ks_test('final_ESNPC_mor_normalized.txt')
    p_fisher = fishers_exact('final_ESNPC_mor_normalized.txt')
    p_LRT = poisson_mle_lrt('final_ESNPC_mor_normalized.txt')

    pyplot.figure()
    pyplot.hist(p_ttest, 20, normed=True)
    pyplot.title('Relative Frequency Distribution for t-test pvalues')
    pyplot.xlim(0, 1)
    pyplot.show()

    pyplot.figure()
    pyplot.hist(p_fisher, 20, normed=True)
    pyplot.title('Relative Frequency Distribution for Fisher\'s exact test pvalues')
    pyplot.xlim(0, 1)
    pyplot.show()

    pyplot.figure()
    pyplot.hist(p_KS, 20, normed=True)
    pyplot.title('Relative Frequency Distribution for K-S test pvalues')
    pyplot.xlim(0, 1)
    pyplot.show()

    pyplot.figure()
    pyplot.hist(p_LRT, 20, normed=True)
    pyplot.title('Relative Frequency Distribution for LRT pvalues')
    pyplot.xlim(0, 1)
    pyplot.show()

    rej_ttest, p_corr_ttest = mtc.multipletests(p_ttest, alpha=0.05, method='fdr_bh')[:2]
    rej_fisher, p_corr_fisher = mtc.multipletests(p_fisher, alpha=0.05, method='fdr_bh')[:2]
    rej_KS, p_corr_KS = mtc.multipletests(p_KS, alpha=0.05, method='fdr_bh')[:2]
    rej_LRT, p_corr_LRT = mtc.multipletests(p_LRT, alpha=0.05, method='fdr_bh')[:2]

    filtered_ttest_corr = np.delete(p_corr_ttest, np.where(p_corr_ttest > 0.05))
    filtered_fisher_corr = np.delete(p_corr_fisher, np.where(p_corr_fisher > 0.05))
    filtered_KS_corr = np.delete(p_corr_KS, np.where(p_corr_KS > 0.05))
    filtered_LRT_corr = np.delete(p_corr_LRT, np.where(p_corr_LRT > 0.05))

    print "Number of significant genes for t-test after correction: ", np.size(filtered_ttest_corr), "out of", np.size(p_corr_ttest)
    print "Number of significant genes for Fisher's exact test after correction: ", np.size(filtered_fisher_corr), "out of", np.size(p_corr_fisher)
    print "Names of genes that have significant p-values from Fisher: ", new_gene_order[np.where(p_corr_fisher < 0.05)]
    print "Number of significant genes for K-S test after correction: ", np.size(filtered_KS_corr), "out of", np.size(p_corr_KS)
    print "Number of significant genes for Poisson LRT after correction: ", np.size(filtered_LRT_corr), "out of", np.size(p_corr_LRT)
    print "Names of genes that have significant p-values from Poisson LRT: ", new_gene_order[np.where(p_corr_LRT < 0.05)]

if __name__ == '__main__':
    main()
