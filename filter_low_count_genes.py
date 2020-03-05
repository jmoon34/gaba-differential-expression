import numpy as np

def filter_low_count_genes(gene_order,counts_matrix,count_threshold=0):
    """
    Filters out genes with low counts across all libraries

    Parameters
    ----------
    gene_order : List[str]
        A list of genes that correspond to the rows of the counts matrix
    counts_matrix : np.ndarray
        A 2 x 2 array where the rows are gene counts across libraries and the
        columns are the different libraries. 
    count_threshold : int
        Genes with inter-library counts below this threshold will be removed
        from the counts matrix and the gene order list. 

    Returns
    -------
    (List[str],List[str],np.ndarray)
        A tuple containing two lists of the remaing and removed genes and the
        filtered counts matrix.
    """
    gene_count = np.nansum(counts_matrix,axis=1)
    low_count_rows = np.where(gene_count <= count_threshold)[0]
    filtered_matrix = np.delete(counts_matrix,low_count_rows,axis=0)
    filtered_genes = np.delete(gene_order,low_count_rows,axis=0)
    removed_genes = [gene_order[i] for i in low_count_rows]
    return filtered_genes,removed_genes,filtered_matrix
