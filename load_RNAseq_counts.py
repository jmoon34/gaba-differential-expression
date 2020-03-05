import numpy as np

def load_RNAseq_counts(filename,delimiter='\t'):
    """
    Loads a RNA-seq counts file from disk.

    Parameters
    ----------
    filename : str
        The name of the RNA-seq counts file to parse
    delimiter : str
        The character by which the parser will separate the columns. 
        It is recommended to pass '\t' if your file is tab-delimited,
        or ' ' if your file is space-delimited.

    Returns
    -------
    (np.ndarray,List[str],List[str])
        Returns a tuple containing a 2 x 2 np.ndarray of RNA-seq counts, where
        each column represents a library and each row contains gene counts for a
        given gene across all libraries, and two lists containing gene and replicate names. 
    """
    rep_order = []
    with open(filename,'r') as handle:
        columns = zip(*[line.split(delimiter) for line in handle])
        columns = [list(column) for column in columns]
        for column in columns:
            name = column[0]
            column.remove(name)
            for i in range(len(column)):
                try:
                    column[i] = int(column[i].strip('\n'))
                except ValueError:
                    column[i] = column[i].strip('\n')
            if 'gene' in name.strip('#').strip('\n').lower():
                gene_order = column
                continue
            else:
                rep_order.append(name.strip('#').strip('\n').lower())

    return np.array([columns[i] for i in range(1,len(columns))]).T,gene_order,rep_order

