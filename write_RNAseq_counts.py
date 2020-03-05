def write_RNAseq_counts(gene_order,lib_order,counts_matrix,filename,delimiter='\t'):
    """
    Writes an RNA-seq counts dict back to disk as a counts file.
    
    Parameters
    ----------
    gene_order : List[str]
        The list of gene names that correspond to each row
    lib_order : List[str]
        The list of library names that correspond to each column
    counts_matrix : Dict[str,List[int]]
        The RNA-seq counts matrix where each row represents gene counts across
        libraries and the columns are the different libraries. 
    filename : str
        The directory and name of the file to write to. (i.e.
        'my_folder/RNA_seq.txt')
    delimiter : str
        Character by which the writer will separate the columns.
        It is recommended to pass '\t' if the desired file output is tab-delimited,
        or ' ' if your desired file output is space-delimited.
    """
    with open(filename,'wb') as handle:
        #write header
        handle.write('#gene_name')
        for lib_name in lib_order:
            if lib_name == lib_order[-1]:
                handle.write(delimiter+'%s\n'%(lib_name))
            else:
                handle.write(delimiter+'%s'%(lib_name))
        #write counts data
        for i in range(len(counts_matrix[:,0])):
            handle.write('%s'%(gene_order[i]))
            for j in range(len(counts_matrix[0,:])):
                if lib_order[j] == lib_order[-1]:
                    handle.write(delimiter+'%s\n'%(counts_matrix[i,j]))
                else:
                    handle.write(delimiter+'%s'%(counts_matrix[i,j]))
