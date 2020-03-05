def load_gene_info(ucsc_file):
    """
    Loads a file output from UCSC genome browser from disk.

    Parameters
    ----------
    ucsc_file : str
        The name of the file to parse. This file must be tab-delimited and have
        the following columns: 
            name
            chrom
            strand
            txStart
            txEnd
            cdsStart
            cdsEnd
            alignID
            mRNA
            geneSymbol   
    Returns
    -------
    gene_info : Dict[str,Dict[str,any]]
        A superdict containing information about each gene transcript.
    """
    gene_info = {}
    with open(ucsc_file,'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            name = line.split('\t')[0].strip()
            accession = line.split('\t')[8].strip()
            gene_symbol = line.split('\t')[9].strip()
            chr = line.split('\t')[1].strip()
            strand = line.split('\t')[2].strip()
            start = int(line.split('\t')[3].strip())
            end = int(line.split('\t')[4].strip())
            length = end-start

            gene_info[name] = {'name':name,
                               'accession':accession,
                               'gene_symbol':gene_symbol,
                               'chr':chr,
                               'strand':strand,
                               'start':start,
                               'end':end,
                               'length':length}
    handle.close()
    return gene_info
