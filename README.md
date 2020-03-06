*Refer to /Reports/Moon_Final_Project.docx for more information*

# Differential Expression of the SLC12A Family in ES and NPC Cells

## Abstract
GABA is a neurotransmitter that changes role between embryonic development and maturity in the brain of humans.  Using RNA-sequencing data of ES cells and NPCs, the differential expression of a portion of the genes that are responsible for this change was observed.  After adjusting the data and performing various testing methods as well as multiple testing correction, the study has found a total of 11 genes that had significant changes in means.  Although the results agree with previous literature on specific genes that are known to be differentially expressed, the limitations of the testing methods, as well as the data itself, suggests more work is needed to fully understand the changes observed in all the other genes.    

## Hypotheses
Since in an immature brain, there is a down regulation of KCC2 and an upregulation of NKCC1,  the null hypothesis is that there is no difference in mean in counts for the genes that are responsible for KCC2 and NKCC1.  The alternative hypothesis is that for ES cells, the mean for KCC2 will be lower and the mean for NKCC1 will be higher than NPC cells. 

## Results
At the alpha level of 0.05, Fisher’s exact test and Poisson LRT yielded 5 and 7 significant p-values respectively after multiple testing correction.  Five genes that code for SLC12A6, SLC12A7 and SLC12A8 were upregulated in ES cells according to Fisher’s exact test.  Likewise, seven genes that code for SLC12A2 (NKCC1), SLC12A5 (KCC2), SLC12A6, SLC12A7, SLC12A8 and SLC12A9 were significantly changed in ES cells according to the LRT.  Three of them were most likely downregulated with the treatment considering they had 0 counts for both NPC replicates, but four had relatively similar counts in all four samples, making it difficult to discern which direction gene expression took during the treatment.  

## Method
After transformation of data, p-values were obtained using Welch’s t-test, Fisher’s exact test, Kolmogorov-Smirnov and Poisson log-likelihood ratio test.  Finally, multiple testing correction using the Benjamini-Hochberg method was performed to adjust the p-values for more accurate analysis.

![Volcano plot (Welch's t-test)](/Diagrams/volcano_ttest.png?raw=true "Volcano plot (Welch's t-test)")

![Volcano plot (K-S test)](/Diagrams/volcano_KS.png?raw=true "Volcano plot (K-S test)")

![Volcano plot (Fisher's Exact Test)](/Diagrams/volcano_Fisher.png?raw=true "Volcano plot (Fisher's Exact Test)")

![Volcano plot (LRT)](/Diagrams/volcano_LRT.png?raw=true "Volcano plot (Poisson Log-Likelihood Ratio Test)")

## Conclusion
Regardless of the specific direction of regulation from ES cells to NPCs, both Fisher’s exact test and LRT show statistically significant changes in counts in a given portion of the genes that code for the SLC12A family.  Specifically, genes that code for NKCC1 and KCC2 both show significant changes, which supports the claim that there is a difference in intracellular chloride concentration, which in turn illustrates the different role that GABA plays as a neurotransmitter in ES cells and in adult brain cells.   


## Reference
Y. Ben-Ari. The GABA Excitatory/Inhibitory Developmental Sequence: A Personal Journey.   Neuroscience 279 (2014) 187-219
