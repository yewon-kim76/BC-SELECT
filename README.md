BC-SELECT is a computational tool adapted from SELECT (Cell, 2021) to predict drug responses in early-stage breast cancer patients. 
It involves two main steps: (1) a 'training' step that constructs a library of clinically relevant candidate SL/SDL/SR partner genes for a given drug target, 
and (2) a 'testing/validation' step that calculates predicted response scores for patients in unseen studies based on the relative expression of these candidate partner genes.

For targeted therapy, we used three filters to identify SL/SDL candidate partner genes. 
A gene must pass the following sequentially: 1) cell line test, 2) survival test, and 3) evolutionary analysis. 
For the immunotherapy module, we replaced the cell line test with a hypergeometric screen to obtain 
a sufficient number of immune-related SR candidate partner genes.

Gene expression is assessed only within datasets (not across datasets, to minimize batch effect confounding) and divided into tertiles, 
with upregulated genes in the top tertile and downregulated genes in the bottom tertile, based on relative expression levels across patients. 
Using the bin_matrix function, we categorize these into tertiles (q2 = 0, 1, and 2) and measure relative gene expression levels.
