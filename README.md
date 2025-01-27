BC-SELECT is a computational tool adapted from SELECT (Cell, 2021) to predict drug responses in early-stage breast cancer patients using transcriptome data. 

It involves two main steps: (1) a 'gene pair identification' step and (2) 'fine-tuning hyperparameters and validation performance'. 
The first part leverages large-scale cell-line datasets and breast cancer patient cohorts, including Breast Cancer TCGA, METABRIC, and SCAN-B, to identify 
clinically relevant genetic interactions between gene pairs A and B. Here, gene A represents the target of a given therapy (e.g., targeted therapy or immunotherapy), 
while multiple partner genes (Bs) are identified based on their interactions with A. Specifically, drug efficacy is explained through mechanisms of genetic interaction. 
To explore this, we constructed a comprehensive library of clinically relevant candidate partner genes (Bs) involved in synthetic lethality (SL), synthetic dosage lethality (SDL), 
or synthetic rescue (SR) interactions with a given drug target (A). 

The second part focuses on fine-tuning hyperparameters using the results from the identified gene pairs in step (1) 
and three distinct breast cancer tuning cohorts. Once the optimal settings are established, we proceed to the model performance validation phase with unseen evaluation breast cancer cohorts. 
In this phase, the expression levels of partner B genes are used to generate a patient's predicted response score (likelihood) to treatments targeting gene A(s).

To assess drug efficacy for each patient, we calculate the predicted response score based on the expression levels of the identified gene pairs (A(s) and Bs) and 
the curated genes' expression obtained through NanoString.

We use separate modules for targeted therapies and immunotherapies, so the "gene pair identification" and "performance validation" folders are organized accordingly, 
with one set for targeted therapy and the other for immunotherapy.


