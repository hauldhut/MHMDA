# This folder contains all source code used in the manuscript

- **miRNA_M_LOOCV_Params.R**: Run Leave-One-Out cross-validation for miRNA Monoplex/Multiplex Networks
  - *M1, M2, M3*: For miRNA Monoplex Networks, i.e., Mono_miRWalk, Mono_TargetScan and MonoNet_Integrated, respectively
  - *M12*: For miRNA Multiplex Network, i.e., MultiNet_miRNA
 
- **miRNA_MH_LOOCV_ROC_Params.R**: Run Leave-One-Out cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks of Diseases and miRNAs
  - *H1, H2, H3*: For Heterogeneous Networks which connect a Disease Similarity Network with a miRNA Monoplex Network (i.e., M1, M2, M3, respectively)
  - *MH*: For Multiplex-Heterogeneous Networks which connect a Disease Similarity Network with a miRNA Multiplex Network (i.e., M12)

- **miRNA_Summarize_AUC_Params.R**: To investigate the prediction performance in terms of AUC resulted from **miRNA_M_LOOCV_Params.R** and **miRNA_MH_LOOCV_ROC_Params.R** by parameters

- **miRNA_Summarize_AUC_Params.R**: To summarize prediction performance of the RWR-based methods on a miRNA multiplex network and a multiplex-heterogeneous network according to the change of parameters**
  
- **miRNA_Summarize_AUC_ROC.R**: To summarize and draw ROC curves resulted from **miRNA_M_LOOCV_Params.R** and **miRNA_MH_LOOCV_ROC_Params.R**

- **miRNA_MH_Predict_Evidence.R**: To predict and select top 20 highly ranked miRNAs for each disease, then find evidence supporting the promissing disease-miRNA associations from existing databases and literature

