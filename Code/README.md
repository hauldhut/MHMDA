# This folder contains all source code used in the manuscript

- **miRNA_M_LOOCV_Params.R**: Run Leave-One-Out cross-validation for miRNA Monoplex/Multiplex Networks (for RWRMDA, MHMDA-M)
  - RWRMDA: *M1, M2, M3*: For miRNA Monoplex Networks, i.e., MonoNet_miRWalk, MonoNet_TargetScan and MonoNet_Integrated, respectively
  - MHMDA-M: *M12*: For miRNA Multiplex Network, i.e., MultiNet_miRNA
 
- **miRNA_MH_LOOCV_Params.R**: Run Leave-One-Out cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks of Diseases and miRNAs (for RWRHMDA, MHMDA-MH)
  - RWRHMDA: *H1, H2, H3*: For Heterogeneous Networks which connect a Disease Similarity Network with a miRNA Monoplex Network (i.e., M1, M2, M3, respectively)
  - MHMDA-MH: *MH*: For Multiplex-Heterogeneous Networks which connect a Disease Similarity Network with a miRNA Multiplex Network (i.e., M12)

- **miRNA_Summarize_AUROC_Final.R**: To investigate the prediction performance in terms of AUROC resulted from **miRNA_M_LOOCV_Params.R** and **miRNA_MH_LOOCV_Params.R** by parameters
- **miRNA_Summarize_AUPRC_Final.R**: To investigate the prediction performance in terms of AUPR resulted from **miRNA_M_LOOCV_Params.R** and **miRNA_MH_LOOCV_Params.R** by parameters

- **miRNA_Summarize_AUROC_Params.R**: To summarize prediction performance in terms of AUROC of MHMDA-M and MHMDA-MH according to the change of parameters**

- **miRNA_MH_Predict_Evidence.R**: To predict and select top 20 highly ranked miRNAs for each disease, then find evidence supporting the promissing disease-miRNA associations from existing databases and literature

- MAMFGAT_Comparison: For the comparison with MAMFGAT

