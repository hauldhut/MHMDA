# MHMDA
MHMDA: Integration of multiple microRNA similarity networks and a disease similarity network improves the prediction of disease-associated microRNAs

![Network construction](https://github.com/hauldhut/MHMDA/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: To store performance evaluation results
  - **Prediction**: To store collected evidence
- **Figures**: To store generated figures from the simulation results

## How to run
- Install R packages
  - *RandomWalkRestartMH, multiMiR, easyPubMed, igraph, foreach, doParallel, ROCR, ggplot2, Metrics, hash*
- Download the repo
- Follow instructions in the folder **Code** to run
  
- *Note*: For large and complex networks (i.e., miRNA multiplex networks, and multiplex-hetergeneous networks of diseases and miRNAs), it is recommended to run on a multi-core and at least 16 GB RAM computer
