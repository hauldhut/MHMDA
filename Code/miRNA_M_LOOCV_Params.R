Method = "M12"#M1/M2/M3/M12

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(ROCR)

setwd("~/Manuscripts/100MHMDA/Code")

#HomomiRWalkNet.txt
#HomoTargetScanNet.txt

miRNA1 <- read.delim("../Data/Mono_miRWalkNet.txt",header = FALSE)
miRNA1.frame <- data.frame(miRNA1[[1]], miRNA1[[3]])
miRNA1.g <- graph.data.frame(d = miRNA1.frame, directed = FALSE)
miRNA1.weight = miRNA1[[2]]
E(miRNA1.g)$weight <- miRNA1.weight

miRNA2 <- read.delim("../Data/Mono_TargetScan.txt",header = FALSE)
miRNA2.frame <- data.frame(miRNA2[[1]], miRNA2[[3]])
miRNA2.g <- graph.data.frame(d = miRNA2.frame, directed = FALSE)
miRNA2.weight = miRNA2[[2]]
E(miRNA2.g)$weight <- miRNA2.weight

#"per-edge average" integration
#Normalize
# hist(miRNA1.weight)#Already normalized
# hist(miRNA2.weight)#Already normalized
#Integrate
miRNA3 = merge(miRNA1, miRNA2, by=c("V1","V3"), all = TRUE)
miRNA3[is.na(miRNA3$V2.x),3] = 0
miRNA3[is.na(miRNA3$V2.y),4] = 0
miRNA3$V2 = (miRNA3$V2.x+miRNA3$V2.y)/2

miRNA3.frame <- data.frame(miRNA3[[1]], miRNA3[[2]])
miRNA3.g <- graph.data.frame(d = miRNA3.frame, directed = FALSE)
miRNA3.weight = miRNA3[[5]]
E(miRNA3.g)$weight <- miRNA3.weight


if(Method == "M12"){
  miRNA_MultiplexObject <- create.multiplex(list(miRNA1.g,miRNA2.g),Layers_Name = c("miRNA1","miRNA2"))  
  tau1 = 1
  tau2 = 1
  tau <- c(tau1, tau2)
}else if(Method == "M1"){
  miRNA_MultiplexObject <- create.multiplex(list(miRNA1.g),Layers_Name = c("miRNA"))
  tau = c(1)
}else if(Method == "M2"){
  miRNA_MultiplexObject <- create.multiplex(list(miRNA2.g),Layers_Name = c("miRNA"))
  tau = c(1)
}else{
  miRNA_MultiplexObject <- create.multiplex(list(miRNA3.g),Layers_Name = c("miRNA"))
  tau = c(1)
}


delta = 0.5
gamma = 0.5

for (gamma in c(0.1, 0.3, 0.5, 0.7, 0.9)){
  # if (delta == 0.5) next
  cat("delta = ", delta,"gamma = ", gamma,"\n")
  
  
  AdjMatrix_miRNA <- compute.adjacency.matrix(miRNA_MultiplexObject, delta = delta)
  AdjMatrixNorm_miRNA <- normalize.multiplex.adjacency(AdjMatrix_miRNA)
  
  #Add MDRelation
  BipartiteNet="HMDD"#miR2Disease/HMDD
  if(BipartiteNet=="miR2Disease"){
    MD.frame <- read.csv("../Data/Phenotype2miRNAs_miR2Disease.csv", header = TRUE)  
  }else{
    MD.frame <- read.csv("../Data/Phenotype2miRNAs_HMDD.csv", header = TRUE)  
  }
  
  MD.frame <- MD.frame[which(MD.frame$miRNA %in% miRNA_MultiplexObject$Pool_of_Nodes),]
  
  cat("PROCESSING Method=",Method,'BipartiteNet=',BipartiteNet,"\n")
  
  res = NULL
  for (i in 1:length(MD.frame$miRNA)) {
    
    prd_miRNA = MD.frame$miRNA[i]
    seeddisease = MD.frame$disease[i]
  
    cat(i,"/",length(MD.frame$miRNA),":",seeddisease,"\n")
    
    disease_relation = MD.frame[which(MD.frame$disease==seeddisease),]
    SeedEnhancer = disease_relation$miRNA[-c(which(disease_relation$miRNA==prd_miRNA))]
    
  
    if (length(disease_relation$miRNA)>=2) {
  
      RWR_Enhancer_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_miRNA,
                                                            miRNA_MultiplexObject,
                                                            SeedEnhancer,tau = tau, r = gamma)
      tf = RWR_Enhancer_Results$RWRM_Results
      #create labels for miRNA results
      tf$labels <- ifelse(tf$NodeNames==prd_miRNA, 1, 0)
      
      # calculating AUC for each miRNA of a disease
      resultspred = prediction(tf$Score, tf$labels)
  
      pauc.perf = performance(resultspred, measure = "auc")
  
      MD.frame$auc[i] <- pauc.perf@y.values[[1]]
      
      res = rbind(res, data.frame(Scores=tf$Score, Labels=tf$labels))
      
    } else {
      MD.frame$auc[i] <- NA #Diseases have only one known associated miRNA
    }
  }
  
  dim(res)
  
  #Store auc by Trial --> From this, auc by Disease can be obtain by aggregate function when summarizing AUCs
  MD.frame = MD.frame[!is.na(MD.frame$auc),]
  
  Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_ROC.csv")
  cat(Result_byTrialFile,"\n")
  write.csv(MD.frame,Result_byTrialFile, row.names = FALSE, quote = FALSE)
  
  aucavgbyTrial = mean(MD.frame$auc)
  aucavgbyTrial.sd = sd(MD.frame$auc)
  
  #Store Scores and Labels for calculating and drawing AUC and ROC
  
  Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
  Result_by_Score_n_Label_SummaryFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label_Summary.csv")
  cat(Result_by_Score_n_LabelFile,"\n")
  write.csv(res,Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  
  resultspred = prediction(res$Scores, res$Labels)
  auc.perf = performance(resultspred, measure = "auc")
  aucavgbyAll = auc.perf@y.values[[1]]
  
  cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")
  write.csv(data.frame(Method=Method,gamma=gamma,delta=delta,aucavgbyAll=aucavgbyAll,aucavgbyTrial=aucavgbyTrial,aucavgbyTrial.sd=aucavgbyTrial.sd),Result_by_Score_n_Label_SummaryFile, row.names = FALSE, quote = FALSE)
  
  end_time <- Sys.time()
  end_time - start_time
}