Method = "MH"#MH/H1/H2/H3

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)

setwd("~/Manuscripts/100MHMDA/Code")

#Mono_miRWalk.txt
#Mono_TargetScanNet.txt

miRNA1 <- read.delim("../Data/Mono_miRWalk.txt",header = FALSE)
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

if(Method == "MH"){
  miRNA_MultiplexObject <- create.multiplex(list(miRNA1.g,miRNA2.g),Layers_Name = c("miRNA1","miRNA2"))  
  tau1 = 1
  tau2 = 1
  tau <- c(tau1, tau2)
}else{
  tau <-c(1)
  if(Method == "H1"){
    miRNA_MultiplexObject <- create.multiplex(list(miRNA1.g),Layers_Name = c("miRNA"))
  }else if(Method == "H2"){
    miRNA_MultiplexObject <- create.multiplex(list(miRNA2.g),Layers_Name = c("miRNA"))
  }else{
    miRNA_MultiplexObject <- create.multiplex(list(miRNA3.g),Layers_Name = c("miRNA"))
  }
}


#Add disease nw
#DiseaseSimNet_OMIM.txt
#DiseaseSimNet_HPO.sif
#DiseaseSimNet_GeneNet.txt

disease <- read.delim("../Data/DiseaseSimNet_OMIM.txt",header = FALSE)

disease.frame <- data.frame(disease[[1]], disease[[3]])
disease.weight = disease[[2]]

disease.g <- graph.data.frame(d = disease.frame, directed = FALSE)
E(disease.g)$weight <- disease.weight

disease_MultiplexObject <- create.multiplex(list(disease.g),
                                            Layers_Name = c("disease"))

#Add MDRelation
BipartiteNet="miR2Disease"#miR2Disease/HMDD
if(BipartiteNet=="miR2Disease"){
  MD.frame <- read.csv("../Data/Phenotype2miRNAs_miR2Disease.csv", header = TRUE)  
}else{
  MD.frame <- read.csv("../Data/Phenotype2miRNAs_HMDD.csv", header = TRUE)  
}

MD.frame <- MD.frame[which(MD.frame$miRNA %in% miRNA_MultiplexObject$Pool_of_Nodes),]
MD.frame <- MD.frame[which(MD.frame$disease %in% disease_MultiplexObject$Pool_of_Nodes),]

delta = 0.5
gamma = 0.5

for (delta in c(0.1, 0.3, 0.5, 0.7, 0.9)){
  if (delta == 0.5) next
  
  cat("delta = ", delta,"gamma = ", gamma,"\n")

  #func
  do_something <- function(miRNA_MultiplexObject,disease_MultiplexObject,
                           EDRelation,SeedEnhancer, seeddisease, prd_miRNAs) {
    
    #Create multiplex-heterosgenous nw
    
    EDRelation_miRNA <- EDRelation[which(EDRelation$miRNA %in% miRNA_MultiplexObject$Pool_of_Nodes),]
    
    #Create multiplex-heterosgenous nw
    miRNA_Disease_Net <- create.multiplexHet(miRNA_MultiplexObject, disease_MultiplexObject, 
                                                EDRelation_miRNA)
    
    miRNAHetTranMatrix <- compute.transition.matrix(miRNA_Disease_Net,lambda = 0.5, delta1=delta, delta2=0.5)
    
    #compute 
    cat(seeddisease, SeedEnhancer, "\n")
    RWRH_miRNA_Disease_Results <- Random.Walk.Restart.MultiplexHet(miRNAHetTranMatrix,
                                                                      miRNA_Disease_Net,SeedEnhancer,
                                                                      seeddisease, r = gamma)
    
    
    #create labels for miRNA results
    tf = RWRH_miRNA_Disease_Results$RWRMH_Multiplex1
    tf$labels <- ifelse(tf$NodeNames==prd_miRNA, 1, 0)
    
    # calculating AUC for each miRNA of a disease
    resultspred = prediction(tf$Score, tf$labels)
    
    pauc.perf = performance(resultspred, measure = "auc")
    return(list(pauc.perf@y.values[[1]],data.frame(Scores=tf$Score, Labels=tf$labels)))
    # return(pauc.perf@y.values[[1]])
  }
  
  
  no_cores <- 4
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  #loop through
  res <- foreach(i=1:length(MD.frame$miRNA), .combine = rbind) %dopar% {#length(MD.frame$miRNA)730/1451/2172/2893
  
  # res = NULL
  # for(i in 1:5){#length(MD.frame$miRNA)730/1451/2172/2893
    cat(i,"\n")
    
    library(RandomWalkRestartMH)
    library(igraph)
    library(ROCR)
    
    prd_miRNA = MD.frame$miRNA[[i]]
    seeddisease = MD.frame$disease[[i]]
    
    disease_relation = MD.frame[which(MD.frame$disease==seeddisease),]
    
    SeedEnhancer = disease_relation$miRNA[-c(which(disease_relation$miRNA==prd_miRNA))]
    
    EDRelation <- EDRelation <- MD.frame[-with(MD.frame, which(miRNA %in% prd_miRNA & disease %in% seeddisease)),]
    
    res <- do_something(miRNA_MultiplexObject,disease_MultiplexObject,
                      EDRelation,SeedEnhancer, seeddisease, prd_miRNAs)
    
  }
  
  dim(res)
  
  stopCluster(cl)
  
  #Store auc by Trial --> From this, auc by Disease can be obtain by aggregate function when summarizing AUCs
  MD.frame$auc <- unlist(res[,1])
  
  
  Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_ROC.csv")
  cat(Result_byTrialFile,"\n")
  write.csv(MD.frame,Result_byTrialFile, row.names = FALSE, quote = FALSE)
  
  aucavgbyTrial = round(mean(MD.frame$auc),3)
  aucavgbyTrial.sd = round(sd(MD.frame$auc),3)
  
  
  res.final = NULL
  for(i in 1:nrow(res)){
    res.final = rbind(res.final, res[i,2][[1]])
  }
  
  
  Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
  Result_by_Score_n_Label_SummaryFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label_Summary.csv")
  
  cat(Result_by_Score_n_LabelFile,"\n")
  write.csv(res.final,Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  # write.csv(data.frame(Scores = Scores, Labels = Labels),Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  
  resultspred = prediction(res.final$Scores, res.final$Labels)
  auc.perf = performance(resultspred, measure = "auc")
  aucavgbyAll = round(auc.perf@y.values[[1]],3)
  
  cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")
  write.csv(data.frame(Method=Method,gamma=gamma,delta=delta,aucavgbyAll=aucavgbyAll,aucavgbyTrial=aucavgbyTrial,aucavgbyTrial.sd=aucavgbyTrial.sd),Result_by_Score_n_Label_SummaryFile, row.names = FALSE, quote = FALSE)
  
  end_time <- Sys.time()
  end_time - start_time
}
