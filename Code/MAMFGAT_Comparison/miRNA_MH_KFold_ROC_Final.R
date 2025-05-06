Method = "MH"#MH/H1/H2/H3

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)

setwd("~/Manuscripts/100MHMDA/Code")

#HomomiRWalkNet.txt
#HomoTargetScanNet.txt

miRNA1 <- read.delim("../Data/MonoNet_miRWalk.txt",header = FALSE)
miRNA1.frame <- data.frame(miRNA1[[1]], miRNA1[[3]])
miRNA1.g <- graph.data.frame(d = miRNA1.frame, directed = FALSE)
miRNA1.weight = miRNA1[[2]]
E(miRNA1.g)$weight <- miRNA1.weight

miRNA2 <- read.delim("../Data/MonoNet_TargetScan.txt",header = FALSE)
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


#add func for RWR on multiplex-heter nw
do_something <- function(miRNA_MultiplexObject,disease_MultiplexObject,
                         MDRelation,SeedmiRNA, seeddisease, prd_miRNAs) {
  
  #Create multiplex-heterosgenous nw
  
  MDRelation_disease <- MDRelation[which(MDRelation$miRNA %in% miRNA_MultiplexObject$Pool_of_Nodes),]
  
  #Create multiplex-heterosgenous nw
  miRNA_Disease_Net <- create.multiplexHet(miRNA_MultiplexObject, disease_MultiplexObject, 
                                              MDRelation_disease)
  
  miRNA_Disease_Net_TranMatrix <- compute.transition.matrix(miRNA_Disease_Net)
  
  #compute 
  # tau <- c(1,1)
  Ranking_Results <- Random.Walk.Restart.MultiplexHet(miRNA_Disease_Net_TranMatrix,
                                                                    miRNA_Disease_Net,SeedmiRNA,
                                                                    seeddisease, r = gamma)
  
  #create labels for ranking results
  tf = Ranking_Results$RWRMH_Multiplex1
  
  tf$labels <- ifelse(tf$NodeNames %in% prd_miRNAs, 1, 0)
  
  # calculating AUC
  resultspred = prediction(tf$Score, tf$labels)
  
  pauc.perf = performance(resultspred, measure = "auc")
  return(list(pauc.perf@y.values[[1]],data.frame(Scores=tf$Score, Labels=tf$labels)))
}

#count disease for each disease
sub_sum <- aggregate(miRNA~disease, data=MD.frame, FUN=function(x) c(count=length(x)))

#extract disease with only k or more diseases
k=5
sub_sum <- sub_sum[which(sub_sum$miRNA>=k),]
sub_sum$disease_no <- c(1:length(sub_sum$disease))

#extract MD.frame with only disease from sub_sum
MD.frame1 <- MD.frame[which(MD.frame$disease %in% sub_sum$disease),]
rownames(MD.frame1) <- NULL #reset frame index

#func to assign k groups for each set of disease-miRNA, 
#as well as increment group no. for each group (k=3)
assign_group_no <- function(sub_sum,MD.frame1,k) {
  
  #set an empty data frame for a new MD.frame
  mylist.names <- c("disease","miRNA", "disease_no","group_no")
  MD.frame2 <- sapply(mylist.names,function(x) NULL)
  
  
  for (j in 1:length(sub_sum$disease)) {
    
    count = sub_sum$miRNA[[j]]
    
    if(count<k) next
    
    set_no = floor(count/k)
    
    print(paste(k,count, set_no))
    
    group_vec = vector()
    for(gi  in 1:(k-1)){
      group_vec = c(group_vec,rep(gi,set_no))
    }
    group_vec = c(group_vec,rep(k,count-set_no*(k-1)))

    # group_vec <- c(rep(1,set_no), rep(2,set_no), rep(3,set_no), rep(4,set_no),rep(5,set_no), rep(6,set_no),rep(7,set_no), rep(8,set_no), rep(9,set_no),rep(10,(count-set_no*9)))
    
    subset <- MD.frame[which(MD.frame$disease==sub_sum$disease[[j]]),]
    subset$disease_no <- rep(j,count)
    subset$group <- group_vec
    
    MD.frame2 <- rbind(MD.frame2,subset)
  }
  return(MD.frame2)
}

#assign each group of n/k miRNA-disease with an group id
out <- assign_group_no(sub_sum,MD.frame1,k)
out$group <- (out$disease_no-1)*k+out$group
MD.frame2 <- out

#set an empty data frame for a new MD.frame
auc_results <- sapply(c("disease","group","auc"),function(x) NULL)

for (i in 1:max(MD.frame2$group)) {
  seeddisease = unique(MD.frame2$disease[which(MD.frame2$group==i)])
  auc_results$disease[i] <- seeddisease
  auc_results$group[i] <- i
}

#set up paralell processing (adjust the no_cores as per running system)
no_cores <- 6
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#loop through
res <- foreach(i = 1:max(MD.frame2$group), .combine = rbind) %dopar% {
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  prd_miRNAs = MD.frame2$miRNA[which(MD.frame2$group==i)]
  seeddisease = unique(MD.frame2$disease[which(MD.frame2$group==i)])
  
  disease_relation = MD.frame2[which(MD.frame2$disease==seeddisease),]
  SeedmiRNA = disease_relation$miRNA[-c(which(disease_relation$miRNA %in% prd_miRNAs))]
  
  # get bipartite graph without prd_miRNAs - disease linkages
  MDRelation <- MD.frame2[-with(MD.frame2, which(miRNA %in% prd_miRNAs & disease %in% seeddisease)),][1:2]
  
  res <- do_something(miRNA_MultiplexObject,disease_MultiplexObject,
                    MDRelation,SeedmiRNA, seeddisease, prd_miRNAs)
  res <- append(res, seeddisease)
}

dim(res)
stopCluster(cl)


df.res = data.frame(trial = c(1:nrow(res)),auc = unlist(res[,1]))
df.res.final = merge(df.res, MD.frame2[c("disease","group")], by.x="trial", by.y="group", all.x = TRUE)
df.res.final = unique(df.res.final)

Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_ROC_KFold5.csv")
cat(Result_byTrialFile,"\n")
write.csv(df.res.final,Result_byTrialFile, row.names = FALSE, quote = FALSE)

aucavgbyTrial = round(mean(df.res.final$auc),3)
aucavgbyTrial.sd = round(sd(df.res.final$auc),3)


res.final = NULL
for(i in 1:nrow(res)){
  res.final = rbind(res.final, res[i,2][[1]])
}

Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_Score_n_Label_KFold5.csv")
cat(Result_by_Score_n_LabelFile,"\n")
write.csv(res.final,Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
# write.csv(data.frame(Scores = Scores, Labels = Labels),Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)

library(ROCR)
resultspred = prediction(res.final$Scores, res.final$Labels)
auc.perf = performance(resultspred, measure = "auc")
aucavgbyAll = round(auc.perf@y.values[[1]],3)

cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")

# Calculate performance for Precision-Recall curve
pr.perf <- performance(resultspred, "prec", "rec")

# Extract Recall (x-axis) and Precision (y-axis) values for plotting
Recall = pr.perf@x.values[[1]]    # Recall values
Precision = pr.perf@y.values[[1]] # Precision values

# Calculate AUPRC (Area Under Precision-Recall Curve)
aupr.perf = performance(resultspred, measure = "aucpr") # 'aucpr' is the measure for AUPR
aupravgbyAll = round(aupr.perf@y.values[[1]], 3)       # Round to 3 decimal places

# Optional: Print the AUPR value
print(paste("AUPRC:", aupravgbyAll))

end_time <- Sys.time()
timediff = end_time - start_time
print(timediff)
