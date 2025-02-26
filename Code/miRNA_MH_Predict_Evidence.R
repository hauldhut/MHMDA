Method = "MH"#MH/H1/H2/H3

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)

setwd("~/Manuscripts/100MHMDA/Code")

#DiseaseSimNet_OMIM.txt
#DiseaseSimNet_HPO.sif
#DiseaseSimNet_GeneNet.txt

miRNA1 <- read.delim("../Data/MonoNet_miRWalkNet.txt",header = FALSE)
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


#Create multiplex-heterosgenous nw
miRNA_Disease_Net <- create.multiplexHet(miRNA_MultiplexObject, disease_MultiplexObject, MD.frame)

miRNAHetTranMatrix <- compute.transition.matrix(miRNA_Disease_Net, lambda = 0.5)

#pick unique disease
unique_disease <- unique(MD.frame$disease)

h.d2ranking = hash()

for(di in 1:length(unique_disease)){ 
  # di = 1
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  seeddisease = unique_disease[[di]]
  
  cat("==> Disease ", di, "/", length(unique_disease), ": ", seeddisease,"\n")
  
  disease_relation = MD.frame[which(MD.frame$disease==seeddisease),]
  
  SeedmiRNA = disease_relation$miRNA[c(which(disease_relation$disease==seeddisease))]
  
  #compute 
  RWRH_miRNA_Disease_Results <- Random.Walk.Restart.MultiplexHet(miRNAHetTranMatrix,
                                                                    miRNA_Disease_Net,SeedmiRNA,
                                                                    seeddisease, r = 0.5)
  
  tf = RWRH_miRNA_Disease_Results$RWRMH_Multiplex1
  h.d2ranking[[seeddisease]] = tf  
}

#####################################################
Disease_Info <- read.delim("../Data/Phenotype2Genes_Full.txt",header = TRUE)

#Summarize for each k
maxk=20
h.K2TotalEvidencedNewAssoc = hash()
h.K2topKEnhEvidence = hash()

#Store top k predictions
for(k in seq(10, maxk, by = 10)){
  
  cat("k =",k,"\n")
  df.topKEnh <- data.frame(matrix(ncol = 2+k, nrow = 0))
  
  for(mimid in unique_disease){ 
    tf = h.d2ranking[[mimid]]
    df.topKEnh[nrow(df.topKEnh)+1,] = c(as.character(mimid),as.character(Disease_Info[Disease_Info$MIMID==mimid,]$Name),as.vector(tf$NodeNames[1:k]))
  }
  topKfile = paste0("../Results/Prediction/",Method,"_",BipartiteNet,"_predict_top",k,".txt")
  write.table(df.topKEnh, topKfile, na ="", row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
}

##################
library('multiMiR')

# df_combined <- get_multimir(disease.drug="Test", table='disease.drug')@data

# df_combined <- data.frame()


# Initialize an empty data frame with specified column names
col_names <- c("mature_mirna_id", "disease_drug", "database", "paper_pubmedID")
df_combined <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(df_combined) <- col_names

k=20
for(mimid in unique_disease){ #c("MIM211980","MIM603956")/unique_disease
  dname_full1 = Disease_Info[Disease_Info$MIMID==mimid,]$Name #MELANOMA, CUTANEOUS MALIGNANT, SUSCEPTIBILITY TO, 2; CMM2
  dname_full2 = strsplit(dname_full1,";")[[1]][1] #MELANOMA, CUTANEOUS MALIGNANT, SUSCEPTIBILITY TO, 2
  dname = strsplit(dname_full2,",")[[1]][1] #MELANOMA
  cat("==> ", mimid, dname,"\n")
  
  tf = h.d2ranking[[mimid]]
  
  r <- get_multimir(disease.drug=dname, table='disease.drug')
  df_r = r@data
  df_filtered = df_r[df_r$mature_mirna_id %in% tf$NodeNames[1:k],]
  
  if (nrow(df_filtered)>0){
    df_combined = rbind(df_combined, df_filtered[, col_names])
  } 
    
}
df_combined

topK_Evidencefile = paste0("../Results/Prediction/",Method,"_",BipartiteNet,"_predict_top",k,"_Evidence.txt")
write.table(df_combined, topK_Evidencefile, na ="", row.names=FALSE, col.names = TRUE, sep='\t', quote=FALSE)
topK_Evidencefile = paste0("../Results/Prediction/",Method,"_",BipartiteNet,"_predict_top",k,"_Evidence.csv")
write.csv(df_combined, topK_Evidencefile, na ="", row.names=FALSE, quote=FALSE)

##########
library(easyPubMed)
library(hash)
# my_query <- 'Duc-Hau Le[AU] AND "2022"[PDAT]'

h.mim2name = hash()
h.mim2name[["MIM211980"]] = "lung cancer"
h.mim2name[["MIM176807"]] = "prostate cancer"
h.mim2name[["MIM181500"]] = "schizophrenia"
h.mim2name[["MIM600274"]] = "frontotemporal dementia"
h.mim2name[["MIM155255"]] = "medulloblastoma"

mimid_list = hash::keys(h.mim2name)

for (mimid in mimid_list){
  tf = h.d2ranking[[mimid]]
  
  disease = h.mim2name[[mimid]]
  cat("==>", mimid, disease,"\n")
  
  col_names <- c("rank","miRNA", "NuStudy", "PubMedID")
  df_miRNA_AbstractCount <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(df_miRNA_AbstractCount) <- col_names
  
  rank = 0
  for (miRNA in tf$NodeNames[1:k]){
    my_query = paste0("(",disease,"[Title/Abstract])"," AND ","(",miRNA,"[Title/Abstract])")
    my_res <- get_pubmed_ids(my_query)
    
    rank <- rank + 1
    PMIDList <- ""
    abCount <- 0
    if (my_res$Count>0){
      PMIDList = paste0(my_res$IdList, collapse = ", ")
      abCount <- my_res$Count
      cat(miRNA, abCount,PMIDList, "\n")
      # my_abstracts_txt <- fetch_pubmed_data(my_res, format = "abstract")
    }
    df_miRNA_AbstractCount[nrow(df_miRNA_AbstractCount)+1,] = c(rank, miRNA,abCount,PMIDList)
  }
  
  df_miRNA_AbstractCount
  
  topK_AbstractCountfile = paste0("../Results/Prediction/",Method,"_",BipartiteNet,"_predict_top",k,"_AbstractCount_",disease,".txt")
  write.table(df_miRNA_AbstractCount, topK_AbstractCountfile, na ="", row.names=FALSE, col.names = TRUE, sep='\t', quote=FALSE)
  
}



# #Test multiMiR
# mimid = "MIM211980" #LUNG CANCER
# dname = Disease_Info[Disease_Info$MIMID==mimid,]$Name
# cat(mimid, dname)
# tf = h.d2ranking[[mimid]]
# tf$NodeNames[1:10] #test
# 
# 
# r <- get_multimir(disease.drug=dname, table='disease.drug')
# nrow(r@data)
# head(r@data)
# 
# r_df = r@data
# 
# k=10
# r_filtered = r_df[r_df$mature_mirna_id %in% tf$NodeNames[1:k],]
# r_filtered

# #Find evidence for top k
# library(easyPubMed)
# library(RISmed)
# 
# 
# 
# for(k in seq(10, maxk, by = 10)){
#   #Sol 1
#   library('phenoscanner')
#   
#   h.mimid2evid = hash()
#   di=0
#   for(mimid in unique_disease){ 
#     
#     di = di+1
#     # if(di<94) next
#     cat(mimid,"\n")
#     
#     meshterm = MeSHID2TermMap[[mimid]]
#     meshtraitvec = MeSHID2TraitMap[[mimid]]
#     
#     tf = h.d2ranking[[mimid]]
#     
#     snpquery = as.vector(tf$NodeNames[1:k])
#     res <- phenoscanner(snpquery=snpquery)
#     
#     EnhInfo = res$results
#     EnhInfo = EnhInfo[EnhInfo$rsid %in% snpquery,]
#     
#     if(nrow(EnhInfo)==0) next
#     
#     df.topKEnhEvidence = NULL
#     for(i in 1:nrow(EnhInfo)){
#       trait = EnhInfo$trait[i]
#       # cat(trait,"\n")
#       efoid = EnhInfo$efo[i]
#       efoid = gsub("_", ":",efoid)
#       
#       GWASMeSHIDSet = h.efoid2meshid[[efoid]]
#       if(meshid %in% GWASMeSHIDSet){
#         # cat(meshid,"\t",meshterm,"\t",efoid,"\t",trait,"\t",EnhInfo$rsid[i],"(",EnhInfo$p[i],")","\t",EnhInfo$pmid[i],"\n")
#         df.topKEnhEvidence = rbind(df.topKEnhEvidence, data.frame(meshid = meshid,meshterm=meshterm, trait=trait, efoid=efoid,rsid = EnhInfo$rsid[i],p=EnhInfo$p[i], pmid=EnhInfo$pmid[i]))
#       }
#     }
#     
#     nevid = 0
#     if(!is.null(df.topKEnhEvidence)){
#       h.meshid2evid[[meshid]] = df.topKEnhEvidence  
#       nevid = nrow(df.topKEnhEvidence)
#     }
#     cat(di,"/", length(unique_disease),":",meshid,nrow(EnhInfo),nevid,"\n")
#   }
#   length(h.meshid2evid)
#   
#   df.topKEnhEvidence = NULL
#   for(meshid in keys(h.meshid2evid)){
#     # print(h.meshid2evid[[meshid]])
#     df.topKEnhEvidence = rbind(df.topKEnhEvidence, h.meshid2evid[[meshid]])
#   }
#   dim(df.topKEnhEvidence)
#   topKEvidfile = paste0("../Results/Prediction/",Method,"_MeSHID_Net_Phase",phase,"_Chr",chr,"_ld","_predict_top",k,"_evid.txt")
#   write.table(df.topKEnhEvidence, topKEvidfile, na ="", row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
#   
# }
# 

# 
# ## Test
# library(easyPubMed)
# # my_query <- 'Duc-Hau Le[AU] AND "2022"[PDAT]'
# disease = "LUNG CANCER"
# miRNA = "hsa-miR-15a"
# my_query = paste0("(",disease,"[Title/Abstract])"," AND ","(",miRNA,"[Title/Abstract])")
# 
# my_res <- get_pubmed_ids(my_query)
# my_res$IdList
# my_abstracts_txt <- fetch_pubmed_data(my_res, format = "abstract")
# 
# 
# library(RISmed)
# disease = "ALZHEIMER DISEASE; AD"
# miRNA = "hsa-miR-659"
# my_query = paste0("(",disease,"[TIAB])"," AND ","(",miRNA,"[TIAB])")
# res <- EUtilsSummary(my_query, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2022, retmax=500)
