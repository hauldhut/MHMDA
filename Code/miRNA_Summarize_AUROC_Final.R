library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
Method = "M12"#MH/H1/H2/M12/M1/M2

BipartiteNet = "HMDD"#HMDD/miR2Disease

setwd("~/Manuscripts/100MHMDA/Code")

Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_Score_n_Label.csv")
df.ScoreLabel = read.csv(Result_by_Score_n_LabelFile)
Scores = df.ScoreLabel$Scores
Labels = df.ScoreLabel$Labels

resultspred = prediction(Scores, Labels)
auc.perf = performance(resultspred, measure = "auc")
aucavgbyAll = round(auc.perf@y.values[[1]],3)

roc.perf <- performance(resultspred,"tpr","fpr")
# plot(roc.perf@x.values[[1]], roc.perf@y.values[[1]], xlab="FPR", ylab="TPR", main = as.character(aucavgbyAll), sub = "Hello")


Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_ROC.csv")
df.AUCbyTrial = read.csv(Result_byTrialFile)

aucavgbyTrial = round(mean(df.AUCbyTrial$auc),3)
aucavgbyTrial.sd = round(sd(df.AUCbyTrial$auc),3)

# Result_byDiseaseFile = paste0("../Results/",Method,"_byDisease_",BipartiteNet,"_ROC.csv")
# df.AUCbyDisease = read.csv(Result_byDiseaseFile)

df.AUCbyDisease <- aggregate(auc~disease, data=df.AUCbyTrial, FUN=function(x) c(mean=mean(x), count=length(x)))

aucavgbyDisease = round(mean(as.numeric(df.AUCbyDisease$auc[,1])),3)
aucavgbyDisease.sd = round(sd(as.numeric(df.AUCbyDisease$auc[,1])),3)


cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")","aucavgbyDisease=",aucavgbyDisease,"(+-",aucavgbyDisease.sd,")\n")

####################
library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
library('cowplot')
setwd("~/Manuscripts/100MHMDA/Code")
BipartiteNets = c("miR2Disease","HMDD")
h.All = hash()
Methods = c("MH","H1","H2")
# Methods = c("MH")
# Methods = c("M12","M1","M2","M3")
# Methods = c("M12")

for(BipartiteNet in BipartiteNets){
  cat("Calculating AUC for ", BipartiteNet, "\n")
  
  df.All = NULL
  
  for(Method in Methods){
    cat("\t==> Analyzing for ", Method, "\n")
    Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_Score_n_Label.csv")
    # Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_0.5_delta_0.5","_Score_n_Label.csv")
    df.ScoreLabel = read.csv(Result_by_Score_n_LabelFile)
    Scores = df.ScoreLabel$Scores
    Labels = df.ScoreLabel$Labels
    resultspred = prediction(Scores, Labels)
    roc.perf <- performance(resultspred,"tpr","fpr")
    FPR = roc.perf@x.values[[1]]
    TPR = roc.perf@y.values[[1]]
    
    auc.perf = performance(resultspred, measure = "auc")
    aucavgbyAll = round(auc.perf@y.values[[1]],3)
    
    
    if(Method=="MH"){
      Network = rep(paste0("MHMDA-MH (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="H1"){
      Network = rep(paste0("RWRHMDA (MonoNet_miRWalk) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="H2"){
      Network = rep(paste0("RWRHMDA (MonoNet_TargetScan) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="H3"){
      Network = rep(paste0("RWRHMDA (MonoNet_Integrated) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="M12"){
      Network = rep(paste0("MHMDA-M (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="M1"){
      Network = rep(paste0("RWRMDA (MonoNet_miRWalk) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else if(Method=="M2"){
      Network = rep(paste0("RWRMDA (MonoNet_TargetScan) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }else{
      Network = rep(paste0("RWRMDA (MonoNet_Integrated) (AUROC = ",aucavgbyAll,")"),length(FPR))  
    }
    df.All = rbind(df.All, data.frame(FPR = FPR, TPR = TPR, Network = Network))
  }
  h.All[[BipartiteNet]] = df.All
}
length(h.All)

######
FigureFile = paste0("../Figures/Figure_",paste0(Methods,collapse = "-"))

fontsize = 12
l.plot = list()

for(BipartiteNet in keys(h.All)){
  print(paste0("Processing plot for ", BipartiteNet))
  df.O = h.All[[BipartiteNet]]
  if(nrow(df.O)>=10^6){
    df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/100),]  
  }else{
    df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/10),]
  }
  
  p <- ggplot(df.D, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
    # geom_point() +
    geom_line(size=1) +
    ggtitle(BipartiteNet) + #Set plot title
    xlab("\nFalse Positive Rate (FPR)") + #Set x label
    ylab("True Positive Rate (TPR)\n") + #Set y label
    # scale_color_manual(labels = c("H1", "H2", "MH"), values = c("red","green","blue")) +
    theme_light() +
    theme(
      text = element_text(size=fontsize),# All font sizes
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size=fontsize-1),
      legend.position = c(0.55, 0.15),
      legend.title=element_blank(),#Remove legend title (Network)
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize)
      # plot.background = element_blank(),
      # panel.background = element_blank()
    )
  l.plot[[BipartiteNet]] = p
  saveRDS(p,paste0(FigureFile,"_",BipartiteNet,".rdata"))#readRDS
}
library('cowplot')
plot_grid(l.plot[["HMDD"]], l.plot[["miR2Disease"]], labels=c("(a)", "(b)"), ncol = 2, nrow = 1)

ggsave(paste0(FigureFile,"_AUROC.pdf"), width = 11, height = 5.5)


# ####################################
# #Compare the Chart between Origin (O) and Downsampled (D) data
# df.O = h.All[["miR2Disease"]]#miR2Disease/HMDD
# p.O = ggplot(df.O, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
#   # geom_point() +
#   geom_line(size=1) +
#   theme(
#     legend.position = c(0.5, 0.15)
#   )
# 
# # #Sol 1: Using downsample in groupdata2, or downSample in caret. Both down sampling to the number of sample in minority class (Network)
# # library(groupdata2)
# # df.D = downsample(df.O, cat_col = "Network")
# # # library(caret)
# # # df.D = downsample(df.O[,1:2], as.factor(df.O$Network))
# # 
# # p.D = ggplot(df.D, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
# #   geom_point() +
# #   theme(
# #     legend.position = c(0.5, 0.15)
# #   )
# 
# #Sol 2: Use sample() function. More plexible
# df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/10),]
# 
# p.D = ggplot(df.D, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
#   # geom_point() +
#   geom_line(size=1) +
#   theme(
#     legend.position = c(0.5, 0.15)
#   )
# 
# plot_grid(p.O, p.D, labels=c("A", "B"), ncol = 2, nrow = 1)
# 
