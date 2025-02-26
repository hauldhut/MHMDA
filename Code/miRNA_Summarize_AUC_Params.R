library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)

# Methods = c("M1", "M2", "M3", "M12")
# Methods = c("H1","H2","MH")
Method = "M12"
BipartiteNet = "miR2Disease"  #miR2Disease/HMDD

setwd("~/Manuscripts/100MHMDA/Code")

Result_bySLFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_Score_n_Label.csv")
df.ScoreLabel = read.csv(Result_bySLFile)
Scores = df.ScoreLabel$Scores
Labels = df.ScoreLabel$Labels

resultspred = prediction(Scores, Labels)
auc.perf = performance(resultspred, measure = "auc")
aucavgbyAll = round(auc.perf@y.values[[1]],3)

# roc.perf <- performance(resultspred,"tpr","fpr")
# plot(roc.perf@x.values[[1]], roc.perf@y.values[[1]], xlab="FPR", ylab="TPR", main = as.character(aucavgbyAll), sub = "Hello")


Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_ROC.csv")
df.AUCbyTrial = read.csv(Result_byTrialFile)

aucavgbyTrial = round(mean(df.AUCbyTrial$auc),3)
aucavgbyTrial.sd = round(sd(df.AUCbyTrial$auc),3)

# Result_byDiseaseFile = paste0("../Results/",Method,"_byDisease_",BipartiteNet,"_ROC.csv")
# df.AUCbyDisease = read.csv(Result_byDiseaseFile)

# df.AUCbyDisease <- aggregate(auc~disease, data=df.AUCbyTrial, FUN=function(x) c(mean=mean(x), count=length(x)))
# 
# aucavgbyDisease = round(mean(as.numeric(df.AUCbyDisease$auc[,1])),3)
# aucavgbyDisease.sd = round(sd(as.numeric(df.AUCbyDisease$auc[,1])),3)
# cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")","aucavgbyDisease=",aucavgbyDisease,"(+-",aucavgbyDisease.sd,")\n")

cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")

####################
library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
library('cowplot')
setwd("~/Manuscripts/100Net4DMP2/Code")


Methods = c("M12", "MH")
BipartiteNets = c("miR2Disease")
BipartiteNet = "miR2Disease"

Params = c("gamma","delta")#,"delta"
h.config2dfAUC = hash()
for(Param in Params){
  config = paste0(BipartiteNet,"_",Param)
  df.AUC = NULL
  for(Method in Methods){
    for(pv in c(0.1, 0.3, 0.5, 0.7, 0.9)){
      cat("Analyzing for", Param, Method, "\n")
      
      if (Param=="gamma"){
        Result_bySLFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",pv,"_delta_0.5_Score_n_Label.csv")  
      }else{
        Result_bySLFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_0.5_delta_",pv,"_Score_n_Label.csv")  
      }
      
      print(Result_bySLFile)
      df.ScoreLabel = read.csv(Result_bySLFile)
      Scores = df.ScoreLabel$Scores
      Labels = df.ScoreLabel$Labels
      resultspred = prediction(Scores, Labels)
      roc.perf <- performance(resultspred,"tpr","fpr")
      FPR = roc.perf@x.values[[1]]
      TPR = roc.perf@y.values[[1]]
      
      auc.perf = performance(resultspred, measure = "auc")
      aucavgbyAll = round(auc.perf@y.values[[1]],3)
      
      if(Method == "MH"){
        NetworkType = paste0("DiSimNet-MultiNet_miRNA")
      }else{
        NetworkType = "MultiNet_miRNA"  
      }
      df.AUC = rbind(df.AUC, data.frame(param = pv, AU = aucavgbyAll, NetworkType = NetworkType))
    }
  }
  #Store
  SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",config,".rdata")
  saveRDS(df.AUC,SummaryFile)#readRDS
  
  h.config2dfAUC[[config]] = df.AUC
}
length(h.config2dfAUC)

h.config2dfAUC[["miR2Disease_gamma"]]
h.config2dfAUC[["miR2Disease_delta"]]

# # #Load:
# Methods = c("M123", "MH123")
# 
# BipartiteNet = "PREDICT"
# 
# Params = c("r","delta")
# Param = "r"
# config = paste0(BipartiteNet,"_",Param)
# 
# setwd("~/Manuscripts/99HDR2/Code")
# DrugSimNetsSTR = paste0(BipartiteNets,collapse = "-")
# SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",config,".rdata")
# 
# df.AUC = readRDS(SummaryFile)
# h.config2dfAUC = hash()
# h.config2dfAUC[[config]] = df.AUC

######


#Draw
h.config2AUplot = hash()
for(Param in Params){
  for(BipartiteNet in BipartiteNets){
    config = paste0(BipartiteNet,"_",Param)
    df.AU = h.config2dfAUC[[config]]
    
    df.AU$param = as.factor(df.AU$param)
    
    if (Param=="gamma"){
      ParamSymbol = expression(paste("restart-probability (",gamma,")"))
    }else{
      ParamSymbol = bquote("jumping-probability ("~delta~")")
    }
    
    fontsize = 16
    p = ggplot(data=df.AU, aes(x=param, y=AU, group = NetworkType)) +
      geom_line(stat="identity",position=position_dodge(),aes(color=NetworkType)) +
      geom_point(aes(color=NetworkType)) +
      ylim(0,1)
    
    # p = p + labs(x=paste0("\n",ParamSymbol), y = "AUC\n") +
    p = p + labs(x=ParamSymbol, y = "AUC\n") +
      scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
      # scale_fill_brewer(palette="Blues") +
      theme_light() +
      theme(
        text = element_text(size=fontsize),# All font sizes
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=fontsize),
        legend.title=element_blank(),#Remove legend title (Network)
        legend.position = c(0.5, 0.11),
        axis.text = element_text(size = fontsize),
        axis.title = element_text(size = fontsize),
        axis.title.x = element_text(margin = margin(t = 20))
      )
    
    h.config2AUplot[[config]] = p
  }
}

library('cowplot')
plot_grid(h.config2AUplot[["miR2Disease_gamma"]],h.config2AUplot[["miR2Disease_delta"]], labels=c("(a)", "(b)"), label_size = fontsize+2, ncol = 2, nrow = 1)

FigureFile = paste0("../Figures/Figure_Parameter.pdf")
ggsave(FigureFile, width = 14, height = 7)
# embed_fonts(FigureFile)
