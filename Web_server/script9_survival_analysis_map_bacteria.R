library(dplyr)
library(survival)
library(glmnet)
library(forestplot)
library(survminer)
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)
library(pec) ##验证模型
library(rms)  ##拟合生存分析模型
library(survival)  ##生存分析包
library(glmnet)
library(pheatmap)

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}
# 设置原始文件目录
datadir <- '/data/jinchuandi/'
# datadir <- './'
# 参数读取
args <- commandArgs(T)
level = tolower(substr(args[1],1,1))
log.scale = args[2] %>% tolower()
outcomes = status = args[3]
time = paste0(status, ".time") 
Kingdom = args[4]
plot.height = args[5] %>% as.numeric()
plot.width = args[6] %>% as.numeric()

CancerType = prefix = strsplit(args[7],';') %>% unlist()
microbe.name = name = strsplit(args[8],';') %>% unlist()
# clinical = adj = strsplit(args[length(args)],';') %>% unlist() %>% tolower()
clinical = adj = ifelse(args[9]=="None",
                        args[9],
                        unlist(strsplit(args[9],';')))%>%
  tolower()

# 参数示例
# level = "g"
# log.scale = "Yes" %>% tolower()
# outcomes = status = "PFI"
# time = paste0(status, ".time")
# plot.height = 12
# plot.width = 12
# 
# CancerType = prefix = c("LUSC","LUAD")
# microbe.name = name = c("g__Pseudomonas", "g__Trueperella")
# clinical = adj = c("age")

if(Kingdom=="Bacteria"){
  if(log.scale=="yes"){
    phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all_logCPM.RData")
  }else{
    phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
  }
}else{
  if(log.scale=="yes"){
    phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all_logCPM.RData")
  }else{
    phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all.RData")
  }
}

# 检查菌是否存在
# taxall <- read.csv(paste0(datadir, "rawdata/taxaall.csv"),header = T,row.names = 1)
# if(all(name%in%c(taxall$genus,taxall$speices))){
#   print("OK")
# }else{
#   print("Not all bacteria in the list")
# }

data <- list(NULL)
for (i in 1:length(prefix)) {
  data[[i]] <- readRDS(phydata[i])
}
HRatio <- matrix(NA, nrow = length(prefix), ncol = length(name))%>%
  as.data.frame()
rownames(HRatio) <- prefix
colnames(HRatio) <- name
P.value <- matrix(NA, nrow = length(prefix), ncol = length(name))%>%
  as.data.frame()
rownames(P.value) <- prefix
colnames(P.value) <- name
for (i in 1:length(prefix)) {
  if(data[[i]]@otu_table@taxa_are_rows){
    microbe <- data[[i]]@otu_table@.Data %>%
      t() %>%
      as.data.frame()
  }else{
    microbe <- data[[i]]@otu_table@.Data %>%
      as.data.frame()
  }
  metadata <- data.frame(data[[i]]@sam_data)
  metadata$status <- metadata[,status]
  metadata$time <- metadata[,time]
  for (j in 1:length(name)) {
    if(name[j]%in%colnames(microbe)){
      dat <- data.frame(bac = microbe[,name[j]],
                        status = metadata[,status],
                        time = metadata[,time])
      if(all(adj!="none")){
        dat <- cbind(dat, metadata[adj])
      }
      fit.cox <- coxph(Surv(time, status) ~ ., data = dat)
      sum <- summary(fit.cox)
      Hratio <- c("HR" = sum$conf.int[1,1], "p.value"=sum$coefficients[1,5]) %>%
        round(3)
      HRatio[i,j] <- Hratio[1]
      P.value[i,j] <- Hratio[2]
    }else{
      HRatio[i,j] <- NA
      P.value[i,j] <- NA
    }
  }
}

HRatio_log <- log10(HRatio)
P.value.s <- ifelse(P.value<0.001, "***",
                    ifelse(P.value<0.01, "**",
                           ifelse(P.value<0.05, "*","ns")))
mycol<-colorRampPalette(c("blue","white","red"))(200)

pheat<-pheatmap(HRatio_log,scale = "none",cluster_row = TRUE, cluster_col = TRUE, border=NA,
                display_numbers = P.value.s,fontsize_number = 8, number_color = "white",
                main = "Heatmap of ln(HR)",
                color=mycol,angle_col = "315")
pdf(paste0("result/bac_survival_map_", level, "_", status, ".pdf"), height = plot.height, width = plot.width)
pheat
dev.off()
# RUN example
# Rscript script9_survival_analysis_map_bacteria.R g Yes PFI Bacteria 6 6 LUSC;LUAD g__Pseudomonas;g__Trueperella Age;Gender