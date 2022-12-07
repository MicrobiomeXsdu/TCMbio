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
library(vegan)

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
outcomes = status = args[2]
time = paste0(status, ".time") 
Kingdom = args[3]
plot.height = args[4] %>% as.numeric()
plot.width = args[5] %>% as.numeric()
CancerType = prefix = strsplit(args[6],';') %>% unlist()
# Alpha.Index = name = args[args%in%c("Simpson", "Shannon", "Chao1", "Observed")]
Alpha.Index = name = strsplit(args[7],';') %>% unlist()
# clinical = adj = args[args%in%c("Age", "Gender", "Race", "None")] %>% tolower()
clinical = adj = ifelse(args[8]=="None", 
                        args[8], 
                        unlist(strsplit(args[8],';')))%>%
  tolower()
# print(paste0("CancerType = ", CancerType,", Alpha.Index = ",Alpha.Index,", clinical = ", clinical))
# 参数示例
# level = "s"
# outcomes = status = "OS"
# time = paste0(status, ".time")
# Kingdom = "Fungus"
# plot.height = 6
# plot.width = 6
# CancerType = prefix = c("BRCA","CESC")
# Alpha.Index = name = c("Simpson", "Shannon")
# clinical = adj = ifelse("None"=="None", 
#                         "None", 
#                         unlist(strsplit("None",';')))%>%
#   tolower()

if(Kingdom=="Bacteria"){
  phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
}else{
  phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all.RData")
}

data <- list(NULL)
for (i in 1:length(prefix)) {
  data[[i]] <- readRDS(phydata[i])
}

alpha <- list(NULL)
for (i in 1:length(prefix)) {
  if(data[[i]]@otu_table@taxa_are_rows){
    otu <- data[[i]]@otu_table@.Data %>%
      t()
  }else{
    otu <- data[[i]]@otu_table@.Data
  }
  alpha[[i]] <- data.frame(row.names = rownames(otu), 
                           Observed = estimateR(otu)[1, ], 
                           Chao1 = estimateR(otu)[2, ], 
                           Shannon = diversity(otu, index = 'shannon', base = exp(1)), 
                           Simpson = diversity(otu, index = 'simpson'))
  alpha[[i]] <- na.omit(alpha[[i]])
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
  alpha_res <- alpha[[i]]
  metadata <- data.frame(data[[i]]@sam_data)
  metadata$status <- metadata[,status]
  metadata$time <- metadata[,time]
  for (j in 1:length(name)) {
    dat <- data.frame(alpha = alpha_res[,name[j]],
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
                color=mycol,angle_col = "315"
)

pdf(paste0("result/survival_alpha_map_", level, "_", status, ".pdf"), height = plot.height, width = plot.width)
pheat
dev.off()
# RUN example
# Rscript script6_survival_analysis_map_alpha.R s OS Fungus 6 6 CESC;BRCA Simpson;Shannon Age;Gender
