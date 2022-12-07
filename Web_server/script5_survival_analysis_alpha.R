library(rms)
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
library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
# 原始数据
datadir <- '/data/jinchuandi/'
# datadir <- './'
# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}

# 参数读取
args <- commandArgs(T)
CancerType = prefix = args[1]
Alpha.Index = name = args[2]
level = tolower(substr(args[3],1,1))
outcomes = status = args[4]
time = paste0(status, ".time")
Confidence.Interval = CI = ifelse(args[5]=="Yes", TRUE,FALSE)
Hazards.Ratio = HR = ifelse(args[6]=="Yes", TRUE,FALSE)
group.cutoff = ifelse(args[7]=="Median", 0.5, as.numeric(args[7]))
Kingdom = args[8]
plot.height = args[9] %>% as.numeric()
plot.width = args[10] %>% as.numeric()
# clinical = adj = args[10:length(args)] %>% tolower()
clinical = adj = strsplit(args[length(args)],';') %>% unlist() %>% tolower()


# # 参数示例
# CancerType = prefix = "BRCA"
# Alpha.Index = name = "Shannon"
# level = tolower(substr("Species",1,1))
# outcomes = status = "OS"
# time = paste0(status, ".time")
# Confidence.Interval = CI = ifelse("Yes"=="Yes", TRUE,FALSE)
# Hazards.Ratio = HR = ifelse("Yes"=="Yes", TRUE,FALSE)
# group.cutoff = ifelse("Median"=="Median", 50, as.numeric(args[7]))
# Kingdom = "Fungus"
# plot.height = "6" %>% as.numeric()
# plot.width = "6" %>% as.numeric()
# # clinical = adj = args[10:length(args)] %>% tolower()
# clinical = adj = strsplit("None",';') %>% unlist() %>% tolower()

if(Kingdom=="Bacteria"){
  phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
}else{
  phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all.RData")
}
# 读取数据并绘图

data <- readRDS(phydata)

if(data@otu_table@taxa_are_rows==TRUE){
  otu <- data@otu_table@.Data %>%
    t()
}else{
  otu <- data@otu_table@.Data
}

alpha_res <- data.frame(row.names = rownames(otu), 
                        Observed = estimateR(otu)[1, ], 
                        Chao1 = estimateR(otu)[2, ], 
                        Shannon = diversity(otu, index = 'shannon', base = exp(1)), 
                        Simpson = diversity(otu, index = 'simpson'))
alpha_res$alpha_index <- alpha_res[,name]
alpha_res <- na.omit(alpha_res)

metadata <- data.frame(data@sam_data)
metadata$status <- metadata[,status]
metadata$time <- metadata[,time]

data <- data.frame(alpha = alpha_res$alpha_index,
                   status = metadata[,status],
                   time = metadata[,time])
data$group <- ifelse(data$alpha>quantile(data$alpha,group.cutoff), paste0(name, "-High"), paste0(name, "-Low")) %>%
  factor(levels = c(paste0(name, "-High"), paste0(name, "-Low")))

if(all(adj!="none")){
  data <- cbind(data, metadata[adj])
}


fit <- survfit(Surv(time, status) ~ group,  # 创建生存对象 
               data = data) # 数据集来源
fit # 查看拟合曲线信息
summary(fit)
if(HR){
  fit.cox <- coxph(Surv(time, status) ~ ., data = data[-1])
  sum <- summary(fit.cox)
  Hratio <- c("HR" = sum$conf.int[1,1], "HR.up"=sum$conf.int[1,4], "HR.low"=sum$conf.int[1,3]) %>%
    round(2)
  p <- ggsurvplot(fit, data = data,
                  pval = TRUE,
                  surv.median.line = "hv", # 添加中位生存时间线
                  conf.int = CI
  )+
    labs(title = paste0(status, "\n HR: ", Hratio[1], "(", Hratio[3],", ",Hratio[2],")"))+
    xlab("Time (Day)")
}else{
  p <- ggsurvplot(fit, data = data,
                  pval = TRUE,
                  surv.median.line = "hv", # 添加中位生存时间线
                  conf.int = CI
  )+
    labs(title = status)+
    xlab("Time (Day)")
}
pdf(paste0( "result/",prefix, "_", name, "_", level, "_", status, ".pdf"),width=plot.width,height=plot.height)
print(p, newpage = FALSE)
dev.off()
## RUN example
# Rscript script5_survival_analysis_alpha.R BRCA Simpson s PFI Yes Yes Median Fungus 6 6 Age;Gender;Race
# Rscript script5_survival_analysis_alpha.R BRCA Simpson s PFI Yes Yes Median Fungus 6 6 None
# Rscript script5_survival_analysis_alpha.R BRCA Simpson s PFI Yes Yes 0.75 Fungus 6 6 None
