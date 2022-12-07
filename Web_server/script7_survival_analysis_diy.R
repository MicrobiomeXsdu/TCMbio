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

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}
# 设置原始文件目录
datadir <- '/data/jinchuandi/'
# datadir <- './'
# 参数读取
args <- commandArgs(T)
CancerType = prefix = args[1]
level = tolower(substr(args[2],1,1))
microbe.name = name = args[3]
outcomes = status = args[4]
time = paste0(status, ".time")
Hazards.Ratio = HR = ifelse(args[5]=="Yes", TRUE,FALSE)
Confidence.Interval = CI = ifelse(args[6]=="Yes", TRUE,FALSE)
group.cutoff = ifelse(args[7]=="Median", 0.5, as.numeric(args[7]))
Kingdom = args[8]
plot.height = args[9] %>% as.numeric()
plot.width = args[10] %>% as.numeric()
# clinical = adj = args[10:length(args)] %>% tolower()
# clinical = adj = strsplit(args[length(args)],';') %>% unlist() %>% tolower()
clinical = adj = ifelse(args[length(args)]=="None", 
                        args[length(args)], 
                        unlist(strsplit(args[length(args)],';')))%>%
  tolower()
# print(clinical)

# # 参数示例
# CancerType = prefix = "LUSC"
# level = "g"
# microbe.name = name = "g__Pseudomonas"
# outcomes = status = "PFI"
# time = paste0(status, ".time")
# Hazards.Ratio = HR = TRUE
# Confidence.Interval = CI = TRUE
# group.cutoff = ifelse("Median"=="Median", 0.5, as.numeric(args[7]))
# plot.height = 6 %>% as.numeric()
# plot.width = 6 %>% as.numeric()
# clinical = adj = c("age","gender", "race")
# phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
# Kingdom = 'Bacteria'

if(Kingdom=="Bacteria"){
  phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all_logCPM.RData")
}else{
  phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all_logCPM.RData")
}

# 检查菌是否存在
taxall <- read.csv(paste0(datadir, "rawdata/taxaall.csv"),header = T)
if(name%in%c(taxall$genus,taxall$speices)){
  print("OK")
}else{
  print("No such bacteria")
}


# 读取数据并绘图

data <- readRDS(phydata)
if(data@otu_table@taxa_are_rows){
  microbe <- data@otu_table@.Data %>%
    t() %>%
    as.data.frame()
}else{
  microbe <- data@otu_table@.Data %>%
    as.data.frame()
}
metadata <- data.frame(data@sam_data)
metadata$status <- metadata[,status]
metadata$time <- metadata[,time]

if(name%in%colnames(microbe)){
  data <- data.frame(bac = microbe[,name],
                     status = metadata[,status],
                     time = metadata[,time])
  data$group <- ifelse(data$bac>quantile(data$bac,group.cutoff), paste0(name, "-High"), paste0(name, "-Low")) %>%
    factor(c(paste0(name, "-High"), paste0(name, "-Low")))
  # data <- cbind(data, metadata[adj])
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
      labs(title = paste0(status, "\n HR: ", Hratio[1], "(", Hratio[2],", ",Hratio[3],")"))+
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
  try(dev.off(), silent = TRUE)
  pdf(paste0("result/",prefix, "_", name, "_", status, ".pdf"),width=plot.width,height=plot.height)
  print(p, newpage = FALSE)
  dev.off()
}else{
  print("No such bacteria")
}

# RUN example
# Rscript script7_survival_analysis_diy.R BRCA s s__Candida_albicans PFI Yes Yes Median Fungus 6 6 Age;Gender;Race
