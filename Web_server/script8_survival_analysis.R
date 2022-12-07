library(rms)
library(dplyr)
library(survival)
library(glmnet)
library(forestplot)
library(survminer)
library(ROCR)  
library(glmnet) 
library(caret)
library(pec)
library(rms)
library(survival)
library(glmnet)

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
outcomes = status = args[3]
time = paste0(status, ".time")
log.scale = args[4] %>% tolower()
P.value.Cutoff = args[5] %>% as.numeric()
Adjust.Method= padjust = tolower(args[6])
Kingdom = args[7]
plot.height = args[8] %>% as.numeric()
plot.width = args[9] %>% as.numeric()
# clinical = adj = strsplit(args[length(args)],';') %>% unlist() %>% tolower()
clinical = adj = ifelse(args[10]=="None",
                        args[10],
                        unlist(strsplit(args[10],';')))%>%
  tolower()

# 参数示例
# CancerType = prefix = "LUSC"
# level = tolower(substr("Species",1,1))
# outcomes = status = "OS"
# time = paste0(status, ".time")
# log.scale = "Yes" %>% tolower()
# P.value.Cutoff = 0.05 %>% as.numeric()
# Adjust.Method= padjust = tolower("None")
# plot.height = args[7] %>% as.numeric()
# plot.width = args[8] %>% as.numeric()
# # clinical = adj = strsplit(args[length(args)],';') %>% unlist() %>% tolower()
# clinical = adj = ifelse(args[length(args)]=="None",
#                         args[length(args)],
#                         unlist(strsplit(args[length(args)],';')))%>%
#   tolower()
# Kingdom = 'Bacteria'

# phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
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
if(all(adj!="none")){
  metadata <- metadata[,c("status", "time", adj)] %>%
    na.omit()
}else{
  metadata <- metadata[,c("status", "time")]%>%
    na.omit()
}
microbe <- microbe[rownames(metadata),]
spe_res_surv <- data.frame(row.names = colnames(microbe),
                           coef = rep(NA, ncol(microbe)), 
                           HR = rep(NA, ncol(microbe)), 
                           HR.up = rep(NA, ncol(microbe)),
                           HR.low = rep(NA, ncol(microbe)),
                           p.value = rep(NA, ncol(microbe)))

for (i in colnames(microbe)) {
  metadata$bac <- microbe[,i]
  fit <- NULL
  try(fit <- coxph(Surv(metadata$time, metadata$status) ~ ., data = metadata),silent = T)
  if(!is.null(fit)){
    sum <- summary(fit)
    spe_res_surv[i,"coef"] <- sum$coefficients[nrow(sum$coefficients),1]
    spe_res_surv[i,"HR"] <- sum$conf.int[nrow(sum$coefficients),1]
    spe_res_surv[i,"HR.up"] <- sum$conf.int[nrow(sum$coefficients),4]
    spe_res_surv[i,"HR.low"] <- sum$conf.int[nrow(sum$coefficients),3]
    spe_res_surv[i,"p.value"] <- sum$coefficients[nrow(sum$coefficients),5]
  }
  
}
spe_res_surv <- na.omit(spe_res_surv)
spe_res_surv$p.adj <- p.adjust(spe_res_surv$p.value, padjust)
spe_res_surv[spe_res_surv=="Inf"] <- NA
spe_res_surv <- na.omit(spe_res_surv)
spe_surv <- spe_res_surv[spe_res_surv$p.value<P.value.Cutoff & !is.na(spe_res_surv$HR.up),] %>%tibble::rownames_to_column("variable")
spe_surv_surv <- microbe[,colnames(microbe)%in%c(spe_surv$variable)]

names(spe_surv) <- c("Bacteria", "Coefficient", "HR", "HR.up", "HR.low", "P-value", "Adjusted P-value")

write.csv(spe_surv, paste0("result/", prefix, "_spe_cox_coef_", level, "_", status, ".csv"), row.names = F)

##火山图1
library("ggplot2")
#自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
spe_res_surv$sig <- NA
spe_res_surv[which(spe_res_surv$p.adj <= P.value.Cutoff & spe_res_surv$HR < 1),'sig'] <- 'Protective'
spe_res_surv[which(spe_res_surv$p.adj <= P.value.Cutoff & spe_res_surv$HR > 1),'sig'] <- 'Risk'
spe_res_surv[which(spe_res_surv$p.adj > P.value.Cutoff ),'sig'] <- 'None'
spe_res_surv$logp <- -log10(spe_res_surv$p.adj)

p <- ggplot(spe_res_surv, aes(x = HR, y = logp, color = sig)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Risk', 'Protective', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 1, color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  # xlim(0, 3) +
  labs(x = '\n HR', y = '-Log10 (P Value)\n', color = '', title= 'Volcano Plot \n')
# p
ggsave(p, filename=paste0("result/", prefix, "_survival_volcano_", level, "_", status, ".pdf"), width = plot.width, height = plot.height, units = c("in"))

# RUN example
# Rscript script8_survival_analysis.R BRCA s PFI Yes 0.05 None Fungus 6 6 None