library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)

# 参数读取
args <- commandArgs(T)
group = args[1]
if(group!="MultipleCancersComparison"){
  CancerType = prefix = args[2]
  alpha_index=args[3]
  level = tolower(substr(args[4],1,1))
  Kingdom = args[5]
  plot.height = args[6] %>% as.numeric()
  plot.width = args[7] %>% as.numeric()
}else{
  CancerTypeA = prefix = args[2]
  CancerTypeB = prefix2 = args[3]
  alpha_index=args[4]
  level = tolower(substr(args[5],1,1))
  Kingdom = args[6]
  plot.height = args[7] %>% as.numeric()
  plot.width = args[8] %>% as.numeric()
}

# ## 参数示例
# group = "SampleType"
# if(group!="MultipleCancersComparison"){
#   CancerType = prefix = "BRCA"
#   alpha_index="Shannon"
#   level = tolower(substr("Species",1,1))
#   Kingdom = "Fungus"
#   plot.height = "6" %>% as.numeric()
#   plot.width = "7" %>% as.numeric()
# }else{
#   CancerTypeA = prefix = "BRCA"
#   CancerTypeB = prefix2 = "CESC"
#   alpha_index="Shannon"
#   level = tolower(substr("Species",1,1))
#   Kingdom = "Fungus"
#   plot.height = "6" %>% as.numeric()
#   plot.width = "6" %>% as.numeric()
# }

datadir <- '/data/jinchuandi/'
#datadir <- './'
# 
if(Kingdom=="Bacteria"){
  if(group=="SampleType"){
    phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_tumorvsnormal.RData")
  }else{
    if(group=="MultipleCancersComparison"){
      phydata1 = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
      phydata2 = paste0(datadir, "rawdata/", prefix2, "_", level, "_all.RData")
    }else{
      phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all.RData")
    }
  }
}else{
  if(group=="SampleType"){
    phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_tumorvsnormal.RData")
  }else{
    if(group=="MultipleCancersComparison"){
      phydata1 = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all.RData")
      phydata2 = paste0(datadir, "rawdata_fungi/", prefix2, "_", level, "_all.RData")
    }else{
      phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all.RData")
    }
  }
  
}


# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}

if(group=="MultipleCancersComparison"){
  phy_dat1 <- readRDS(phydata1)
  phy_dat1@sam_data$CancerType <- prefix
  phy_dat2 <- readRDS(phydata2)
  phy_dat2@sam_data$CancerType <- prefix2
  otu1 <- as.data.frame(phy_dat1@otu_table@.Data) %>%
    tibble::rownames_to_column("microbe")
  otu2 <- as.data.frame(phy_dat2@otu_table@.Data) %>%
    tibble::rownames_to_column("microbe")
  otu <- merge(otu1, otu2, all = T, by = "microbe")
  otu[is.na(otu)] <- 0
  otu <- otu %>%
    tibble::column_to_rownames("microbe") %>%
    as.matrix()
  otu_tab <- otu_table(otu, taxa_are_rows = T)
  metadata <- data.frame(row.names = c(sample_names(phy_dat1), sample_names(phy_dat2)),
                         CancerType = c(phy_dat1@sam_data$CancerType, phy_dat2@sam_data$CancerType))
  metadata <- sample_data(metadata)
  phy_dat <- phyloseq(otu_tab, metadata)
}else{
  phy_dat <- readRDS(phydata)
}

if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
# try(names(phy_dat@sam_data)[which(names(phy_dat@sam_data)=="SampleType")] <- "TumorvsNon.Tumor", silent = T)
try(names(phy_dat@sam_data)[which(names(phy_dat@sam_data)=="CancerType")] <- "MultipleCancersComparison", silent = T)
phy_dat <- prune_samples(!is.na(phy_dat@sam_data[[group]]), phy_dat)
metadata <- data.frame(phy_dat@sam_data)
alpha_res <- data.frame(row.names = sample_names(phy_dat), Observed = rep(NA, length(sample_names(phy_dat))), Chao1 = rep(NA, length(sample_names(phy_dat))), Shannon = rep(NA, length(sample_names(phy_dat))), Simpson =rep(NA, length(sample_names(phy_dat))),group = metadata[,group])
otu <- otu[sample_names(phy_dat),]
if(max(length(unique(alpha_res$group))<2|is.na(unique(alpha_res$group)))==1){
  stop("The number of groups is less than two")
}else{
  alpha_res$Observed <- rowSums(otu > 0)
  # 计算 Chao 1指数
  alpha_res$Chao1 <- estimateR(otu)[2, ]
  # 计算Shannon 指数,通常使用e作为底数
  alpha_res$Shannon <- diversity(otu, index = 'shannon', base = exp(1))
  # Simpson指数分为经典 Simpson 指数和Gini-Simpson 指数，不过平时常用的 Simpson 指数即为 Gini-Simpson 指数
  # Gini-Simpson 指数代码
  alpha_res$Simpson <- diversity(otu, index = 'simpson')
  alpha_res$alpha_index <- alpha_res[,alpha_index]
  alpha_res <- na.omit(alpha_res)
  alpha_res$group <- paste0(toupper(substr(alpha_res$group, 1,1)), substr(alpha_res$group, 2,nchar(alpha_res$group)))
  combine <- combn(unique(alpha_res$group),2)
  list <- list(NULL)
  for (i in 1:ncol(combine)) {
    list[[i]] <- combine[,i]
  }
  alpha_boxplot <- ggplot2::ggplot(alpha_res, aes(x=group, y=alpha_index, fill=group))+
    geom_boxplot()+labs(x=NULL,y=alpha_index)+
    guides(fill=guide_legend(group)) +
    theme(plot.title=element_text(hjust=0.5))+
    stat_compare_means(comparisons = list,
                       method = "wilcox",
                       label = "p.signif") +
    theme_classic()
  if (group=="MultipleCancersComparison") {
    ggsave(alpha_boxplot,filename=paste0("result/",prefix, "_", prefix2, "_", alpha_index, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
  }else{
    ggsave(alpha_boxplot,filename=paste0("result/",prefix, "_", alpha_index, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
  }
}
# RUN example
# Rscript script1_alphadiversity.R SampleType BRCA Shannon s Fungus 6 6
