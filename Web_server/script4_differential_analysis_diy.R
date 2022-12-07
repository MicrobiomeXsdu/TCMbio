library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}
# 设置原始文件目录
datadir <- '/data/jinchuandi/'
# datadir <- './'
# 参数读取
args <- commandArgs(T)
group = args[1]
if(group!="MultipleCancersComparison"){
  CancerType = prefix = args[2]
  level = tolower(substr(args[3],1,1))
  microbe.name = args[4]
  Kingdom = args[5]
  log.scale = args[6] %>% tolower()
  plot.height = args[7] %>% as.numeric()
  plot.width = args[8] %>% as.numeric()
}else{
  CancerTypeA = prefix = args[2]
  CancerTypeB = prefix2 = args[3]
  level = tolower(substr(args[4],1,1))
  microbe.name = args[5] 
  Kingdom = args[6]
  log.scale = args[7] %>% tolower()
  plot.height = args[8] %>% as.numeric()
  plot.width = args[9] %>% as.numeric()
}
# # # 参数示例
# group = "SampleType"
# if(group!="MultipleCancersComparison"){
#   CancerType = prefix = "CESC"
#   level = tolower(substr("Species",1,1))
#   microbe.name = "s__Cercospora_beticola"
#   Kingdom = "Fungus"
#   log.scale = "Yes" %>% tolower()
#   plot.height = "6" %>% as.numeric()
#   plot.width = "6" %>% as.numeric()
# }else{
#   CancerTypeA = prefix = "CESC"
#   CancerTypeB = prefix2 = "UCEC"
#   level = tolower(substr("Species",1,1))
#   microbe.name = "s__Cercospora_beticola"
#   Kingdom = "Fungus"
#   log.scale = "Yes" %>% tolower()
#   plot.height = "6" %>% as.numeric()
#   plot.width = "6" %>% as.numeric()
# }

if(Kingdom=="Bacteria"){
  if(log.scale=="no"){
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
      phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_tumorvsnormal_logCPM.RData")
    }else{
      if(group=="MultipleCancersComparison"){
        phydata1 = paste0(datadir, "rawdata/", prefix, "_", level, "_all_logCPM.RData")
        phydata2 = paste0(datadir, "rawdata/", prefix2, "_", level, "_all_logCPM.RData")
      }else{
        phydata = paste0(datadir, "rawdata/", prefix, "_", level, "_all_logCPM.RData")
      }
    }
  }
}else{
  if(log.scale=="no"){
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
  }else{
    if(group=="SampleType"){
      phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_tumorvsnormal_logCPM.RData")
    }else{
      if(group=="MultipleCancersComparison"){
        phydata1 = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all_logCPM.RData")
        phydata2 = paste0(datadir, "rawdata_fungi/", prefix2, "_", level, "_all_logCPM.RData")
      }else{
        phydata = paste0(datadir, "rawdata_fungi/", prefix, "_", level, "_all_logCPM.RData")
      }
    }
  }
}

# 检查菌是否存在
taxall <- read.csv(paste0(datadir, "rawdata/taxaall.csv"),header = T)
if(microbe.name%in%c(taxall$genus,taxall$speices)){
  print("OK")
}else{
  print("No such bacteria")
}

# 读取数据

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
                         CancerType = c(phy_dat1@sam_data$CancerType, phy_dat2@sam_data$CancerType)) %>%
    sample_data()
  try(names(metadata)[which(names(metadata)=="CancerType")] <- "MultipleCancersComparison", silent = T)
  phy_dat <- phyloseq(otu_tab, metadata)
}else{
  phy_dat <- readRDS(phydata)
  # try(names(phy_dat@sam_data)[which(names(phy_dat@sam_data)=="SampleType")] <- "TumorvsNon.Tumor", silent = T)
}
try(phy_dat <- prune_samples(!is.na(phy_dat@sam_data[[group]]), phy_dat),silent = T)
if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
metadata <- data.frame(phy_dat@sam_data)

if(microbe.name%in%colnames(otu)){
  res <- data.frame(row.names = sample_names(phy_dat), bac = otu[,microbe.name],group = metadata[,group])
  if(max(length(unique(res$group))<2|is.na(unique(res$group)))==1){
    print("The number of groups is less than two")
  }else{
    res <- na.omit(res)
    
    res$group <- paste0(toupper(substr(res$group, 1,1)), substr(res$group, 2,nchar(res$group)))
    combine <- combn(unique(res$group),2)
    list <- list(NULL)
    for (i in 1:ncol(combine)) {
      list[[i]] <- combine[,i]
    }
    boxplot <- ggplot2::ggplot(res, aes(x=group, y=bac, fill=group))+
      geom_boxplot()+labs(x=group, y=microbe.name)+
      guides(fill=guide_legend(group)) +
      theme(plot.title=element_text(hjust=0.5))+
      stat_compare_means(comparisons = list,
                         method = "wilcox",
                         label = "p.signif") +
      theme_classic()
    if (group=="MultipleCancersComparison") {
      ggsave(boxplot,filename=paste0("result/",prefix, "_", prefix2, "_", microbe.name, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
    }else{
      ggsave(boxplot,filename=paste0("result/",prefix, "_", microbe.name, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
    }
  }
}else{
  print(paste0("No such bacteria in the cancer selected"))
}
# RUN example
# Rscript script4_differential_analysis_diy.R MultipleCancersComparison UCEC CESC S s__Cercospora_beticola Fungus Yes 6 6
# Rscript script4_differential_analysis_diy.R SampleType CESC S s__Saccharomyces_eubayanus Fungus Yes 6 6