library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(phyloseq)
pacman::p_load(tidyverse,ape,ggsignif,patchwork,multcomp)
require("RColorBrewer")

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}
# 设置原始文件目录
datadir <- '/data/jinchuandi/'
# datadir <- './'
# 准备数据
# 参数读取
args <- commandArgs(T)
group = args[1]
if(group!="MultipleCancersComparison"){
  CancerType = prefix = args[2]
  Differential.Method = args[3]
  level = tolower(substr(args[4],1,1))
  Kingdom = args[5]
  if(Differential.Method=="LEfSe"){
    strict.mode = ifelse(args[6]=="Yes", TRUE, FALSE)
    Adjust.Method= tolower(args[7])
    P.value.Cutoff = args[8] %>% as.numeric()
    LDA.Cutoff = args[9] %>% as.numeric()
    plot.height = args[10] %>% as.numeric()
    plot.width = args[11] %>% as.numeric()
  }
  if(Differential.Method=="Wilcoxon Test"){
    P.value.Cutoff = args[6] %>% as.numeric()
    logFCcutoff = args[7] %>% as.numeric()
    Adjust.Method= tolower(args[8])
    plot.height = args[9] %>% as.numeric()
    plot.width = args[10] %>% as.numeric()
  }
  if(Differential.Method=="Random Forest"){
    top.features = args[6] %>% as.numeric()
    plot.height = args[7] %>% as.numeric()
    plot.width = args[8] %>% as.numeric()
  }
  if(Differential.Method=="DESeq2"){
    P.value.Cutoff = args[6] %>% as.numeric()
    Adjust.Method= tolower(args[7])
    logFCcutoff = args[8] %>% as.numeric()
    plot.height = args[9] %>% as.numeric()
    plot.width = args[10] %>% as.numeric()
  }
}else{
  CancerTypeA = prefix = args[2]
  CancerTypeB = prefix2 = args[3]
  Differential.Method = args[4]
  level = tolower(substr(args[5],1,1))
  Kingdom = args[6]
  if(Differential.Method=="LEfSe"){
    strict.mode = ifelse(args[7]=="Yes", TRUE, FALSE)
    Adjust.Method= tolower(args[8])
    P.value.Cutoff = args[9] %>% as.numeric()
    LDA.Cutoff = args[10] %>% as.numeric()
    plot.height = args[11] %>% as.numeric()
    plot.width = args[12] %>% as.numeric()
  }
  if(Differential.Method=="Wilcoxon Test"){
    P.value.Cutoff = args[7] %>% as.numeric()
    logFCcutoff = args[8] %>% as.numeric()
    Adjust.Method= tolower(args[9])
    plot.height = args[10] %>% as.numeric()
    plot.width = args[11] %>% as.numeric()
  }
  if(Differential.Method=="Random Forest"){
    top.features = args[7] %>% as.numeric()
    plot.height = args[8] %>% as.numeric()
    plot.width = args[9] %>% as.numeric()
  }
  if(Differential.Method=="DESeq2"){
    Adjust.Method= tolower(args[7])
    P.value.Cutoff = args[8] %>% as.numeric()
    logFCcutoff = args[9] %>% as.numeric()
    plot.height = args[10] %>% as.numeric()
    plot.width = args[11] %>% as.numeric()
  }
}

# # # 参数示例
# group = "MultipleCancersComparison"
# if(group!="MultipleCancersComparison"){
#   CancerType = prefix = "BRCA"
#   Differential.Method = "Wilcoxon Test" #"LEfSe" "Wilcoxon Test" "Random Forest" "DESeq2"
#   level = tolower(substr("Species",1,1))
#   Kingdom = "Fungus"
#   if(Differential.Method=="LEfSe"){
#     strict.mode = ifelse("Yes"=="Yes", TRUE, FALSE)
#     Adjust.Method= tolower("FDR")
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     LDA.Cutoff = "3" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="Wilcoxon Test"){
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     logFCcutoff = "1" %>% as.numeric()
#     Adjust.Method= tolower("FDR")
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="Random Forest"){
#     top.features = "20" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="DESeq2"){
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     Adjust.Method= tolower("FDR")
#     logFCcutoff = "1" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
# }else{
#   CancerTypeA = prefix = "CESC"
#   CancerTypeB = prefix2 = "UCEC"
#   Differential.Method = "DESeq2"
#   level = tolower(substr("Species",1,1))
#   Kingdom = "Fungus"
#   if(Differential.Method=="LEfSe"){
#     strict.mode = ifelse("Yes"=="Yes", TRUE, FALSE)
#     Adjust.Method= tolower("FDR")
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     LDA.Cutoff = "3" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="Wilcoxon Test"){
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     logFCcutoff = "1" %>% as.numeric()
#     Adjust.Method= tolower("FDR")
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="Random Forest"){
#     top.features = "FDR" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
#   if(Differential.Method=="DESeq2"){
#     P.value.Cutoff = "0.05" %>% as.numeric()
#     Adjust.Method= tolower("FDR")
#     logFCcutoff = "1" %>% as.numeric()
#     plot.height = "6" %>% as.numeric()
#     plot.width = "6" %>% as.numeric()
#   }
# }

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

# 读取数据并绘图

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
  phy_dat <- phyloseq(otu_tab, metadata)
}else{
  phy_dat <- readRDS(phydata)
}
# try(names(phy_dat@sam_data)[which(names(phy_dat@sam_data)=="SampleType")] <- "SampleType", silent = T)
try(names(phy_dat@sam_data)[which(names(phy_dat@sam_data)=="CancerType")] <- "MultipleCancersComparison", silent = T)

phy_dat <- prune_samples(!is.na(phy_dat@sam_data[[group]])&sample_sums(phy_dat)>0, phy_dat)
group_name <- c("group1", "group2", "group3", "group4")
phy_dat@sam_data[["group"]] <- NA
for (i in 1:length(unique(phy_dat@sam_data[[group]]))) {
  phy_dat@sam_data[["group"]] <- ifelse(phy_dat@sam_data[[group]] == unique(phy_dat@sam_data[[group]])[i], group_name[i],  phy_dat@sam_data[["group"]])
}


###### LEfSe #####

if(Differential.Method=="LEfSe"){
  diffres <- NULL
  try(diffres <- diff_analysis(phy_dat, classgroup = "group", padjust = Adjust.Method, 
                               firstalpha = P.value.Cutoff, strictmod = strict.mode, 
                               ldascore = LDA.Cutoff),silent = TRUE)
  if(!is.null(diffres)){
    result <- NULL
    for (i in 1:length(diffres@secondvars)) {
      result <- c(result, as.character(diffres@secondvars[[i]][,2]))
    }
    if(!all(result=="FALSE")){
      for (i in names(diffres@secondvars)) {
        colnames(diffres@secondvars[[i]])[2] <- i
      }
      lefres <- diffres@secondvars[[names(diffres@secondvars)[1]]][c(1,2)]
      try(
        for (i in names(diffres@secondvars)[-1]) {
          lefres <- merge(lefres, diffres@secondvars[[i]][c(1,2)], by = "f", all = T)
        }
      )
      lefres <- left_join(lefres, diffres@mlres)
      for (i in 2:(ncol(lefres)-3)) {
        lefres[,i] <- as.logical(lefres[,i])
      }
      
      lefres[is.na(lefres)] <- FALSE
      lefres$group <- rep(NA, nrow(lefres)) %>% as.character()
      for (i in names(diffres@secondvars)) {
        j <- substr(i, nchar(i), nchar(i)) %>% as.numeric()
        lefres$group <- ifelse(lefres[,i]==TRUE, unique(phy_dat@sam_data[[group]])[j], lefres$group)
      }
      lefres$group <- paste0(toupper(substr(lefres$group, 1,1)), substr(lefres$group, 2,nchar(lefres$group)))
      lefres <- lefres[order(lefres$group,lefres$LDAmean),] %>% na.omit()
      lefres$f <- factor(lefres$f, levels = lefres$f)
      lefres <- lefres[c(1,(ncol(lefres)-3):ncol(lefres))]
      names(lefres) <- c("Bateria", "LDAupper", "LDAmean", "LDAlower", "Group")
      lefres <- lefres[lefres$LDAmean>LDA.Cutoff,]
      p1 <- ggplot(lefres, aes(Bateria, LDAmean, fill = Group)) + 
        geom_bar(stat = 'identity',colour="black", position=position_dodge(0.3), width=1) +
        scale_fill_manual(values = c('#167C3B', '#E61F19', 	"#3b167c", "#7c3b16")[1:length(unique(phy_dat@sam_data[[group]]))]) +
        xlab("") +
        ylab("LDA")+
        theme_bw() + 
        theme(panel.grid =element_blank()) + 
        theme(panel.border = element_blank()) + 
        theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) + 
        coord_flip()
      if (group=="MultipleCancersComparison") {
        ggsave(p1,filename=paste0("result/",prefix, "_", prefix2, "_LEfSe_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
        write.csv(lefres, paste0("result/", prefix, "_", prefix2, "_LEfSe_", level, "_", group, ".csv"), row.names = F)
      }else{
        ggsave(p1,filename=paste0("result/", prefix, "_LEfSe_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
        write.csv(lefres, paste0("result/", prefix, "_LEfSe_", level, "_", group, ".csv"), row.names = F)
      }
      
    }else{
      print("There is no significant differential microbes between groups.")
    }
  }else{
    print("There is no significant differential microbes between groups.")
  }
}


###### wilcoxon ######

if(Differential.Method=="Wilcoxon Test"){
  if(length(unique(phy_dat@sam_data[[group]]))==2){
    if(phy_dat@otu_table@taxa_are_rows){
      bac_table <- t(phy_dat@otu_table@.Data)%>%
        as.data.frame()
    }else{
      bac_table <- phy_dat@otu_table@.Data%>%
        as.data.frame()
    }
    bac_table$group <- phy_dat@sam_data[[group]]
    bac_table$group <- bac_table$group <- paste0(toupper(substr(bac_table$group, 1,1)), 
                                                 substr(bac_table$group, 2,nchar(bac_table$group)))
    
    p.value <- data.frame(p = rep(NA, length(colnames(bac_table))-1),
                          p.adj = rep(NA, length(colnames(bac_table))-1),
                          logFC = rep(NA, length(colnames(bac_table))-1)
    )
    row.names(p.value) <- colnames(bac_table)[-ncol(bac_table)]
    for (bac in row.names(p.value)) {
      test <- wilcox.test(bac_table[,bac]~bac_table$group)
      p.value[bac, 1] <- test$p.value
      group1 <- bac_table[bac_table$group==unique(bac_table$group)[1], bac]
      group2 <- bac_table[bac_table$group==unique(bac_table$group)[2], bac]
      p.value[bac, 3] <- log2(mean(group1+1)/mean(group2+1))
    }
    if(Adjust.Method == "none"){
      p.value$p.adj <- p.value$p
    }
    if(Adjust.Method == "fdr"){
      p.value$p.adj <- p.adjust(p.value$p, 'fdr')
    }
    if(Adjust.Method == "bonferroni"){
      p.value$p.adj <- p.adjust(p.value$p, 'bonferroni')
    }
    p.value$neg_logQ <- c(-log10(p.value$p.adj))
    p.value$neg_logQ <- c(-log10(p.value$p.adj))
    
    p.value[which(p.value$p.adj <= P.value.Cutoff & p.value$logFC < 0-logFCcutoff),'sig'] <- 'Down'
    p.value[which(p.value$p.adj <= P.value.Cutoff & p.value$logFC > logFCcutoff),'sig'] <- 'Up'
    p.value[which(p.value$p.adj > P.value.Cutoff | p.value$logFC < logFCcutoff &
                    p.value$logFC > 0-logFCcutoff),'sig'] <- 'None'
    
    p <- ggplot(p.value, aes(x = logFC, y = neg_logQ, color = sig)) +
      geom_point(alpha = 0.6, size = 3) +
      scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
      theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
      geom_vline(xintercept = logFCcutoff, color = 'gray', size = 0.3) +
      geom_vline(xintercept = 0-logFCcutoff, color = 'gray', size = 0.3) +
      geom_hline(yintercept = -log(P.value.Cutoff, 10), color = 'gray', size = 0.3) +
      labs(x = 'Log2 Fold Change', y = '-Log(P Value)', color = '', title= paste0(group, ": ", unique(bac_table$group)[1] ," v.s.", unique(bac_table$group)[2]))
    
    p.value <- p.value[,-which(colnames(p.value)=="neg_logQ")]%>%
      tibble::rownames_to_column("Bacteria")
    names(p.value) <- c("Bacteria", "P-value", "Adjusted P-value", "log(FC)", "Significance")
    
    if (group=="CancerType") {
      ggsave(p,filename=paste0("result/",prefix, "_", prefix2, "_wilcoxon_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
      write.csv(p.value, paste0("result/", prefix, "_", prefix2, "_wilcoxon_", level, "_", group, ".csv"),row.names = F)
    }else{
      ggsave(p,filename=paste0("result/", prefix, "_wilcoxon_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
      write.csv(p.value, paste0("result/", prefix, "_wilcoxon_", level, "_", group, ".csv"),row.names = F)
    }
  }else{
    stop("The number of groups is not 2.")
  }
}




###### random forest ######

if(Differential.Method=="Random Forest"){
  library(randomForest)
  if(length(unique(phy_dat@sam_data[[group]]))==2){
    if(phy_dat@otu_table@taxa_are_rows){
      bac_table <- t(phy_dat@otu_table@.Data)%>%
        as.data.frame()
    }else{
      bac_table <- phy_dat@otu_table@.Data%>%
        as.data.frame()
    }
    bac_table$group <- factor(phy_dat@sam_data[["group"]])
    
    colnames(bac_table) <- gsub("[", ".", colnames(bac_table), fixed = T)
    colnames(bac_table) <- gsub("]", ".", colnames(bac_table), fixed = T)
    colnames(bac_table) <- gsub("-", ".", colnames(bac_table), fixed = T)
    colnames(bac_table) <- gsub("/", ".", colnames(bac_table), fixed = T)
    colnames(bac_table) <- gsub("(", ".", colnames(bac_table), fixed = T)
    colnames(bac_table) <- gsub(")", ".", colnames(bac_table), fixed = T)
    set.seed(1234)
    Groups.rf = randomForest(group ~ ., data=bac_table, ntree=500, importance=TRUE, proximity=TRUE)
    datares <- round(importance(Groups.rf), 2) %>%
      as.data.frame()
    Gini_top <- rownames(datares[order(datares$MeanDecreaseGini,decreasing = T),])[1:top.features]
    Acc_top <- rownames(datares[order(datares$MeanDecreaseAccuracy,decreasing = T),])[1:top.features]
    top <- unique(c(Gini_top, Acc_top))
    table <- data.frame(Bacteria = top,
                        Gini = datares[top,"MeanDecreaseGini"],
                        Accuracy = datares[top,"MeanDecreaseAccuracy"])
    if (group=="MultipleCancersComparison") {
      pdf(paste0("result/",prefix, "_", prefix2, "_random_forest_", level, "_", group, ".pdf"),height = plot.height, width = plot.width)
      varImpPlot(Groups.rf,sort = TRUE, n.var = top.features)
      dev.off()
      write.csv(table, paste0("result/", prefix, "_", prefix2, "_random_forest_", level, "_", group, ".csv"),row.names = F)
    }else{
      pdf(paste0("result/",prefix, "_random_forest_", level, "_", group, ".pdf"),height = plot.height, width = plot.width)
      varImpPlot(Groups.rf,sort = TRUE, n.var = top.features)
      dev.off()
      write.csv(table, paste0("result/", prefix, "_random_forest_", level, "_", group, ".csv"),row.names = F)
    }
    
    
  }else{
    print("The number of groups is not 2.")
  }
}

###### DESeq2 ######
if(Differential.Method=="DESeq2"){
  library(DESeq2)
  if(length(unique(phy_dat@sam_data[[group]]))==2){
    if(phy_dat@otu_table@taxa_are_rows){
      bac_table <- phy_dat@otu_table@.Data+1
    }else{
      bac_table <- phy_dat@otu_table@.Data+1%>%
        t()
    }
    sample <- data.frame(phy_dat@sam_data)
    sample$group <- paste0(toupper(substr(sample[,group], 1,1)), 
                           substr(sample[,group], 2,nchar(sample[,group])))
    sample$group <- factor(sample$group,levels = unique(sample$group))
    
    sample <- subset(sample, select = group)
    dds <- DESeqDataSetFromMatrix(countData=bac_table, colData=sample, design= ~group)
    dds2 <- DESeq(dds)
    res <- results(dds2, contrast=c("group",as.character(unique(sample[,"group"]))),
                   pAdjustMethod = Adjust.Method)
    # resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    res <- as.data.frame(res)
    res$padj[is.na(res$padj)] <- 1
    res$g=ifelse(res$pvalue>P.value.Cutoff,'stable',
                 ifelse(res$log2FoldChange > logFCcutoff,'up',
                        ifelse(res$log2FoldChange < 0-logFCcutoff,'down','stable')))
    res$neg_log_Q <- 0-log10(res$padj)
    table(res$g)
    p <- ggplot(res, aes(x = log2FoldChange, y = neg_log_Q, color = g)) +
      geom_point(alpha = 0.6, size = 3) +
      scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c("up", "down", "stable")) +
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
      theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
      geom_vline(xintercept = logFCcutoff, color = 'gray', size = 0.3) +
      geom_vline(xintercept = 0-logFCcutoff, color = 'gray', size = 0.3) +
      geom_hline(yintercept = -log(P.value.Cutoff, 10), color = 'gray', size = 0.3) +
      labs(x = 'Log2 Fold Change', y = '-Log10(Q Value)', color = '', title= paste0(group, ": ", as.character(unique(sample[,"group"]))[1] ,
                                                                                    " v.s.", as.character(unique(sample[,"group"]))[2]))
    res <- res[-which(colnames(res)%in%c("baseMean", "lfcSE", "stat",
                                         "neg_log_Q"))]%>%
      tibble::rownames_to_column("Bacteria")
    names(res) <- c("Bacteria", "log2(FoldChange)", "P-value", 
                    "Adjusted P-value", "Group")
    if (group=="MultipleCancersComparison") {
      ggsave(p,filename=paste0("result/",prefix, "_", prefix2, "_DESeq2_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
      write.csv(res, paste0("result/", prefix, "_", prefix2, "_DESeq2_", level, "_", group, ".csv"),row.names = F)
    }else{
      ggsave(p,filename=paste0("result/", prefix, "_DESeq2_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"))
      write.csv(res, paste0("result/", prefix, "_DESeq2_", level, "_", group, ".csv"),row.names = F)
    }
  }else{
    print("Groups are more than 2.")
  }
}
# RUN example
# Rscript script3_differential_analysis.R SampleType CESC LEfSe Species Fungus Yes FDR 0.05 3 6 6
# Rscript script3_differential_analysis.R SampleType CESC DESeq2 Species Fungus 0.05 FDR 1 6 6
# Rscript script3_differential_analysis.R MultipleCancersComparison CESC UCEC DESeq2 Species Fungus FDR 0.05 1 6 6
# Rscript script3_differential_analysis.R MultipleCancersComparison CESC UCEC LEfSe Species Fungus Yes FDR 0.05 3 6 6