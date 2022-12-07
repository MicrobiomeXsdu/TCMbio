options(warn = -1)
suppressMessages(library(rjson))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}

if (!dir.exists('result/function_result_output/')) {
  dir.create('result/function_result_output/')
}

args <- commandArgs(T)

cancer_type_a = args[1]
cancer_type_b = args[2]
type = args[3]
p_value <- args[4] %>% as.numeric()
plot.width <- args[5] %>% as.numeric()
plot.height <- args[6] %>% as.numeric()
Kingdom <- args[7]
Download_single_plot <- args[8]
datadir <- '/data/jinchuandi/'
# #demo
# setwd("D:/anlysis-data-r/pancacer_rna/rdata_decontaminate")
# cancer_type_a = "LUSC"
# cancer_type_b = "PAAD"
# type = "gene"
# p_value <- 0.01
# plot.width <- 5
# plot.height <- 4
# Kingdom <- "Fungus"
# Download_single_plot <- "No"

if (Kingdom == "Fungus") {
  filename <- "rawdata_fungi/"
} else {
  if (Kingdom == "Bacteria"){
    filename <- "rawdata/"
  } else {
    print(paste0("No functional data for this type ",cancer_type_a,"_",cancer_type_b,"_"," of ",Kingdom))
    q()
  }
}

obj <- try(dataset_a <- read.csv(paste0(datadir,filename,cancer_type_a,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1),silent=TRUE)
if (is(obj, "try-error")){
  q()
} else{
}
dataset_a <- read.csv(paste0(datadir,filename,cancer_type_a,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1)

if (type=="gene"){
  data_suba <- dataset_a %>%
    select(c(starts_with("K",ignore.case=F),SampleType))
} else {
  data_suba <- dataset_a %>%
    select(c(starts_with("ko",ignore.case=F),starts_with("map",ignore.case=F),SampleType))
}

obj <- try(dataset_b <- read.csv(paste0(datadir,filename,cancer_type_b,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1),silent=TRUE)
if (is(obj, "try-error")){
  q()
} else{
}

dataset_b <- read.csv(paste0(datadir,filename,cancer_type_b,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1)

if (type=="gene"){
  data_subb <- dataset_b %>%
    select(c(starts_with("K",ignore.case=F),SampleType))
} else {
  data_subb <- dataset_b %>%
    select(c(starts_with("ko",ignore.case=F),starts_with("map",ignore.case=F),SampleType))
}

name_same <- intersect(colnames(data_suba),colnames(data_subb))

data_suba2 <- data_suba[,colnames(data_suba) %in% name_same]
data_suba2$cancer_type <- cancer_type_a
data_subb2 <- data_subb[,colnames(data_subb) %in% name_same]
data_subb2$cancer_type <- cancer_type_b

data_all <- rbind(data_suba2,data_subb2)
data_all$SampleType <- as.factor(data_all$SampleType)
if (cancer_type_a == "LAML" | cancer_type_b == "LAML"){
  data_all2 <- data_all[data_all$SampleType=="1"|data_all$SampleType=="3",]
} else {
  data_all2 <- data_all[data_all$SampleType=="1",]
}

length  <- length(colnames(data_all2))-2
result_table <- matrix(ncol = 3,nrow = length(colnames(data_all2))-2)
result_table <- data.frame(result_table)
colnames(result_table) <- c("name","p_value","adj_p_value")

result_table$name <- colnames(data_all2)[1:length]
for (i in 1:length) {
  res_result <- wilcox.test(data_all2[,i]~data_all2$cancer_type)
  result_table[i,2] <- res_result[["p.value"]]
  if (Download_single_plot == "Yes") {
    if(res_result[["p.value"]]<p_value){
      pic <- ggplot(data_all2, aes(x=cancer_type,
                                   y=colnames(data_all2)[i],
                                   fill=cancer_type)) +
        geom_boxplot()+
        stat_boxplot(geom = "errorbar",
                     lwd=0.5,
                     width=0.2)+
        geom_jitter(color="black", size=0.8, alpha=0.9)+
        theme_minimal()+
        scale_fill_nejm()+
        theme(legend.position ="none")+
        stat_compare_means(label = "p.signif")+
        labs(y=NULL)
      pic
      #写出plot
      ggsave(pic,filename=paste0("result/function_result_output/function_result_",cancer_type_a,"_",cancer_type_b,"_",colnames(data_all2)[i],"_",type,"_","_",p_value,"_single_barplot_result.pdf"),width=plot.width,height=plot.height,units=c("in"))
    } else {
      next
    }
  } else {
    print("Not download single plot")
  }
}
result_table$adj_p_value <- p.adjust(result_table$p_value)

#写出list
write.csv(result_table,file = paste0("result/function_result_",cancer_type_a,"_",cancer_type_b,"_",type,"_list.csv"))
bar_data <-   aggregate(data_all2[,1:length], by = list(data_all2$cancer_type),FUN=mean)
colnames(bar_data)[1] <- "group"
rownames(bar_data) <- bar_data$group
bar_data_sig <- bar_data[,colnames(bar_data) %in% result_table$name[which(result_table$adj_p_value<p_value)]]
bar_data_sig$group <- rownames(bar_data_sig)

if (!ncol(bar_data_sig) == 1) {
  bar_data_sig2 <- melt(bar_data_sig)
  bar_plot <-ggplot(bar_data_sig2, aes(x=variable,y=value)) +
    geom_bar(aes(fill = group), stat = "identity",color="black",size=0.8,
             position = position_dodge(0.5), width = 0.5)+
    scale_fill_manual(values = c("#C77CFF","#00BFC4")) +
    labs(title = NULL, x = NULL, y = 'value', fill = NULL) +
    theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
    theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) 
  bar_plot
  
  #写出plot
  ggsave(bar_plot, filename=paste0("result/function_result_",cancer_type_a,"_",cancer_type_b,"_",type,"_",p_value,"_barplot_result.pdf"), width = plot.width, height = plot.height, units = c("in"))
} else{
  print(paste0("There is no significant ", type," under this P value.")) 
}

# Rscript function_analysis_02.R LUSC PAAD gene 0.01 5 4 Fungus No