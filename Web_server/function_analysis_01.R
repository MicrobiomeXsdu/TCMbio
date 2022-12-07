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

# 输入参数
cancer_type = args[1]
type = args[2]
p_value <- args[3] %>% as.numeric()
group <- args[4]
sub1 <- args[5]
sub2 <- args[6]
plot.width <- args[7] %>% as.numeric()
plot.height <- args[8] %>% as.numeric()
Kingdom <- args[9]
Download_single_plot <- args[10]
datadir <- '/data/jinchuandi/'
# #demo
# setwd("D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/")
# datadir <- "D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/"
# datadir <- '/data/jinchuandi/'
# cancer_type = "LUSC"
# type = "gene"
# type = "pathway"
# p_value <- 0.30
# group <- "gender"
# group <- "Immune.Landscape"
# group <- "Conserved.TME.subtype"
# # #gender的时候，sub1和2是NULL
# sub1 <- NULL
# sub2 <- NULL
# # #Immune.Landscape的时候，sub1和2是C1:C6,六选2
# sub1 <-"C1"
# sub2 <- "C2"
# # #Conserved.TME.subtype，sub1和2如下。
# sub1 <-"IE"
# sub2 <- "F"
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
    print(paste0("No functional data for this type ",cancer_type," of ",Kingdom))
    q()
  }
}

obj <- try(dataset <- read.csv(paste0(datadir,filename,cancer_type,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1),silent=TRUE)
if (is(obj, "try-error")){
  q()
} else{
}

#name <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
dataset <- read.csv(paste0(datadir,filename,cancer_type,"_eggnog_function_data.csv"),check.names = F,header = T,row.names = 1)

if (type=="gene"){
  data_sub <- dataset %>%
    select(c(starts_with("K",ignore.case=F),group,SampleType))
} else {
  data_sub <- dataset %>%
    select(c(starts_with("ko",ignore.case=F),starts_with("map",ignore.case=F),group,SampleType))
}
data_sub$SampleType <- as.factor(data_sub$SampleType)


if (group == "gender"){
  data_sub2 <- data_sub
  data_sub2$gender <- as.factor(data_sub2$gender)
  data_sub2 <- subset(data_sub2,data_sub2$SampleType=="1")
  data_sub2 <- data_sub2[,!colnames(data_sub2) %in% "SampleType"]
  data_sub2 <- data_sub2[which(!is.na(data_sub2$gender)),]
  if (nrow(data_sub2)==0 | length(unique(data_sub2$gender))==1) {
    print("Under this condition, the sample size is insufficient.")
  } else {
    print("Under this condition, the sample size is sufficient.")
    result_table <- matrix(ncol = 3,nrow = length(colnames(data_sub2))-1)
    result_table <- data.frame(result_table)
    colnames(result_table) <- c("name","p_value","adj_p_value")
    result_table$name <- colnames(data_sub2)[1:length(colnames(data_sub2))-1]
    for (i in 1:(length(colnames(data_sub2))-1)) {
      res_result <- wilcox.test(data_sub2[,i]~data_sub2$gender)
      result_table[i,2] <- res_result[["p.value"]]
      if (Download_single_plot == "Yes") {
        if(res_result[["p.value"]]<p_value){
          pic <- ggplot(data_sub2, aes(x=cancer_type,
                                       y=colnames(data_sub2)[i],
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
          ggsave(pic,filename=paste0("result/function_result_output/function_result_",cancer_type,"_",colnames(data_sub2)[i],"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_single_barplot_result.pdf"),width=plot.width,height=plot.height,units=c("in"))
        } else {
          next
        }
      } else {
        print("Not download single plot")
      }
    }
    result_table$adj_p_value <- p.adjust(result_table$p_value)
    
    #写出list
    write.csv(result_table,file = paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_table_result.csv"))
    bar_data <-   aggregate(data_sub2[,1:length(colnames(data_sub2))-1], by = list(data_sub2$gender),FUN=mean)
    colnames(bar_data)[1] <- "group"
    rownames(bar_data) <- bar_data$group
    bar_data_sig <- bar_data[,colnames(bar_data) %in% result_table$name[which(result_table$adj_p_value<p_value)]]
    bar_data_sig <- data.frame(bar_data_sig)
    colnames(bar_data_sig) <- result_table$name[which(result_table$adj_p_value<p_value)]
    bar_data_sig$group <- rownames(bar_data)
    
    if (!ncol(bar_data_sig) == 1) {
      bar_data_sig2 <- reshape2::melt(bar_data_sig)
      bar_plot <-ggplot(bar_data_sig2, aes(x=variable,y=value)) +
        geom_bar(aes(fill = group), stat = "identity",color="black",size=0.8,
                 position = position_dodge(0.5), width = 0.5)+
        scale_fill_manual(values = c("#C77CFF","#00BFC4")) +
        labs(title = NULL, x = NULL, y = 'value', fill = NULL) +
        theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
        theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) 
      bar_plot
      
      #写出plot
      ggsave(bar_plot, filename=paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_barplot_result.pdf"), width = plot.width, height = plot.height, units = c("in"))
    } else{
      print(paste0("There is no significant ", type," under this P value in ",cancer_type," cancer.")) 
    }
  }
} else {
  if (group == "Immune.Landscape"){
    data_sub2 <- data_sub
    data_sub2$Immune.Landscape <- as.factor(data_sub2$Immune.Landscape)
    data_sub2 <- subset(data_sub2,data_sub2$SampleType=="1")
    data_sub2 <- data_sub2[,!colnames(data_sub2) %in% "SampleType"]
    data_sub2 <- subset(data_sub2,data_sub2$Immune.Landscape==sub1 | data_sub2$Immune.Landscape==sub2)
    data_sub2 <- data_sub2[which(!is.na(data_sub2$Immune.Landscape)),]
    if (nrow(data_sub2)==0 | length(unique(data_sub2$Immune.Landscape))==1) {
      print("Under this condition, the sample size is insufficient.")
    } else {
      print("Under this condition, the sample size is sufficient.")
      result_table <- matrix(ncol = 3,nrow = length(colnames(data_sub2))-1)
      result_table <- data.frame(result_table)
      colnames(result_table) <- c("name","p_value","adj_p_value")
      result_table$name <- colnames(data_sub2)[1:length(colnames(data_sub2))-1]
      for (i in 1:(length(colnames(data_sub2))-1)) {
        res_result <- wilcox.test(data_sub2[,i]~data_sub2$Immune.Landscape)
        result_table[i,2] <- res_result[["p.value"]]
        if (Download_single_plot == "Yes") {
          if(res_result[["p.value"]]<p_value){
            pic <- ggplot(data_sub2, aes(x=cancer_type,
                                         y=colnames(data_sub2)[i],
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
            ggsave(pic,filename=paste0("result/function_result_output/function_result_",cancer_type,"_",colnames(data_sub2)[i],"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_single_barplot_result.pdf"),width=plot.width,height=plot.height,units=c("in"))
          } else {
            next
          }
        } else {
          print("Not download single plot")
        }
      }
      result_table$adj_p_value <- p.adjust(result_table$p_value)
      
      #写出list
      write.csv(result_table,file = paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_table_result.csv"))
      bar_data <-   aggregate(data_sub2[,1:length(colnames(data_sub2))-1], by = list(data_sub2$Immune.Landscape),FUN=mean)
      colnames(bar_data)[1] <- "group"
      rownames(bar_data) <- bar_data$group
      bar_data_sig <- bar_data[,colnames(bar_data) %in% result_table$name[which(result_table$adj_p_value<p_value)]]
      bar_data_sig <- data.frame(bar_data_sig)
      colnames(bar_data_sig) <- result_table$name[which(result_table$adj_p_value<p_value)]
      bar_data_sig$group <- rownames(bar_data)
      
      if (!ncol(bar_data_sig) == 1) {
        bar_data_sig2 <- reshape2::melt(bar_data_sig)
        bar_plot <-ggplot(bar_data_sig2, aes(x=variable,y=value)) +
          geom_bar(aes(fill = group), stat = "identity",color="black",size=0.8,
                   position = position_dodge(0.5), width = 0.5)+
          scale_fill_manual(values = c("#C77CFF","#00BFC4")) +
          labs(title = NULL, x = NULL, y = 'value', fill = NULL) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
          theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) 
        bar_plot
        
        #写出plot
        ggsave(bar_plot, filename=paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_barplot_result.pdf"), width = plot.width, height = plot.height, units = c("in"))
      } else{
        print(paste0("There is no significant ", type," under this P value in ",cancer_type," cancer.")) 
      }
    }
  } else {
    if(group == "Conserved.TME.subtype"){
      data_sub2 <- data_sub
      data_sub2$Conserved.TME.subtype <- as.factor(data_sub2$Conserved.TME.subtype)
      data_sub2 <- subset(data_sub2,data_sub2$SampleType=="1")
      data_sub2 <- data_sub2[,!colnames(data_sub2) %in% "SampleType"]
      data_sub2 <- subset(data_sub2,data_sub2$Conserved.TME.subtype==sub1 | data_sub2$Conserved.TME.subtype==sub2)
      data_sub2 <- data_sub2[which(!is.na(data_sub2$Conserved.TME.subtype)),]
      if (nrow(data_sub2)==0 | length(unique(data_sub2$Conserved.TME.subtype))==1) {
        print("Under this condition, the sample size is insufficient.")
      } else {
        print("Under this condition, the sample size is sufficient.")
        result_table <- matrix(ncol = 3,nrow = length(colnames(data_sub2))-1)
        result_table <- data.frame(result_table)
        colnames(result_table) <- c("name","p_value","adj_p_value")
        result_table$name <- colnames(data_sub2)[1:length(colnames(data_sub2))-1]
        if (Download_single_plot == "Yes") {
          if(res_result[["p.value"]]<p_value){
            pic <- ggplot(data_sub2, aes(x=cancer_type,
                                         y=colnames(data_sub2)[i],
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
            ggsave(pic,filename=paste0("result/function_result_output/function_result_",cancer_type,"_",colnames(data_sub2)[i],"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_single_barplot_result.pdf"),width=plot.width,height=plot.height,units=c("in"))
          } else {
            next
          }
        } else {
          print("Not download single plot")
        }
      }
      result_table$adj_p_value <- p.adjust(result_table$p_value)
      
      #写出list
      write.csv(result_table,file = paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_table_result.csv"))
      bar_data <-   aggregate(data_sub2[,1:length(colnames(data_sub2))-1], by = list(data_sub2$Conserved.TME.subtype),FUN=mean)
      colnames(bar_data)[1] <- "group"
      rownames(bar_data) <- bar_data$group
      bar_data_sig <- bar_data[,colnames(bar_data) %in% result_table$name[which(result_table$adj_p_value<p_value)]]
      bar_data_sig <- data.frame(bar_data_sig)
      colnames(bar_data_sig) <- result_table$name[which(result_table$adj_p_value<p_value)]
      bar_data_sig$group <- rownames(bar_data)
      
      if (!ncol(bar_data_sig) == 1) {
        bar_data_sig2 <- reshape2::melt(bar_data_sig)
        bar_plot <-ggplot(bar_data_sig2, aes(x=variable,y=value)) +
          geom_bar(aes(fill = group), stat = "identity",color="black",size=0.8,
                   position = position_dodge(0.5), width = 0.5)+
          scale_fill_manual(values = c("#C77CFF","#00BFC4")) +
          labs(title = NULL, x = NULL, y = 'value', fill = NULL) +
          theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
          theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) 
        bar_plot
        
        #写出plot
        ggsave(bar_plot, filename=paste0("result/function_result_",cancer_type,"_",type,"_",group,"_",sub1,"_",sub2,"_",p_value,"_barplot_result.pdf"), width = plot.width, height = plot.height, units = c("in"))
      } else{
        print(paste0("There is no significant ", type," under this P value in ",cancer_type," cancer.")) 
      } 
      
    } else {
      print("error")
    }
  }
}



# Rscript function_analysis_01.R ACC function 0.30 gender NULL NULL 5 4 Fungus No