options(warn = -1)
suppressMessages(library(dplyr))
suppressMessages(library(phyloseq))
suppressMessages(library(MicrobiotaProcess))
suppressMessages(library(vegan))
suppressMessages(library(picante))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(scales))
suppressMessages(library(ggsci))

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}

args <- commandArgs(T)

#对应参数Taxa
taxa <- args[1]
tax_name <- args[2]
compare <- args[3]
plot.width <- args[4] %>% as.numeric()
plot.height <- args[5] %>% as.numeric()
Kingdom <- args[6]

datadir <- "/data/jinchuandi/";
# demo
# setwd("D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/")
# datadir <-"D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/"
# taxa <- "genus"
# tax_name <- "g__Helicobacter"
# compare <- "tumor_vs_non_tumor"
# plot.width <- 10
# plot.height <- 5
# Kingdom <- "Bacteria"



if (Kingdom == "Fungus") {
  filename <- "rawdata_fungi/"
} else {
  filename <- "rawdata/"
}


taxa_name_list <- read.csv(paste0(datadir, "rawdata/taxaall.csv"),header = T,check.names = F)

if (tax_name %in% taxa_name_list[,which(colnames(taxa_name_list)==taxa)]){
  print("The entered name is recognizable")
  phy_data_tax <- readRDS(paste0(datadir,filename, "phy_data_",taxa,"_logcpm_all_cancer.RData"))
  phy_data_tax_sub <- prune_taxa(rownames(phy_data_tax@otu_table)==tax_name,phy_data_tax)
  otu <- data.frame(phy_data_tax_sub@otu_table) %>%
    t()%>%
    data.frame()
  otu$name <- rownames(otu)
  otu$name <- gsub(".","-",otu$name,fixed = T)
  meta <- data.frame(phy_data_tax_sub@sam_data)
  meta$name <- rownames(meta)
  meta$name <- gsub(".","-",meta$name,fixed = T)
  mic_dat <- inner_join(otu,meta)
  colnames(mic_dat)[1] <- "micro_name"
  
  filter <- aggregate(mic_dat$micro_name,by=list(mic_dat$cancer_type),sum)
  name_cancer <- filter$Group.1[!filter$x==0]
  mic_dat_sub <- mic_dat[mic_dat$cancer_type %in% name_cancer,]
  mic_dat_sub$log_mic <- log(mic_dat_sub$micro_name+1,2)
  
  if (compare == "tumor_vs_non_tumor") {
    pan_mic_boxplot <- ggplot2::ggplot(mic_dat_sub, aes(x=cancer_type, y=log_mic, fill= SampleType))+
      geom_boxplot()+labs(x="Cancer Type", y="Log2 Microbial Abundance")+
      scale_fill_igv()+
      theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=.5))+stat_compare_means(label = "p.signif")
    pan_mic_boxplot
    colnames(mic_dat_sub) <- c("Microbiota Abundance","Sample ID","Sample Type","Cancer Type","Log2 Microbial Abundance")
    # 写出结果list和plot
    ggsave(pan_mic_boxplot,filename=paste0("result/pancancer_microbial_abundance_",tax_name,"_",compare,"_plot.pdf"),
           width=plot.width,height=plot.height,units=c("in"),device = cairo_pdf)
    write.csv(mic_dat_sub,file=paste0("result/pancancer_microbial_abundance_",tax_name,"_",compare,"_list.csv"))
  } else {
    if (compare == "only_tumor"){
      mic_dat_sub2 <- mic_dat_sub[mic_dat_sub$SampleType=="tumor",]
      pan_mic_boxplot2 <- ggplot2::ggplot(mic_dat_sub2, aes(x=cancer_type, y=log_mic, fill= cancer_type))+
        geom_boxplot()+labs(x="Cancer Type", y="Log2 Microbial Abundance")+
        scale_fill_igv()+
        theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=.5))+stat_compare_means(label = "p.signif")
      pan_mic_boxplot2
      colnames(mic_dat_sub2) <- c("Microbiota Abundance","Sample ID","Sample Type","Cancer Type","Log2 Microbial Abundance")
      # 写出结果list和plot
      ggsave(pan_mic_boxplot2,filename=paste0("result/pancancer_microbial_abundance_",tax_name,"_",compare,"_plot.pdf"),
             width=plot.width,height=plot.height,units=c("in"),device = cairo_pdf)
      write.csv(mic_dat_sub2,file=paste0("result/pancancer_microbial_abundance_",tax_name,"_",compare,"_list.csv"))
    }else{
      print("Some required options are not selected.")
    }
  }
}else {
  print("The entered name is not recognizable maybe due to spelling mistakes")
}

# RUN Example
# Rscript script_pancancer_microbial_abundance.R genus g__Helicobacter tumor_vs_non_tumor 10 5 Bacteria
