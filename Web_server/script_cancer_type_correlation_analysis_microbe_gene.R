options(warn = -1)
suppressMessages(library(phyloseq))
suppressMessages(library(MicrobiotaProcess))
suppressMessages(library(vegan))
suppressMessages(library(data.table))
suppressMessages(library(limma))
suppressMessages(library(reshape2))
suppressMessages(library(picante))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(ggExtra))
suppressMessages(library(ggpubr))

# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}
if (!dir.exists('result/temp')) {
  dir.create('result/temp')
}

args <- commandArgs(T)

#对应参数
cancer_type = args[1]
taxaa <- args[2]
gene <- args[3]
correlation_coef <- args[4]
plot.height <- args[5] %>% as.numeric()
plot.width <- args[6] %>% as.numeric()
Kingdom <- args[7]

datadir <- '/data/jinchuandi/'

# demo
# setwd("D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/")
# cancer_type = "LUSC"
# taxaa <- "g__Aspergillus"
# gene <- "DPM1"
# correlation_coef <- "pearson"
# plot.height=6
# plot.width=6
# Kingdom <- "Fungus"

if (Kingdom == "Fungus") {
  filename <- "rawdata_fungi/"
} else {
  filename <- "rawdata/"
}

phydata_g = paste0(datadir, filename, cancer_type, "_g", "_all_logCPM.RData")
phydata_s = paste0(datadir, filename, cancer_type, "_s", "_all_logCPM.RData")
phy_dat_g <- readRDS(phydata_g)
phy_dat_s <- readRDS(phydata_s)
if(phy_dat_g@otu_table@taxa_are_rows==TRUE){
  otu_g <- phy_dat_g@otu_table@.Data %>%
    t() %>%
    data.frame()
}else{
  otu_g <- phy_dat_g@otu_table@.Data %>%
    data.frame()
}

if(phy_dat_s@otu_table@taxa_are_rows==TRUE){
  otu_s <- phy_dat_s@otu_table@.Data %>%
    t() %>%
    data.frame()
}else{
  otu_s <- phy_dat_s@otu_table@.Data
}

metadata <- data.frame(phy_dat_g@sam_data)

taxaname <- data.frame(cbind(colnames(otu_g),colnames(otu_s)))
colnames(taxaname) <- c("genus","species")

if (taxaa %in% taxaname$genus | taxaa %in% taxaname$species ){
  print(paste0("Taxa A is ok in ",cancer_type," cancer"))
} else {
  print(paste0("Taxa A is error in ",cancer_type," cancer"))
  q()
}
deg_list <- fread(paste0(datadir, "rawdata/rna_deg_all.csv"),header = T,check.names = F)
deg_list_cancer <- deg_list[deg_list$cancer==cancer_type,gene_name]


gene_list <- read.csv(paste0(datadir, "rawdata/",cancer_type,"_gene_name_list.csv"),header = T,check.names = F)
if (is.null(gene)) {
  print("Not input Gene")
} else {
  if (gene %in% deg_list_cancer){
    print(paste0("Gene is DEGs in ",cancer_type," cancer"))
  } else {
    if (gene %in% gene_list$gene_name) {
      print(paste0("The entered gene name is recognized in"," ",cancer_type," cancer"))
    } else {
      print(paste0("The entered gene name is not recognized in"," ",cancer_type," cancer"))
      q()
    }}
}

rna_data = fread(paste0(datadir, "rawdata/TCGA_rna_counts/tpm/",cancer_type,"_tpm.csv"))

otu_g$name <- rownames(otu_g)
otu_s$name <- rownames(otu_s)
otu <- inner_join(otu_g,otu_s)
rownames(otu) <- substr(otu$name,1,15)

rna_data <- data.frame(t(rna_data))
colnames(rna_data) <- rna_data[1,]
rna_data <- rna_data[-1,]
ann <- fread(paste0(datadir, "rawdata/TCGA_rna_counts/annotation_gene.csv"),header = T)
colnames(rna_data) <- ann$gene_name[match(colnames(rna_data),ann$gene_id)]
rownames(rna_data) <- gsub(".","-",rownames(rna_data),fixed = T)
rna_data$name <- rownames(rna_data)
rna_data2 <- rna_data[!duplicated(rownames(rna_data)),]
rna_gene <-rna_data2[rownames(otu),]
rna_gene <- cbind(rna_gene,otu)

if (!is.null(gene)) {
  x=as.numeric(rna_gene[,taxaa])
  y=as.numeric(rna_gene[,gene])
  obj <- try(cor  <- cor.test(x, y, method=correlation_coef),silent=TRUE)
  if (is(obj, "try-error")){
    q()
  } else{
  }
  cor  <- cor.test(x, y, method=correlation_coef)
  outTab=data.frame()
  outVector=cbind(Xlab=taxaa, Ylab=gene, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  outFile=paste0("result/cancer_type_correlation_analysis_",cancer_type,"_",taxaa,"_cor_", gene,"_",correlation_coef,".pdf")
  df1=as.data.frame(cbind(x,y))
  p1=ggplot(df1, aes(x, y)) + 
    xlab(paste0(taxaa)) + ylab(gene)+
    geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = correlation_coef, aes(x =x, y =y))
  p2=ggMarginal(p1, type="density",  xparams=list(fill = "orange"), 
                yparams=list(fill = "blue"))
  
  #写出list和plot
  ggsave(p2,filename=outFile,width=4,height=4,units=c("in"))
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")
  file.exists("result/Rplots.pdf")
  file.remove("result/Rplots.pdf")
  write.table(outTab,file=paste0("result/temp/",cancer_type,"_",taxaa,"_",gene,"_",correlation_coef,"_cor_","result.txt"),
              sep="\t",row.names=F,quote=F)
} else {
  print("error")
  q()
}
## RUN Example 
# Rscript script_cancer_type_correlation_analysis_microbe_gene.R 'LUSC' 'g__Aspergillus' 'DPM1' 'pearson' '6' '6' 'Fungus'
