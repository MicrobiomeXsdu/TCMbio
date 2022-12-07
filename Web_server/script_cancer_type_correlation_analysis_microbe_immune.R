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
signatures <- args[3]
correlation_coef <- args[4]
plot.height <- args[5] %>% as.numeric()
plot.width <- args[6] %>% as.numeric()
Kingdom <- args[7]

# demo
# setwd("D:/anlysis-data-r/pancacer_rna/rdata_decontaminate/10_24/")
datadir <- '/data/jinchuandi/'

# cancer_type = "LUSC"
# taxaa <- "g__Yarrowia"
# signatures <- "T.cells.CD4.naive"
# correlation_coef <- "pearson"
# plot.height=6
# plot.width=6
# Kingdom <- "Fungus"



if (Kingdom == "Fungus") {
  filename <- "rawdata_fungi/"
} else {
  filename <- "rawdata/"
}

filename <- paste0(datadir, filename)
phydata_g = paste0(filename, cancer_type, "_g", "_all_logCPM.RData")
phydata_s = paste0(filename, cancer_type, "_s", "_all_logCPM.RData")
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
cibersort = read.csv(paste0(datadir, "rawdata/TCGA-cibersort-relative.csv"),header=T,check.names=F)
cibersort_sub <- cibersort[cibersort$type==cancer_type,]
rm(cibersort)

cibersort_sub$name <- substr(cibersort_sub$SampleID,1,15)
cibersort_sub$name2 <- substr(cibersort_sub$SampleID,1,16)

if (is.null(signatures)) {
  print("Not input Signatures")
} else {
  if (signatures %in% colnames(cibersort_sub)[3:24]){
    print("The entered signature is recognized")
  } else {
    print("The entered signature is not recognized")
    q()
  }
}

otu_g$name <- rownames(otu_g)
otu_s$name <- rownames(otu_s)
otu <- inner_join(otu_g,otu_s)
rownames(otu) <- substr(otu$name,1,15)

micro_immune <- inner_join(cibersort_sub,otu)
micro_immune <- micro_immune[!duplicated(micro_immune$name2),]
rownames(micro_immune) <- micro_immune$name2

signatures_a <- gsub("."," ",signatures,fixed = T)

if (!is.null(signatures)){
  x=as.numeric(micro_immune[,taxaa])
  y=as.numeric(micro_immune[,signatures])
  obj <- try(cor  <- cor.test(x, y, method=correlation_coef),silent=TRUE)
  if (is(obj, "try-error")){
    q()
  } else{
  }
  cor  <- cor.test(x, y, method=correlation_coef)
  outTab=data.frame()
  outVector=cbind(Xlab=taxaa, Ylab=signatures, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  outFile=paste0("result/cancer_type_correlation_analysis_",cancer_type,"_",taxaa,"_cor_", signatures,"_",correlation_coef,".pdf")
  df1=as.data.frame(cbind(x,y))
  p1=ggplot(df1, aes(x, y)) + 
    xlab(paste0(taxaa)) + ylab(signatures_a)+
    geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = correlation_coef, aes(x =x, y =y))
  p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), 
                yparams=list(fill = "blue"))
  
  #写出list和plot
  ggsave(p2,filename=outFile,width=4,height=4,units=c("in"))
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")
  file.exists("result/Rplots.pdf")
  file.remove("result/Rplots.pdf")
  write.table(outTab,file=paste0("result/temp/",cancer_type,"_",taxaa,"_",signatures,"_",correlation_coef,"_cor_","result.txt"),sep="\t",row.names=F,quote=F)
} else { 
  print("error")
  q()
}
# RUN Example
# Rscript script_cancer_type_correlation_analysis_microbe_immune.R 'LUSC' 'g__Yarrowia' 'T.cells.CD4.naive' 'pearson' 6 6 'Fungus'