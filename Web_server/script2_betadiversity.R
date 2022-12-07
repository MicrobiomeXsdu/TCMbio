library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(phyloseq)
pacman::p_load(tidyverse,ape,ggsignif,patchwork,multcomp)
library(umap)
require("RColorBrewer")
library(tsne)
# 参数读取
args <- commandArgs(T)
group = args[1]
if(group!="MultipleCancersComparison"){
  CancerType = prefix = args[2]
  level = tolower(substr(args[3],1,1))
  visualization.method = args[4]
  Kingdom = args[5]
  plot.height =  args[6]%>% as.numeric()
  plot.width = args[7] %>% as.numeric()
}else{
  CancerTypeA = prefix = args[2]
  CancerTypeB = prefix2 = args[3]
  level = tolower(substr(args[4],1,1))
  visualization.method = args[5]
  Kingdom = args[6]
  plot.height =  args[7]%>% as.numeric()
  plot.width =  args[8]%>% as.numeric()
}


# 参数示例
# group = "clinical_2_stage"
# if(group!="MultipleCancersComparison"){
#   CancerType = prefix = "CESC"
#   level = tolower(substr("Species",1,1))
#   visualization.method ="NMDS"
#   Kingdom = "Fungus"
#   plot.height =  "6" %>% as.numeric()
#   plot.width = "6" %>% as.numeric()
# }else{
#   CancerTypeA = prefix = "BRCA"
#   CancerTypeB = prefix2 = "CESC"
#   level = tolower(substr("Species",1,1))
#   visualization.method = "Species"
#   Kingdom = "Fungus"
#   plot.height = "6"  %>% as.numeric()
#   plot.width =  "6" %>% as.numeric()
# }

datadir <- '/data/jinchuandi/'
# datadir <- './'
# 创建结果目录
if (!dir.exists('result')) {
  dir.create('result')
}

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

if(visualization.method=="PCoA"){
if(max(length(unique(phy_dat@sam_data[[group]]))<2|is.na(unique(phy_dat@sam_data[[group]])))==1){
  stop("The number of groups is less than two")
}else{

phy_dat <- prune_taxa(taxa_sums(phy_dat)>0, phy_dat)

if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
data.h=decostand(otu,method="hellinger")
groups <- phy_dat@sam_data %>%
  as.list()


#bray距离
bray_dist <- vegdist(data.h,method = "bray")
  pcoa <- pcoa(bray_dist, correction = "none", rn = NULL)
  PC1 = pcoa$vectors[,1]
  PC2 = pcoa$vectors[,2]
  pcoadata <- data.frame(rownames(pcoa$vectors),
                         PC1,PC2,groups[group])
  colnames(pcoadata) <-c("sample","PC1","PC2","group")
  
  pcoadata$group <- paste0(toupper(substr(pcoadata$group, 1,1)), substr(pcoadata$group, 2,nchar(pcoadata$group)))
  pcoadata$group <- factor(pcoadata$group)
  
  str(pcoadata)
  yf <- pcoadata
  yd1 <- yf %>% group_by(group) %>% summarise(Max = max(PC1))
  yd2 <- yf %>% group_by(group) %>% summarise(Max = max(PC2))
  yd1$Max <- yd1$Max + max(yd1$Max)*0.1
  yd2$Max <- yd2$Max + max(yd2$Max)*0.1
  
  res1 <- aov(PC1~group,data = pcoadata) %>% 
    glht(linfct=mcp(group="Tukey")) %>% cld(alpah=0.05)
  res2 <- aov(PC2~group,data = pcoadata) %>% 
    glht(linfct=mcp(group="Tukey")) %>% cld(alpah=0.05)
  
  test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                     yd1 = yd1$Max,yd2 = yd2$Max,group = yd1$group)
  test$group <- as.factor(test$group)
  p1 <- ggplot(pcoadata, aes(PC1, PC2)) +
    geom_point(aes(colour=group),size=2.5)+
    stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+
    labs(x=(floor(pcoa$values$Relative_eig[1]*100)) %>% 
           paste0("PCo1 ( ", ., "%", " )"),
         y=(floor(pcoa$values$Relative_eig[2]*100)) %>% 
           paste0("PCo2 ( ", ., "%", " )")) +
    theme(text=element_text(size=12))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          axis.title.x=element_text(colour='black', size=12),
          axis.title.y=element_text(colour='black', size=12),
          axis.text=element_text(colour='black',size=12),
          legend.title=element_blank(),
          legend.key.height=unit(0.6,"cm"),
          legend.position = c(0.75, 0.95),legend.direction = "horizontal")
  p2 <- ggplot(pcoadata,aes(group,PC1)) +
    geom_boxplot(aes(fill = group))+
    # scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))+
    geom_jitter(shape=16,size=1.5,position=position_jitter(0.2))+
    geom_text(data = test,aes(x = group,y = yd1,label = PC1),
              size = 5,color = "black",fontface = "plain")+
    theme(panel.background = element_rect(fill='white',
                                          colour='black'))+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(colour='black',size=12,face = "plain"),
          axis.text.x=element_blank(),
          legend.position = "none")+coord_flip()
  p3 <- ggplot(pcoadata,aes(group,PC2)) +
    geom_boxplot(aes(fill = group)) +
    # scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))+
    geom_jitter(shape=16,size=1.5,position=position_jitter(0.2))+
    geom_text(data = test,aes(x = group,y = yd2,label = PC2),
              size = 5,color = "black",fontface = "plain")+
    theme(panel.background = element_rect(fill='white',
                                          colour='black'))+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=12,angle = 35,
                                   vjust = 1,hjust = 0.5,face = "plain"),
          axis.text.y=element_blank(),
          legend.position = "none")
  
  set.seed(1234)
  otu.adonis=adonis2(bray_dist~group,data = pcoadata)
  p4 <- ggplot(pcoadata,
               aes(PC1, PC2))+
    geom_text(aes(x = -0.5,
                  y = 0.6,
                  label = paste("PERMANOVA:\ndf = ",
                                otu.adonis$Df[1],"\nR2 = ",
                                round(otu.adonis$R2[1],4),
                                "\np-value = ",
                                otu.adonis$`Pr(>F)`[1],
                                sep = "")),size = 4) +theme_bw() +
    xlab(NULL) + ylab(NULL) +
    theme(panel.grid=element_blank(), 
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  betaplot <- p2+p4+p1+p3 + 
          plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  if(group=="MultipleCancersComparison"){
  ggsave(betaplot,filename=paste0("result/",prefix, "_", prefix2, "_bray_dist_", visualization.method, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"),limitsize = FALSE)
}else{
  ggsave(betaplot,filename=paste0("result/",prefix,"_bray_dist_",visualization.method,"_", level, "_", group, ".pdf"),width=plot.width, height=plot.height,  units=c("in"),limitsize = FALSE)
}
}
}

# NMDS
if(visualization.method=="NMDS"){
  if(max(length(unique(phy_dat@sam_data[[group]]))<2|is.na(unique(phy_dat@sam_data[[group]])))==1){
  stop("The number of groups is less than two")
}else{
  
phy_dat <- prune_taxa(taxa_sums(phy_dat)>0, phy_dat)

if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
data.h=decostand(otu,method="hellinger")
groups <- phy_dat@sam_data %>%
  as.list()

  data.h.nmds <- metaMDS(data.h,distance="bray")
  # stressplot(data.h.nmds, main = "Shepard 图")
  map<- scores(data.h.nmds)#提取前两轴坐标/data.h.nmds$points#查看前两轴物种坐标
  #ggplot2绘图
  library(ggplot2)
  #我们将nmds_dis_site数据导出来，然后加上分组信息再读入
  # map<-data.frame(map)
  map<-data.frame(map$sites)
  #也可以将excel表中的分组信息(需要导入R中)与前两轴坐标进行合并(map = merge(group, map,by="row.names",all=F))
  map$name <- rownames(map)
  map$group <- phy_dat@sam_data[[group]]
  head(map)
  #正式绘图
  betaplot <- ggplot(data = map, aes(NMDS1, NMDS2))+
    geom_point(size = 2,#绘制散点图
               aes(color = group, shape = group))+
    stat_ellipse(aes(fill = group),#加上置信椭圆
                 geom = 'polygon',
                 level = 0.95, #置信区间95%
                 alpha = 0.1, #透明度
                 segments = 51, #椭圆类型
                 show.legend = F)+
    theme_bw()+#主题
    geom_vline(xintercept = 0, color = 'black', size = 1, linetype = 3)+ #添加过原点线(y)
    geom_hline(yintercept = 0, color = 'black', size = 1, linetype = 3)+ #添加过原点线(x)
    xlab("NMDS1")+#x轴标题
    ylab("NMDS2")+#y轴标题
    theme(axis.title.y = element_text(size = 14))+ #y轴标题大小
    theme(axis.title.x = element_text(size = 14))+ #x轴标题大小
    geom_text(aes(x=-0.8,y=0.55,label=paste0("stress = ", round(data.h.nmds$stress,4))),size = 6)+ #设置文字位置、内容和大小
    # geom_text(aes(x=-0.8,y=0.30,label=paste0("R^2 = ", data.h.nmds$r)),size = 6)+ #添加R^2
    theme(axis.text.x=element_text(size=14,angle=0,color="Black"),#设置x和y轴字体大小和颜色
          axis.text.y=element_text(size=14,angle=0,color="Black"))+
    # theme(legend.position=c(0.95,0.9))+#设置图例位置
    theme(legend.text=element_text(size=14))+ #设置图例字体大小
    theme(legend.title = element_text(face = "bold", size = 14))#设置图例标题字和大小
if(group=="MultipleCancersComparison"){
  ggsave(betaplot,filename=paste0("result/",prefix, "_", prefix2, "_bray_dist_", visualization.method, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"),limitsize = FALSE)
}else{
  ggsave(betaplot,filename=paste0("result/",prefix,"_bray_dist_",visualization.method,"_", level, "_", group, ".pdf"),width=plot.width, height=plot.height,  units=c("in"),limitsize = FALSE)
}
}
}
 
# UMAP
if(visualization.method=="UMAP"){
  if(max(length(unique(phy_dat@sam_data[[group]]))<2|is.na(unique(phy_dat@sam_data[[group]])))==1){
  print("The number of groups is less than two")
}else{
  
phy_dat <- prune_taxa(taxa_sums(phy_dat)>0, phy_dat)

if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
data.h=decostand(otu,method="hellinger")
groups <- phy_dat@sam_data %>%
  as.list()

# 使用umap函数进行UMAP降维分析
umap = umap(otu)
umap_obj <- as.data.frame(umap$layout)

names(umap_obj) <- c("UMAP_1", "UMAP_2")
umap_obj$group <- phy_dat@sam_data[[group]]
umap_obj$group <- paste0(toupper(substr(umap_obj$group, 1,1)), substr(umap_obj$group, 2,nchar(umap_obj$group)))
# 可视化UMAP的结果
betaplot <- ggplot(umap_obj, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour=group),size=2.5)+
    labs(x="UMAP_1",
         y="UMAP_2",
         title = visualization.method) +
    theme(text=element_text(size=12))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          axis.title.x=element_text(colour='black', size=12),
          axis.title.y=element_text(colour='black', size=12),
          axis.text=element_text(colour='black',size=12),
          legend.title=element_blank(),
          legend.key.height=unit(0.6,"cm"),
          legend.position = c(0.55, 0.95),legend.direction = "horizontal")
if(group=="MultipleCancersComparison"){
  ggsave(betaplot,filename=paste0("result/",prefix, "_", prefix2, "_bray_dist_", visualization.method, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"),limitsize = FALSE)
}else{
  ggsave(betaplot,filename=paste0("result/",prefix,"_bray_dist_",visualization.method,"_", level, "_", group, ".pdf"),width=plot.width, height=plot.height,  units=c("in"),limitsize = FALSE)
}
}
}

# TSNE

if(visualization.method=="T-SNE"){
if(max(length(unique(phy_dat@sam_data[[group]]))<2|is.na(unique(phy_dat@sam_data[[group]])))==1){
  print("The number of groups is less than two")
}else{
  
phy_dat <- prune_taxa(taxa_sums(phy_dat)>0, phy_dat)

if(phy_dat@otu_table@taxa_are_rows==TRUE){
  otu <- phy_dat@otu_table@.Data %>%
    t()
}else{
  otu <- phy_dat@otu_table@.Data
}
data.h=decostand(otu,method="hellinger")
groups <- phy_dat@sam_data %>%
  as.list()


colors = rainbow(length(unique(iris$Species)))
names(colors) = unique(iris$Species)
head(colors)
##      setosa  versicolor   virginica 
## "#FF0000FF" "#00FF00FF" "#0000FFFF"

# 使用tsne函数进行tSNE降维分析
tsne = tsne(otu,k=2,perplexity=50)

tsne_obj <- as.data.frame(tsne)
names(tsne_obj) <- c("TSNE_1", "TSNE_2")
tsne_obj$group <- phy_dat@sam_data[[group]]
tsne_obj$group <- paste0(toupper(substr(tsne_obj$group, 1,1)), substr(tsne_obj$group, 2,nchar(tsne_obj$group)))

# 可视化UMAP的结果
betaplot <- ggplot(tsne_obj, aes(TSNE_1, TSNE_2)) +
    geom_point(aes(colour=group),size=2.5)+
    labs(x="TSNE_1",
         y="TSNE_2",
         title = visualization.method) +
    theme(text=element_text(size=12))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          axis.title.x=element_text(colour='black', size=12),
          axis.title.y=element_text(colour='black', size=12),
          axis.text=element_text(colour='black',size=12),
          legend.title=element_blank(),
          legend.key.height=unit(0.6,"cm"),
          legend.direction = "horizontal")
if(group=="MultipleCancersComparison"){
  ggsave(betaplot,filename=paste0("result/",prefix, "_", prefix2, "_bray_dist_", visualization.method, "_", level, "_", group, ".pdf"),width=plot.width,height=plot.height,units=c("in"),limitsize = FALSE)
}else{
  ggsave(betaplot,filename=paste0("result/",prefix,"_bray_dist_",visualization.method,"_", level, "_", group, ".pdf"),width=plot.width, height=plot.height,  units=c("in"),limitsize = FALSE)
}
}
}

# RUN example
# Rscript script2_betadiversity.R clinical_2_stage CESC Species NMDS Fungus 8 8
