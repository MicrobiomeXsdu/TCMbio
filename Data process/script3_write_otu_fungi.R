setwd("/data_200t/shengdahsuang2/TCGA_RNA/CESC/temp/")
library(dplyr)

name <- read.table("ls.txt")

name$V2 <- paste0("s", c(2001:(2000+length(name$V1))))

names(name) <- c("file_id", "kresult_id")

kresult <- read.table(paste0("bracken_report_fungi/",name[1,1],".new.report"))
names(kresult)[1] <- "taxa"
colnames(kresult)[2] <- name[1,2]

#length(name$V1)是样本量的个数
for (i in 2:nrow(name)) {
  a <- read.table(paste0("bracken_report_fungi/",name[i,1],".new.report"))
  names(a)[1] <- "taxa"
  kresult <- merge(kresult, a, all = T)
  colnames(kresult)[i+1] <- name[i,2]
  print(setdiff(row.names(a), row.names(kresult)))
}
kresult <- kresult %>% tibble::column_to_rownames("taxa")
kresult <- kresult[-grep("k__Metazoa",rownames(kresult)),]


kresult[is.na(kresult)] <- 0

write.csv(kresult, "../result/composition_fungi/CESC_kresult.csv", quote = F)

write.csv(name, "../result/composition_fungi/CESC_kresult_name.csv",row.names = F, quote = F)
